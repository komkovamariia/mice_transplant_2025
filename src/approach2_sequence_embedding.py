"""Approach 2: sequence embedding, local density enrichment, and repertoire geometry."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from mir.embedding.tcremp import TCREmp
from mir.embedding.presets import get_preset
from mir.density import DensitySpace, calibrate_radius, neighbor_enrichment, _embed
from mir.repertoire import _make_rff, SampleEmbedding, mmd_matrix, _hill

from .figures import save_figure
from .runtime import available_cpus, permanova_parallel
from .strata import STRATUM_LABELS, ensure_groups, repository_root, select_stratum

APPROACH = "02_sequence_embedding"


def _output_dir() -> Path:
    path = repository_root() / "outputs" / "tables" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _cache_dir() -> Path:
    path = repository_root() / "outputs" / "cache" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _embed_frame(df: pd.DataFrame) -> pl.DataFrame:
    if "v_germ" in df.columns:
        v_call = df["v_germ"].where(df["v_germ"].notna(), df["v_gene"]).astype(str)
    else:
        v_call = df["v_gene"].astype(str)
    return pl.DataFrame({"v_call": v_call.to_list(), "j_call": ["TRAJ0*01"] * len(df),
                         "junction_aa": df["cdr3"].astype(str).to_list()})


def prepare_embedding(repertoire: pd.DataFrame, *, seed: int = 0) -> dict:
    """Fit one global coordinate basis; all biological strata are tested in this same space."""
    n_workers = available_cpus()
    preset = get_preset("mouse", "TRA")
    model = TCREmp.from_defaults("mouse", "TRA", mode="cdr123", threads=n_workers)
    cols = ["ckey", "cdr3", "v_gene"] + (["v_germ"] if "v_germ" in repertoire.columns else [])
    uni = repertoire[cols].drop_duplicates("ckey").reset_index(drop=True)

    cache = _cache_dir(); coords_path = cache / "union_coords.npy"; meta_path = cache / "union_meta.parquet"
    basis_path = cache / "embedding_basis.joblib"
    if coords_path.exists() and meta_path.exists() and basis_path.exists():
        from joblib import load
        cached = pd.read_parquet(meta_path)
        if len(cached) == len(uni) and cached["ckey"].equals(uni["ckey"]):
            basis = load(basis_path); coords = np.load(coords_path)
            space = DensitySpace(model=model, space="full", scaler=basis["scaler"], pca=basis["pca"])
            return {"repertoire": repertoire, "uni": uni, "coords": coords, "model": model,
                    "preset": preset, "space": space, "scaler": basis["scaler"], "pca": basis["pca"]}

    rng = np.random.default_rng(seed); fit_n = min(60_000, len(uni)); fit_idx = rng.choice(len(uni), fit_n, replace=False)
    xfit = _embed(model, _embed_frame(uni.iloc[fit_idx]), "full").astype(np.float32)
    scaler = StandardScaler().fit(xfit)
    n_components = min(int(preset.n_components), xfit.shape[0] - 1, xfit.shape[1])
    pca = PCA(n_components=n_components, random_state=seed).fit(scaler.transform(xfit))
    coords = np.empty((len(uni), n_components), dtype=np.float32)
    chunk = 50_000
    for start in range(0, len(uni), chunk):
        raw = _embed(model, _embed_frame(uni.iloc[start:start + chunk]), "full").astype(np.float32)
        coords[start:start + len(raw)] = pca.transform(scaler.transform(raw)).astype(np.float32)

    from joblib import dump
    np.save(coords_path, coords); uni.to_parquet(meta_path, index=False); dump({"scaler": scaler, "pca": pca}, basis_path)
    space = DensitySpace(model=model, space="full", scaler=scaler, pca=pca)
    return {"repertoire": repertoire, "uni": uni, "coords": coords, "model": model,
            "preset": preset, "space": space, "scaler": scaler, "pca": pca}


def _group_umi(sub: pd.DataFrame) -> pd.DataFrame:
    wide = sub.groupby(["ckey", "group"], as_index=False)["umi"].sum().pivot(index="ckey", columns="group", values="umi").fillna(0)
    for g in ("g1", "g5", "g6"):
        if g not in wide: wide[g] = 0.0
    return wide[["g1", "g5", "g6"]]


def _plot_density(ranking: pd.DataFrame, enriched: pd.DataFrame, stratum: str) -> None:
    top = ranking.head(15).sort_values("n_enriched")
    if not top.empty:
        fig, ax = plt.subplots(figsize=(8, 6)); ax.barh(top["v_gene"], top["n_enriched"])
        ax.set_xlabel("Locally enriched aaV clonotypes"); ax.set_ylabel("V segment")
        ax.set_title(f"Embedding-density ranking: {STRATUM_LABELS[stratum]}")
        save_figure(fig, APPROACH, stratum, "density_v_segment_ranking")
    if not enriched.empty:
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        ax.scatter(np.log2(enriched["fold"].clip(lower=1e-12)),
                   -np.log10(enriched["qvalue"].clip(lower=np.finfo(float).tiny)), s=10, alpha=0.55)
        ax.set_xlabel("log2 local enrichment"); ax.set_ylabel("-log10 q-value")
        ax.set_title(f"Local sequence-space enrichment: {STRATUM_LABELS[stratum]}")
        save_figure(fig, APPROACH, stratum, "density_enrichment_scatter")


def _classical_mds(distance: np.ndarray) -> np.ndarray:
    n = len(distance); h = np.eye(n) - np.ones((n, n)) / n; b = -0.5 * h @ (distance ** 2) @ h
    vals, vecs = np.linalg.eigh(b); order = np.argsort(vals)[::-1][:2]; vals = np.maximum(vals[order], 0)
    return vecs[:, order] * np.sqrt(vals)


def _sample_geometry(context: dict, sub: pd.DataFrame, stratum: str):
    uni = context["uni"]; coords = context["coords"]; key_to_row = dict(zip(uni["ckey"], np.arange(len(uni))))
    length_scale = calibrate_radius(context["space"], seed=0); rff = _make_rff(coords.shape[1], 2048, length_scale, 0)
    persamp = sub.groupby(["sample_id", "ckey"], as_index=False)["umi"].sum()
    smeta = sub.groupby("sample_id").agg(group=("group", "first"), mouse_id=("mouse_id", "first")).reset_index()
    means, neffs, divs, kept = [], [], [], []
    for sid, block in persamp.groupby("sample_id"):
        block = block[block["ckey"].isin(key_to_row)].copy(); rows = np.array([key_to_row[k] for k in block["ckey"]], dtype=int)
        if len(rows) < 2: continue
        abundance = block["umi"].to_numpy(float); weights = np.log1p(abundance)
        if weights.sum() <= 0: continue
        weights /= weights.sum(); z = coords[rows].astype(np.float64)
        psi = np.sqrt(2.0 / 2048) * np.cos(z @ rff.omega + rff.b)
        means.append((weights @ psi).astype(np.float32)); neffs.append(float(1.0 / np.sum(weights * weights)))
        freq = abundance / abundance.sum(); d0, d1, d2 = _hill(freq); divs.append([np.log(d0), np.log(d1), np.log(d2)]); kept.append(sid)
    means = np.asarray(means); neffs = np.asarray(neffs); divs = np.asarray(divs)
    meta = smeta.set_index("sample_id").loc[kept].reset_index()
    embs = [SampleEmbedding(mean=means[i].astype(np.float64), diversity=divs[i], second=None, n_eff=float(neffs[i])) for i in range(len(kept))]
    d2 = np.maximum(mmd_matrix(embs, unbiased=True), 0); distance = np.sqrt(d2)
    perma_F, perma_p = permanova_parallel(distance, meta["group"].to_numpy(), n_perm=9999, seed=0)
    perma = {"F": float(perma_F), "p": float(perma_p), "n_perm": 9999}
    xy = _classical_mds(distance); fig, ax = plt.subplots(figsize=(7, 6))
    for group in sorted(meta["group"].unique()):
        mask = meta["group"].eq(group).to_numpy(); ax.scatter(xy[mask, 0], xy[mask, 1], label=group, s=35)
    ax.set_xlabel("PCoA 1"); ax.set_ylabel("PCoA 2"); ax.set_title(f"RFF-MMD repertoire geometry: {STRATUM_LABELS[stratum]}")
    ax.legend(frameon=False); save_figure(fig, APPROACH, stratum, "repertoire_mmd_pcoa")
    meta.to_csv(_output_dir() / f"{stratum}_sample_metadata.csv", index=False); np.save(_output_dir() / f"{stratum}_mmd_matrix.npy", distance)
    return meta, perma, means, rff


def _witness(context: dict, sub: pd.DataFrame, stratum: str, meta: pd.DataFrame, means: np.ndarray, rff) -> pd.DataFrame:
    groups = meta["group"].to_numpy(); g1 = means[groups == "g1"]; allo = means[np.isin(groups, ["g5", "g6"])]
    if len(g1) == 0 or len(allo) == 0: return pd.DataFrame()
    witness = g1.mean(0) - allo.mean(0); uni = context["uni"]; coords = context["coords"]; membership = _group_umi(sub)
    keys = membership.index[membership["g1"].gt(0)]; key_to_row = dict(zip(uni["ckey"], np.arange(len(uni))))
    keys = [k for k in keys if k in key_to_row]; rows = np.array([key_to_row[k] for k in keys], dtype=int)
    scores = np.empty(len(rows), dtype=float); chunk = 100_000
    for start in range(0, len(rows), chunk):
        z = coords[rows[start:start + chunk]].astype(np.float64); psi = np.sqrt(2.0 / 2048) * np.cos(z @ rff.omega + rff.b)
        scores[start:start + len(z)] = psi @ witness
    lookup = uni.set_index("ckey"); result = lookup.loc[keys, ["cdr3", "v_gene"]].reset_index(); result["witness_score"] = scores
    result = result.sort_values("witness_score", ascending=False).reset_index(drop=True); result["stratum"] = stratum
    result.to_csv(_output_dir() / f"{stratum}_witness_clonotypes.csv", index=False)
    top = result.head(500)["v_gene"].value_counts().head(15).sort_values()
    if not top.empty:
        fig, ax = plt.subplots(figsize=(8, 6)); ax.barh(top.index, top.values)
        ax.set_xlabel("Clonotypes among top 500 witness scores"); ax.set_ylabel("V segment")
        ax.set_title(f"Witness V-segment representation: {STRATUM_LABELS[stratum]}")
        save_figure(fig, APPROACH, stratum, "witness_v_segment_ranking")
    return result


def _motif_summary(enriched: pd.DataFrame, stratum: str) -> pd.DataFrame:
    rows = []
    for v_gene, block in enriched.groupby("v_gene"):
        seqs = block["cdr3"].dropna().astype(str)
        if len(seqs) < 3: continue
        modal_len = int(seqs.str.len().mode().iloc[0]); seqs = seqs[seqs.str.len().eq(modal_len)]
        if len(seqs) < 3: continue
        ent, consensus = [], []
        for pos in range(modal_len):
            f = seqs.str[pos].value_counts(normalize=True); p = f.to_numpy(float)
            ent.append(float(-(p * np.log2(p)).sum())); consensus.append(str(f.index[0]))
        rows.append({"stratum": stratum, "v_gene": v_gene, "n_sequences": len(seqs), "modal_length": modal_len,
                     "mean_position_entropy": float(np.mean(ent)), "consensus_cdr3": "".join(consensus)})
    out = pd.DataFrame(rows)
    if not out.empty: out = out.sort_values(["mean_position_entropy", "n_sequences"], ascending=[True, False])
    out.to_csv(_output_dir() / f"{stratum}_motif_summary.csv", index=False); return out


def run_stratum(context: dict, stratum: str) -> dict:
    sub = select_stratum(context["repertoire"], stratum); ensure_groups(sub, ("g1", "g5", "g6"), f"Approach 2 / {stratum}")
    membership = _group_umi(sub); uni = context["uni"]; pos = uni.set_index("ckey")
    common = membership.index.intersection(pos.index); membership = membership.loc[common]; rows = pos.index.get_indexer(common); coords = context["coords"][rows]
    obs_mask = membership["g1"].to_numpy() > 0; bg_mask = (membership["g5"].to_numpy() > 0) | (membership["g6"].to_numpy() > 0)
    if obs_mask.sum() < 2 or bg_mask.sum() < 2: raise ValueError(f"Approach 2 / {stratum}: insufficient g1 or allogeneic clonotypes.")
    result = neighbor_enrichment(coords[obs_mask], coords[bg_mask], radius=None, lambda0=3.0, test="poisson", calibrate="median",
                                 abundance=membership.loc[common[obs_mask], "g1"].to_numpy(float), weight="log1p", orphan=True, backend="kdtree")
    obs_keys = common[obs_mask]; obs = pos.loc[obs_keys, ["cdr3", "v_gene"]].reset_index(); obs["g1_umi"] = membership.loc[obs_keys, "g1"].to_numpy(float)
    obs["fold"] = result.fold; obs["qvalue"] = result.qvalue; obs["density_score"] = result.score; obs["stratum"] = stratum
    enriched = obs[(obs["qvalue"] < 0.05) & (obs["fold"] > 1)].copy()
    ranking = enriched.groupby("v_gene", as_index=False).agg(n_enriched=("ckey", "size"), median_fold=("fold", "median"), max_fold=("fold", "max"), total_g1_umi=("g1_umi", "sum"))
    totals = obs.groupby("v_gene").size().rename("n_total").reset_index(); ranking = ranking.merge(totals, on="v_gene", how="outer").fillna({"n_enriched": 0})
    ranking["enrichment_rate"] = ranking["n_enriched"] / ranking["n_total"].replace(0, np.nan)
    ranking = ranking.sort_values(["n_enriched", "median_fold", "total_g1_umi"], ascending=[False, False, False], na_position="last").reset_index(drop=True)
    ranking["rank"] = np.arange(1, len(ranking) + 1); ranking["score"] = ranking["n_enriched"].astype(float); ranking["stratum"] = stratum
    out = _output_dir(); obs.to_csv(out / f"{stratum}_density_clonotypes.csv", index=False); enriched.to_csv(out / f"{stratum}_density_enriched_clonotypes.csv", index=False)
    ranking.to_csv(out / f"{stratum}_v_gene_ranking.csv", index=False); _plot_density(ranking, enriched, stratum)
    meta, perma, means, rff = _sample_geometry(context, sub, stratum); witness = _witness(context, sub, stratum, meta, means, rff); motifs = _motif_summary(enriched, stratum)
    conclusion = {"stratum": stratum, "label": STRATUM_LABELS[stratum], "n_samples": int(sub["sample_id"].nunique()),
                  "n_density_enriched_clonotypes": int(len(enriched)), "top_v_genes": ranking.head(5)["v_gene"].dropna().tolist(),
                  "permanova_F": float(perma["F"]), "permanova_p": float(perma["p"]), "permanova_permutations": int(perma["n_perm"]),
                  "top_witness_v_genes": witness.head(500)["v_gene"].value_counts().head(5).index.tolist() if not witness.empty else [],
                  "lowest_entropy_motifs": motifs.head(5)["consensus_cdr3"].tolist() if not motifs.empty else []}
    (out / f"{stratum}_conclusion.json").write_text(json.dumps(conclusion, indent=2)); return conclusion


def compile_summary(strata) -> pd.DataFrame:
    rows = []; out = _output_dir()
    for stratum in strata:
        path = out / f"{stratum}_conclusion.json"
        if path.exists(): rows.append(json.loads(path.read_text()))
    summary = pd.DataFrame(rows); summary.to_csv(out / "stratum_summary.csv", index=False); return summary
