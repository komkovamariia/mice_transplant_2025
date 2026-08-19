"""Approach 2: sequence-space enrichment and repertoire geometry."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import polars as pl
from mir.density import (
    DensitySpace,
    _embed,
    calibrate_radius,
    neighbor_enrichment,
)
from mir.embedding.presets import get_preset
from mir.embedding.tcremp import TCREmp
from mir.repertoire import (
    SampleEmbedding,
    _hill,
    _make_rff,
    mmd_matrix,
)
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .figures import save_figure
from .reporting import format_gene_list, save_conclusion
from .runtime import available_cpus, permanova_parallel
from .strata import (
    STRATUM_LABELS,
    analysis_metadata,
    ensure_groups,
    repository_root,
    select_stratum,
)

APPROACH = "02_sequence_embedding"


def _output_dir() -> Path:
    path = repository_root() / "results" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _cache_dir() -> Path:
    path = repository_root() / "results" / "cache" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _embed_frame(df: pd.DataFrame) -> pl.DataFrame:
    if "v_germ" in df.columns:
        v_call = df["v_germ"].where(df["v_germ"].notna(), df["v_gene"]).astype(str)
    else:
        v_call = df["v_gene"].astype(str)
    return pl.DataFrame(
        {
            "v_call": v_call.to_list(),
            "j_call": ["TRAJ0*01"] * len(df),
            "junction_aa": df["cdr3"].astype(str).to_list(),
        }
    )


def _fingerprint(frame: pd.DataFrame) -> str:
    hashed = pd.util.hash_pandas_object(frame, index=False).values.tobytes()
    return hashlib.sha256(hashed).hexdigest()


def prepare_embedding(
    repertoire: pd.DataFrame,
    *,
    seed: int = 0,
) -> dict:
    """Fit one technical coordinate basis shared by all stratum-specific tests."""
    n_workers = available_cpus()
    preset = get_preset("mouse", "TRA")
    model = TCREmp.from_defaults(
        "mouse",
        "TRA",
        mode="cdr123",
        threads=n_workers,
    )

    columns = ["ckey", "cdr3", "v_gene"]
    if "v_germ" in repertoire.columns:
        columns.append("v_germ")
    unique = (
        repertoire[columns]
        .drop_duplicates("ckey")
        .sort_values("ckey")
        .reset_index(drop=True)
    )
    if len(unique) < 3:
        raise ValueError("At least three unique clonotypes are required for embedding.")

    cache = _cache_dir()
    coordinates_path = cache / "union_coordinates.npy"
    metadata_path = cache / "union_metadata.parquet"
    basis_path = cache / "embedding_basis.joblib"
    manifest_path = cache / "embedding_manifest.json"
    fingerprint = _fingerprint(unique)

    cache_files = [
        coordinates_path,
        metadata_path,
        basis_path,
        manifest_path,
    ]
    if all(path.exists() for path in cache_files):
        from joblib import load

        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        cached = pd.read_parquet(metadata_path)
        if (
            manifest.get("fingerprint") == fingerprint
            and manifest.get("seed") == seed
            and cached["ckey"].equals(unique["ckey"])
        ):
            basis = load(basis_path)
            coordinates = np.load(coordinates_path)
            space = DensitySpace(
                model=model,
                space="full",
                scaler=basis["scaler"],
                pca=basis["pca"],
            )
            return {
                "repertoire": repertoire,
                "unique": unique,
                "coordinates": coordinates,
                "model": model,
                "preset": preset,
                "space": space,
                "scaler": basis["scaler"],
                "pca": basis["pca"],
                "fingerprint": fingerprint,
            }

    rng = np.random.default_rng(seed)
    fit_size = min(60_000, len(unique))
    fit_indices = rng.choice(len(unique), fit_size, replace=False)
    fit_embedding = _embed(
        model,
        _embed_frame(unique.iloc[fit_indices]),
        "full",
    ).astype(np.float32)
    scaler = StandardScaler().fit(fit_embedding)
    n_components = min(
        int(preset.n_components),
        fit_embedding.shape[0] - 1,
        fit_embedding.shape[1],
    )
    pca = PCA(
        n_components=n_components,
        random_state=seed,
    ).fit(scaler.transform(fit_embedding))

    coordinates = np.empty(
        (len(unique), n_components),
        dtype=np.float32,
    )
    chunk_size = 50_000
    for start in range(0, len(unique), chunk_size):
        block = unique.iloc[start : start + chunk_size]
        raw = _embed(model, _embed_frame(block), "full").astype(np.float32)
        coordinates[start : start + len(raw)] = pca.transform(
            scaler.transform(raw)
        ).astype(np.float32)

    from joblib import dump

    np.save(coordinates_path, coordinates)
    unique.to_parquet(metadata_path, index=False)
    dump({"scaler": scaler, "pca": pca}, basis_path)
    manifest_path.write_text(
        json.dumps(
            {
                "fingerprint": fingerprint,
                "seed": seed,
                "n_clonotypes": len(unique),
                "n_components": n_components,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    space = DensitySpace(
        model=model,
        space="full",
        scaler=scaler,
        pca=pca,
    )
    return {
        "repertoire": repertoire,
        "unique": unique,
        "coordinates": coordinates,
        "model": model,
        "preset": preset,
        "space": space,
        "scaler": scaler,
        "pca": pca,
        "fingerprint": fingerprint,
    }


def _group_umi(data: pd.DataFrame) -> pd.DataFrame:
    wide = (
        data.groupby(["ckey", "group"], as_index=False)["umi"]
        .sum()
        .pivot(index="ckey", columns="group", values="umi")
        .fillna(0)
    )
    for group in ("g1", "g5", "g6"):
        if group not in wide:
            wide[group] = 0.0
    return wide[["g1", "g5", "g6"]]


def _plot_density(
    ranking: pd.DataFrame,
    enriched: pd.DataFrame,
    stratum: str,
) -> None:
    top = ranking.head(15).sort_values("n_enriched")
    if not top.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.barh(top["v_gene"], top["n_enriched"])
        ax.set_xlabel("Locally enriched aaV clonotypes")
        ax.set_ylabel("V segment")
        ax.set_title(f"Sequence-space density ranking: {STRATUM_LABELS[stratum]}")
        save_figure(
            fig,
            APPROACH,
            stratum,
            "density_v_segment_ranking",
        )

    if not enriched.empty:
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        ax.scatter(
            enriched["effect_g1_vs_allogeneic"],
            -np.log10(enriched["qvalue"].clip(lower=np.finfo(float).tiny)),
            s=10,
            alpha=0.55,
        )
        ax.set_xlabel("Local log2 enrichment: g1 / allogeneic")
        ax.set_ylabel("-log10 q-value")
        ax.set_title(f"Local sequence-space enrichment: {STRATUM_LABELS[stratum]}")
        save_figure(
            fig,
            APPROACH,
            stratum,
            "density_enrichment_scatter",
        )


def _classical_mds(distance: np.ndarray) -> np.ndarray:
    n_samples = len(distance)
    centering = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
    gram = -0.5 * centering @ (distance**2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1][:2]
    eigenvalues = np.maximum(eigenvalues[order], 0)
    return eigenvectors[:, order] * np.sqrt(eigenvalues)


def _sample_geometry(
    context: dict,
    data: pd.DataFrame,
    stratum: str,
):
    unique = context["unique"]
    coordinates = context["coordinates"]
    row_by_key = dict(zip(unique["ckey"], np.arange(len(unique))))

    length_scale = calibrate_radius(context["space"], seed=0)
    rff = _make_rff(
        coordinates.shape[1],
        2048,
        length_scale,
        0,
    )
    per_unit = data.groupby(
        ["analysis_unit", "ckey"],
        as_index=False,
    )["umi"].sum()
    metadata = analysis_metadata(data)

    means = []
    effective_sizes = []
    diversities = []
    retained_units = []
    for analysis_unit, block in per_unit.groupby("analysis_unit"):
        block = block[block["ckey"].isin(row_by_key)].copy()
        rows = np.array(
            [row_by_key[key] for key in block["ckey"]],
            dtype=int,
        )
        if len(rows) < 2:
            continue

        abundance = block["umi"].to_numpy(float)
        weights = np.log1p(abundance)
        if weights.sum() <= 0:
            continue
        weights /= weights.sum()

        embedded = coordinates[rows].astype(np.float64)
        features = np.sqrt(2.0 / 2048) * np.cos(embedded @ rff.omega + rff.b)
        means.append((weights @ features).astype(np.float32))
        effective_sizes.append(float(1.0 / np.sum(weights * weights)))

        frequencies = abundance / abundance.sum()
        d0, d1, d2 = _hill(frequencies)
        diversities.append([np.log(d0), np.log(d1), np.log(d2)])
        retained_units.append(analysis_unit)

    if len(retained_units) < 4:
        raise ValueError(
            f"Approach 2 / {stratum}: fewer than four analysis units "
            "contain at least two clonotypes."
        )

    means = np.asarray(means)
    effective_sizes = np.asarray(effective_sizes)
    diversities = np.asarray(diversities)
    metadata = metadata.loc[retained_units].reset_index()
    metadata["contrast_group"] = metadata["group"].replace(
        {"g5": "allogeneic", "g6": "allogeneic"}
    )
    group_counts = metadata["contrast_group"].value_counts()
    if group_counts.get("g1", 0) < 2 or group_counts.get("allogeneic", 0) < 2:
        raise ValueError(
            f"Approach 2 / {stratum}: PERMANOVA requires at least two "
            "independent units in g1 and allogeneic groups."
        )

    embeddings = [
        SampleEmbedding(
            mean=means[index].astype(np.float64),
            diversity=diversities[index],
            second=None,
            n_eff=float(effective_sizes[index]),
        )
        for index in range(len(retained_units))
    ]
    squared_mmd = np.maximum(
        mmd_matrix(embeddings, unbiased=True),
        0,
    )
    distance = np.sqrt(squared_mmd)
    statistic, p_value = permanova_parallel(
        distance,
        metadata["contrast_group"].to_numpy(),
        n_perm=9999,
        seed=0,
    )
    permanova = {
        "statistic": float(statistic),
        "p_value": float(p_value),
        "n_permutations": 9999,
        "contrast": "g1_vs_g5_plus_g6",
    }

    coordinates_2d = _classical_mds(distance)
    fig, ax = plt.subplots(figsize=(7, 6))
    for group in sorted(metadata["group"].unique()):
        mask = metadata["group"].eq(group).to_numpy()
        ax.scatter(
            coordinates_2d[mask, 0],
            coordinates_2d[mask, 1],
            label=group,
            s=35,
        )
    ax.set_xlabel("PCoA 1")
    ax.set_ylabel("PCoA 2")
    ax.set_title(f"RFF-MMD repertoire geometry: {STRATUM_LABELS[stratum]}")
    ax.legend(frameon=False)
    save_figure(
        fig,
        APPROACH,
        stratum,
        "repertoire_mmd_pcoa",
    )

    output = _output_dir()
    metadata.to_csv(
        output / f"{stratum}_analysis_unit_metadata.csv",
        index=False,
    )
    np.save(output / f"{stratum}_mmd_matrix.npy", distance)
    pd.DataFrame([permanova]).to_csv(
        output / f"{stratum}_permanova.csv",
        index=False,
    )
    return metadata, permanova, means, rff


def _witness(
    context: dict,
    data: pd.DataFrame,
    stratum: str,
    metadata: pd.DataFrame,
    means: np.ndarray,
    rff,
) -> pd.DataFrame:
    groups = metadata["contrast_group"].to_numpy()
    g1 = means[groups == "g1"]
    allogeneic = means[groups == "allogeneic"]
    if len(g1) == 0 or len(allogeneic) == 0:
        return pd.DataFrame()

    witness = g1.mean(0) - allogeneic.mean(0)
    unique = context["unique"]
    coordinates = context["coordinates"]
    membership = _group_umi(data)
    keys = membership.index[membership["g1"].gt(0)]
    row_by_key = dict(zip(unique["ckey"], np.arange(len(unique))))
    keys = [key for key in keys if key in row_by_key]
    rows = np.array([row_by_key[key] for key in keys], dtype=int)

    scores = np.empty(len(rows), dtype=float)
    chunk_size = 100_000
    for start in range(0, len(rows), chunk_size):
        embedded = coordinates[rows[start : start + chunk_size]].astype(np.float64)
        features = np.sqrt(2.0 / 2048) * np.cos(embedded @ rff.omega + rff.b)
        scores[start : start + len(embedded)] = features @ witness

    lookup = unique.set_index("ckey")
    result = lookup.loc[keys, ["cdr3", "v_gene"]].reset_index()
    result["witness_score"] = scores
    result["stratum"] = stratum
    result = result.sort_values(
        "witness_score",
        ascending=False,
    ).reset_index(drop=True)
    result.to_csv(
        _output_dir() / f"{stratum}_witness_clonotypes.csv",
        index=False,
    )

    top = result.head(500)["v_gene"].value_counts().head(15).sort_values()
    if not top.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.barh(top.index, top.values)
        ax.set_xlabel("Clonotypes among the top 500 witness scores")
        ax.set_ylabel("V segment")
        ax.set_title(f"Witness V-segment representation: {STRATUM_LABELS[stratum]}")
        save_figure(
            fig,
            APPROACH,
            stratum,
            "witness_v_segment_ranking",
        )
    return result


def _motif_summary(
    enriched: pd.DataFrame,
    stratum: str,
) -> pd.DataFrame:
    rows = []
    for v_gene, block in enriched.groupby("v_gene"):
        sequences = block["cdr3"].dropna().astype(str)
        if len(sequences) < 3:
            continue
        modal_length = int(sequences.str.len().mode().iloc[0])
        sequences = sequences[sequences.str.len().eq(modal_length)]
        if len(sequences) < 3:
            continue

        entropy = []
        consensus = []
        for position in range(modal_length):
            frequencies = sequences.str[position].value_counts(normalize=True)
            probabilities = frequencies.to_numpy(float)
            entropy.append(float(-(probabilities * np.log2(probabilities)).sum()))
            consensus.append(str(frequencies.index[0]))
        rows.append(
            {
                "stratum": stratum,
                "v_gene": v_gene,
                "n_sequences": len(sequences),
                "modal_length": modal_length,
                "mean_position_entropy": float(np.mean(entropy)),
                "consensus_cdr3": "".join(consensus),
            }
        )

    output = pd.DataFrame(rows)
    if not output.empty:
        output = output.sort_values(
            ["mean_position_entropy", "n_sequences"],
            ascending=[True, False],
        )
    output.to_csv(
        _output_dir() / f"{stratum}_motif_summary.csv",
        index=False,
    )
    return output


def run_stratum(context: dict, stratum: str) -> dict:
    data = select_stratum(context["repertoire"], stratum)
    ensure_groups(data, ("g1", "g5", "g6"), f"Approach 2 / {stratum}")

    membership = _group_umi(data)
    unique = context["unique"]
    indexed_unique = unique.set_index("ckey")
    common = membership.index.intersection(indexed_unique.index)
    membership = membership.loc[common]
    rows = indexed_unique.index.get_indexer(common)
    coordinates = context["coordinates"][rows]

    observed_mask = membership["g1"].to_numpy() > 0
    background_mask = (membership["g5"].to_numpy() > 0) | (
        membership["g6"].to_numpy() > 0
    )
    if observed_mask.sum() < 2 or background_mask.sum() < 2:
        raise ValueError(
            f"Approach 2 / {stratum}: insufficient g1 or allogeneic clonotypes."
        )

    density = neighbor_enrichment(
        coordinates[observed_mask],
        coordinates[background_mask],
        radius=None,
        lambda0=3.0,
        test="poisson",
        calibrate="median",
        abundance=membership.loc[
            common[observed_mask],
            "g1",
        ].to_numpy(float),
        weight="log1p",
        orphan=True,
        backend="kdtree",
    )

    observed_keys = common[observed_mask]
    observed = indexed_unique.loc[
        observed_keys,
        ["cdr3", "v_gene"],
    ].reset_index()
    observed["g1_umi"] = membership.loc[
        observed_keys,
        "g1",
    ].to_numpy(float)
    observed["fold"] = density.fold
    observed["qvalue"] = density.qvalue
    observed["density_score"] = density.score
    observed["effect_g1_vs_allogeneic"] = np.log2(
        observed["fold"].clip(lower=np.finfo(float).tiny)
    )
    observed["stratum"] = stratum

    enriched = observed[observed["qvalue"].lt(0.05) & observed["fold"].gt(1)].copy()
    ranking = enriched.groupby("v_gene", as_index=False).agg(
        n_enriched=("ckey", "size"),
        median_effect_g1_vs_allogeneic=(
            "effect_g1_vs_allogeneic",
            "median",
        ),
        max_effect_g1_vs_allogeneic=(
            "effect_g1_vs_allogeneic",
            "max",
        ),
        total_g1_umi=("g1_umi", "sum"),
    )
    totals = observed.groupby("v_gene").size().rename("n_total").reset_index()
    ranking = ranking.merge(
        totals,
        on="v_gene",
        how="outer",
    ).fillna({"n_enriched": 0})
    ranking["enrichment_rate"] = ranking["n_enriched"] / ranking["n_total"].replace(
        0, np.nan
    )
    ranking = ranking.sort_values(
        [
            "n_enriched",
            "median_effect_g1_vs_allogeneic",
            "total_g1_umi",
        ],
        ascending=[False, False, False],
        na_position="last",
    ).reset_index(drop=True)
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    ranking["score"] = ranking["n_enriched"].astype(float)
    ranking["stratum"] = stratum

    output = _output_dir()
    observed.to_csv(
        output / f"{stratum}_density_clonotypes.csv",
        index=False,
    )
    enriched.to_csv(
        output / f"{stratum}_density_enriched_clonotypes.csv",
        index=False,
    )
    ranking.to_csv(
        output / f"{stratum}_v_gene_ranking.csv",
        index=False,
    )
    _plot_density(ranking, enriched, stratum)

    metadata, permanova, means, rff = _sample_geometry(
        context,
        data,
        stratum,
    )
    witness = _witness(
        context,
        data,
        stratum,
        metadata,
        means,
        rff,
    )
    motifs = _motif_summary(enriched, stratum)

    top_genes = ranking[ranking["n_enriched"].gt(0)].head(5)["v_gene"].dropna().tolist()
    top_witness_genes = (
        witness.head(500)["v_gene"].value_counts().head(5).index.tolist()
        if not witness.empty
        else []
    )
    payload = {
        "stratum": stratum,
        "label": STRATUM_LABELS[stratum],
        "n_input_samples": int(data["sample_id"].nunique()),
        "n_analysis_units": int(data["analysis_unit"].nunique()),
        "n_density_enriched_clonotypes": len(enriched),
        "top_v_genes": top_genes,
        "permanova_statistic": permanova["statistic"],
        "permanova_p_value": permanova["p_value"],
        "permanova_permutations": permanova["n_permutations"],
        "top_witness_v_genes": top_witness_genes,
        "lowest_entropy_motifs": (
            motifs.head(5)["consensus_cdr3"].tolist() if not motifs.empty else []
        ),
    }
    return save_conclusion(
        output,
        stratum,
        f"Approach 2 conclusion: {STRATUM_LABELS[stratum]}",
        payload,
        [
            (
                f"Approach 2 analyzed {payload['n_analysis_units']} independent "
                f"units and identified {len(enriched)} locally g1-enriched clonotypes "
                "in local TCR sequence space at q < 0.05."
            ),
            (
                "The two-group RFF-MMD PERMANOVA for g1 versus g5+g6 "
                f"produced pseudo-F = {permanova['statistic']:.3g} and "
                f"permutation p = {permanova['p_value']:.4g}."
            ),
            (f"The density-supported V segments were {format_gene_list(top_genes)}."),
        ],
    )


def compile_summary(strata) -> pd.DataFrame:
    rows = []
    output = _output_dir()
    for stratum in strata:
        path = output / f"{stratum}_conclusion.json"
        if path.exists():
            rows.append(json.loads(path.read_text(encoding="utf-8")))

    summary = pd.DataFrame(rows)
    summary.to_csv(output / "stratum_summary.csv", index=False)
    return summary
