"""Cross-approach comparison of V-segment rankings across biological strata."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from .figures import save_figure
from .strata import STRATUM_LABELS, repository_root

APPROACH = "04_cross_approach"
SOURCES = {
    "set_count": "01_set_count",
    "sequence_embedding": "02_sequence_embedding",
    "clone_alloreactivity": "03_clone_alloreactivity",
}


def _output_dir() -> Path:
    path = repository_root() / "outputs" / "tables" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _ranking_path(folder: str, stratum: str) -> Path:
    return repository_root() / "outputs" / "tables" / folder / f"{stratum}_v_gene_ranking.csv"


def load_rankings(stratum: str) -> dict[str, pd.DataFrame]:
    out = {}
    missing = []
    for name, folder in SOURCES.items():
        path = _ranking_path(folder, stratum)
        if not path.exists():
            missing.append(str(path))
            continue
        df = pd.read_csv(path)
        required = {"v_gene", "rank"}
        if not required.issubset(df.columns):
            raise ValueError(f"{path} must contain {sorted(required)}.")
        out[name] = df[["v_gene", "rank"]].dropna().drop_duplicates("v_gene")
    if missing:
        raise FileNotFoundError(
            "Cross-approach comparison requires all three upstream approaches. Missing: "
            + "; ".join(missing)
        )
    return out


def _pairwise_metrics(rankings: dict[str, pd.DataFrame], stratum: str) -> pd.DataFrame:
    rows = []
    for a, b in combinations(rankings, 2):
        merged = rankings[a].merge(rankings[b], on="v_gene", suffixes=("_a", "_b"))
        rho, p = spearmanr(merged["rank_a"], merged["rank_b"]) if len(merged) >= 3 else (np.nan, np.nan)
        top10_a = set(rankings[a].nsmallest(10, "rank")["v_gene"])
        top10_b = set(rankings[b].nsmallest(10, "rank")["v_gene"])
        rows.append({
            "stratum": stratum,
            "approach_a": a,
            "approach_b": b,
            "n_shared_v_genes": len(merged),
            "spearman_rho": rho,
            "spearman_p": p,
            "top10_overlap": len(top10_a & top10_b),
        })
    return pd.DataFrame(rows)


def _consensus(rankings: dict[str, pd.DataFrame], stratum: str) -> pd.DataFrame:
    all_v = sorted(set().union(*(set(df["v_gene"]) for df in rankings.values())))
    table = pd.DataFrame({"v_gene": all_v})
    for name, df in rankings.items():
        maximum = max(float(df["rank"].max()), 1.0)
        tmp = df.rename(columns={"rank": f"{name}_rank"}).copy()
        tmp[f"{name}_rank_fraction"] = tmp[f"{name}_rank"] / maximum
        table = table.merge(tmp, on="v_gene", how="left")
    frac_cols = [c for c in table if c.endswith("_rank_fraction")]
    table["n_approaches"] = table[frac_cols].notna().sum(axis=1)
    table["mean_rank_fraction"] = table[frac_cols].mean(axis=1, skipna=True)
    table = table.sort_values(["n_approaches", "mean_rank_fraction"], ascending=[False, True]).reset_index(drop=True)
    table["consensus_rank"] = np.arange(1, len(table) + 1)
    table["stratum"] = stratum
    return table


def _plot_consensus(consensus: pd.DataFrame, stratum: str) -> None:
    top = consensus[consensus["n_approaches"].eq(3)].head(15).sort_values("mean_rank_fraction", ascending=False)
    if top.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(top["v_gene"], 1 - top["mean_rank_fraction"])
    ax.set_xlabel("Cross-approach rank concordance")
    ax.set_ylabel("V segment")
    ax.set_title(f"Consensus V-segment ranking: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "consensus_v_segment_ranking")


def run_stratum(stratum: str) -> dict:
    rankings = load_rankings(stratum)
    metrics = _pairwise_metrics(rankings, stratum)
    consensus = _consensus(rankings, stratum)
    out = _output_dir()
    metrics.to_csv(out / f"{stratum}_pairwise_metrics.csv", index=False)
    consensus.to_csv(out / f"{stratum}_consensus_ranking.csv", index=False)
    _plot_consensus(consensus, stratum)

    top = consensus[consensus["n_approaches"].eq(3)].head(5)
    conclusion = {
        "stratum": stratum,
        "label": STRATUM_LABELS[stratum],
        "top_consensus_v_genes": top["v_gene"].tolist(),
        "pairwise_spearman": {
            f"{row.approach_a}_vs_{row.approach_b}": (
                None if pd.isna(row.spearman_rho) else float(row.spearman_rho)
            )
            for row in metrics.itertuples()
        },
        "pairwise_top10_overlap": {
            f"{row.approach_a}_vs_{row.approach_b}": int(row.top10_overlap)
            for row in metrics.itertuples()
        },
    }
    (out / f"{stratum}_conclusion.json").write_text(json.dumps(conclusion, indent=2))
    return conclusion


def compile_summary(strata: list[str] | tuple[str, ...]) -> pd.DataFrame:
    out = _output_dir()
    rows = []
    for stratum in strata:
        p = out / f"{stratum}_conclusion.json"
        if p.exists():
            rows.append(json.loads(p.read_text()))
    summary = pd.DataFrame(rows)
    summary.to_json(out / "cross_stratum_conclusions.json", orient="records", indent=2)

    top_rows = []
    for stratum in strata:
        p = out / f"{stratum}_consensus_ranking.csv"
        if not p.exists():
            continue
        d = pd.read_csv(p)
        d = d[d["n_approaches"].eq(3)].head(5)
        top_rows.extend({"stratum": stratum, "v_gene": v} for v in d["v_gene"])
    freq = pd.DataFrame(top_rows)
    if not freq.empty:
        overall = (
            freq.groupby("v_gene", as_index=False)
            .agg(n_strata_top5=("stratum", "nunique"))
            .sort_values(["n_strata_top5", "v_gene"], ascending=[False, True])
        )
    else:
        overall = pd.DataFrame(columns=["v_gene", "n_strata_top5"])
    overall.to_csv(out / "overall_consensus_v_genes.csv", index=False)

    if not overall.empty:
        top = overall.head(15).sort_values("n_strata_top5")
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.barh(top["v_gene"], top["n_strata_top5"])
        ax.set_xlabel("Biological strata with top-five consensus support")
        ax.set_ylabel("V segment")
        ax.set_title("Cross-stratum consensus across all three analytical approaches")
        save_figure(fig, APPROACH, None, "overall_cross_stratum_consensus")
    return summary
