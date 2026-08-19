"""Approach 3: clone-level alloreactivity prioritization from independent evidence layers."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .figures import save_figure
from .strata import STRATUM_LABELS, ensure_groups, repository_root, select_stratum

APPROACH = "03_clone_alloreactivity"

HYDROPHOBICITY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8, "G": -0.4,
    "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5,
    "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8, "T": -0.7, "V": 4.2,
    "W": -0.9, "Y": -1.3,
}


def _output_dir() -> Path:
    path = repository_root() / "outputs" / "tables" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _sequence_features(cdr3: pd.Series) -> pd.DataFrame:
    def hydro(seq: str) -> float:
        vals = [HYDROPHOBICITY.get(a, 0.0) for a in seq]
        return float(np.mean(vals)) if vals else np.nan

    def charge(seq: str) -> float:
        return float(sum(a in "KR" for a in seq) - sum(a in "DE" for a in seq))

    return pd.DataFrame({
        "cdr3_length": cdr3.str.len().astype(float),
        "mean_hydrophobicity": cdr3.map(hydro),
        "net_charge": cdr3.map(charge),
    }, index=cdr3.index)


def _convergence_degree(features: pd.DataFrame) -> pd.Series:
    """Count same-V, same-length CDR3 neighbors at Hamming distance <= 1."""
    degree = pd.Series(0, index=features.index, dtype=int)
    for (_, length), block in features.groupby(["v_gene", features["cdr3"].str.len()]):
        pattern_members = defaultdict(set)
        for idx, seq in zip(block.index, block["cdr3"]):
            for pos in range(int(length)):
                pattern_members[seq[:pos] + "*" + seq[pos + 1:]].add(idx)
        neighbors = defaultdict(set)
        for members in pattern_members.values():
            if len(members) < 2:
                continue
            for idx in members:
                neighbors[idx].update(members - {idx})
        for idx, members in neighbors.items():
            degree.loc[idx] = len(members)
    return degree


def _plot_candidate_scores(candidates: pd.DataFrame, stratum: str) -> None:
    top = candidates.head(20).sort_values("integrated_score")
    if not top.empty:
        labels = top["v_gene"] + " | " + top["cdr3"].str.slice(0, 14)
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.barh(labels, top["integrated_score"])
        ax.set_xlabel("Integrated alloreactivity score")
        ax.set_ylabel("Candidate clonotype")
        ax.set_title(f"Clone-level candidate ranking: {STRATUM_LABELS[stratum]}")
        save_figure(fig, APPROACH, stratum, "integrated_candidate_ranking")

    by_v = (
        candidates.groupby("v_gene", as_index=False)
        .agg(
            candidate_count=("ckey", "size"),
            mean_integrated_score=("integrated_score", "mean"),
            max_integrated_score=("integrated_score", "max"),
        )
        .sort_values(["candidate_count", "mean_integrated_score"], ascending=False)
        .head(15)
        .sort_values("candidate_count")
    )
    if not by_v.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.barh(by_v["v_gene"], by_v["candidate_count"])
        ax.set_xlabel("Prioritized clonotypes")
        ax.set_ylabel("V segment")
        ax.set_title(f"Clone-level V-segment support: {STRATUM_LABELS[stratum]}")
        save_figure(fig, APPROACH, stratum, "v_segment_candidate_support")


def run_stratum(repertoire: pd.DataFrame, stratum: str) -> dict:
    sub = select_stratum(repertoire, stratum)
    ensure_groups(sub, ("g1", "g5", "g6"), f"Approach 3 / {stratum}")
    sub = sub[sub["group"].isin(["g1", "g5", "g6"])].copy()

    per_mouse = (
        sub.groupby(["ckey", "cdr3", "v_gene", "group", "mouse_id"], as_index=False)["umi"]
        .sum()
    )
    group_summary = (
        per_mouse.groupby(["ckey", "cdr3", "v_gene", "group"], as_index=False)
        .agg(mice=("mouse_id", "nunique"), umi=("umi", "sum"))
    )
    mouse_totals = sub.groupby("group")["mouse_id"].nunique().to_dict()

    pivot_mice = group_summary.pivot_table(
        index=["ckey", "cdr3", "v_gene"], columns="group", values="mice", fill_value=0
    )
    pivot_umi = group_summary.pivot_table(
        index=["ckey", "cdr3", "v_gene"], columns="group", values="umi", fill_value=0
    )
    for frame in (pivot_mice, pivot_umi):
        for g in ("g1", "g5", "g6"):
            if g not in frame:
                frame[g] = 0

    stats = pivot_mice[["g1", "g5", "g6"]].rename(
        columns={"g1": "g1_mice", "g5": "g5_mice", "g6": "g6_mice"}
    )
    stats = stats.join(
        pivot_umi[["g1", "g5", "g6"]].rename(
            columns={"g1": "g1_umi", "g5": "g5_umi", "g6": "g6_umi"}
        )
    ).reset_index()
    stats["allo_mice"] = stats["g5_mice"] + stats["g6_mice"]
    stats["allo_umi"] = stats["g5_umi"] + stats["g6_umi"]
    n_g1 = max(int(mouse_totals.get("g1", 0)), 1)
    n_allo = max(int(mouse_totals.get("g5", 0) + mouse_totals.get("g6", 0)), 1)
    stats["g1_mouse_fraction"] = stats["g1_mice"] / n_g1
    stats["allo_mouse_fraction"] = stats["allo_mice"] / n_allo
    stats["prevalence_delta"] = stats["allo_mouse_fraction"] - stats["g1_mouse_fraction"]
    stats["log2_umi_fc"] = np.log2((stats["allo_umi"] + 0.5) / (stats["g1_umi"] + 0.5))

    seq = _sequence_features(stats["cdr3"])
    stats = pd.concat([stats, seq], axis=1)
    stats["convergence_degree"] = _convergence_degree(stats)

    # Empirical percentiles place heterogeneous evidence layers on a common scale.
    stats["score_prevalence"] = stats["prevalence_delta"].rank(pct=True)
    stats["score_abundance"] = stats["log2_umi_fc"].rank(pct=True)
    stats["score_convergence"] = stats["convergence_degree"].rank(pct=True)
    stats["integrated_score"] = (
        0.45 * stats["score_prevalence"]
        + 0.35 * stats["score_abundance"]
        + 0.20 * stats["score_convergence"]
    )

    # Observation in at least two allogeneic mice is required. Absence from g1 is
    # supporting evidence, not a mandatory definition of alloreactivity.
    candidates = stats[stats["allo_mice"] >= 2].copy()
    candidates = candidates.sort_values(
        ["integrated_score", "allo_mice", "prevalence_delta"],
        ascending=False,
    ).reset_index(drop=True)
    candidates["rank"] = np.arange(1, len(candidates) + 1)
    candidates["score"] = candidates["integrated_score"]
    candidates["stratum"] = stratum

    v_rank = (
        candidates.groupby("v_gene", as_index=False)
        .agg(
            candidate_count=("ckey", "size"),
            mean_integrated_score=("integrated_score", "mean"),
            max_integrated_score=("integrated_score", "max"),
            mean_prevalence_delta=("prevalence_delta", "mean"),
            mean_log2_umi_fc=("log2_umi_fc", "mean"),
            mean_convergence_degree=("convergence_degree", "mean"),
        )
        .sort_values(
            ["candidate_count", "mean_integrated_score", "max_integrated_score"],
            ascending=False,
        )
        .reset_index(drop=True)
    )
    v_rank["rank"] = np.arange(1, len(v_rank) + 1)
    v_rank["score"] = v_rank["candidate_count"].astype(float)
    v_rank["stratum"] = stratum

    out = _output_dir()
    stats.to_csv(out / f"{stratum}_all_clonotype_evidence.csv", index=False)
    candidates.to_csv(out / f"{stratum}_candidate_clonotypes.csv", index=False)
    v_rank.to_csv(out / f"{stratum}_v_gene_ranking.csv", index=False)
    _plot_candidate_scores(candidates, stratum)

    top = candidates.head(20)
    conclusion = {
        "stratum": stratum,
        "label": STRATUM_LABELS[stratum],
        "n_samples": int(sub["sample_id"].nunique()),
        "n_mice_g1": int(n_g1),
        "n_mice_allogeneic": int(n_allo),
        "n_candidate_clonotypes": int(len(candidates)),
        "n_top20_absent_from_g1": int(top["g1_mice"].eq(0).sum()),
        "top_v_genes": v_rank.head(5)["v_gene"].dropna().tolist(),
        "top_candidate_clonotypes": top["ckey"].head(5).tolist(),
    }
    (out / f"{stratum}_conclusion.json").write_text(json.dumps(conclusion, indent=2))
    return conclusion


def compile_summary(strata: list[str] | tuple[str, ...]) -> pd.DataFrame:
    rows = []
    out = _output_dir()
    for stratum in strata:
        path = out / f"{stratum}_conclusion.json"
        if path.exists():
            rows.append(json.loads(path.read_text()))
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "stratum_summary.csv", index=False)
    return summary
