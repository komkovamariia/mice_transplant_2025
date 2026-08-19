"""Approach 3: clone-level prioritization from independent evidence layers."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .figures import save_figure
from .reporting import format_gene_list, save_conclusion
from .strata import (
    STRATUM_LABELS,
    ensure_groups,
    repository_root,
    select_stratum,
)

APPROACH = "03_clone_alloreactivity"

HYDROPHOBICITY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}


def _output_dir() -> Path:
    path = repository_root() / "outputs" / "tables" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _sequence_features(cdr3: pd.Series) -> pd.DataFrame:
    def mean_hydrophobicity(sequence: str) -> float:
        values = [HYDROPHOBICITY.get(residue, 0.0) for residue in sequence]
        return float(np.mean(values)) if values else np.nan

    def approximate_charge(sequence: str) -> float:
        positive = sum(residue in "KR" for residue in sequence)
        negative = sum(residue in "DE" for residue in sequence)
        return float(positive - negative)

    return pd.DataFrame(
        {
            "cdr3_length": cdr3.str.len().astype(float),
            "mean_hydrophobicity": cdr3.map(mean_hydrophobicity),
            "approximate_net_charge": cdr3.map(approximate_charge),
        },
        index=cdr3.index,
    )


def _convergence_degree(features: pd.DataFrame) -> pd.Series:
    """Count same-V, same-length CDR3 neighbors at Hamming distance one."""
    degree = pd.Series(0, index=features.index, dtype=int)
    for (_, length), block in features.groupby(["v_gene", features["cdr3"].str.len()]):
        wildcard_members = defaultdict(set)
        for index, sequence in zip(block.index, block["cdr3"]):
            for position in range(int(length)):
                pattern = sequence[:position] + "*" + sequence[position + 1 :]
                wildcard_members[pattern].add(index)

        neighbors = defaultdict(set)
        for members in wildcard_members.values():
            if len(members) < 2:
                continue
            for index in members:
                neighbors[index].update(members - {index})
        for index, members in neighbors.items():
            degree.loc[index] = len(members)
    return degree


def _plot_candidate_scores(
    candidates: pd.DataFrame,
    stratum: str,
) -> None:
    top = candidates.head(20).sort_values("integrated_score")
    if not top.empty:
        labels = top["v_gene"] + " | " + top["cdr3"].str.slice(0, 14)
        fig, ax = plt.subplots(figsize=(9, 7))
        ax.barh(labels, top["integrated_score"])
        ax.set_xlabel("Integrated clone-evidence score")
        ax.set_ylabel("Candidate clonotype")
        ax.set_title(f"Clone-level candidate ranking: {STRATUM_LABELS[stratum]}")
        save_figure(
            fig,
            APPROACH,
            stratum,
            "integrated_candidate_ranking",
        )

    by_v_gene = (
        candidates.groupby("v_gene", as_index=False)
        .agg(
            candidate_count=("ckey", "size"),
            mean_integrated_score=("integrated_score", "mean"),
            max_integrated_score=("integrated_score", "max"),
        )
        .sort_values(
            ["candidate_count", "mean_integrated_score"],
            ascending=False,
        )
        .head(15)
        .sort_values("candidate_count")
    )
    if not by_v_gene.empty:
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.barh(
            by_v_gene["v_gene"],
            by_v_gene["candidate_count"],
        )
        ax.set_xlabel("Prioritized clonotypes")
        ax.set_ylabel("V segment")
        ax.set_title(f"Clone-level V-segment support: {STRATUM_LABELS[stratum]}")
        save_figure(
            fig,
            APPROACH,
            stratum,
            "v_segment_candidate_support",
        )


def _group_pivot(
    summary: pd.DataFrame,
    value: str,
) -> pd.DataFrame:
    pivot = summary.pivot_table(
        index=["ckey", "cdr3", "v_gene"],
        columns="group",
        values=value,
        fill_value=0,
    )
    for group in ("g1", "g5", "g6"):
        if group not in pivot:
            pivot[group] = 0
    return pivot[["g1", "g5", "g6"]]


def run_stratum(
    repertoire: pd.DataFrame,
    stratum: str,
) -> dict:
    data = select_stratum(repertoire, stratum)
    ensure_groups(data, ("g1", "g5", "g6"), f"Approach 3 / {stratum}")
    data = data[data["group"].isin(["g1", "g5", "g6"])].copy()

    per_mouse = data.groupby(
        ["ckey", "cdr3", "v_gene", "group", "mouse_id"],
        as_index=False,
    )["umi"].sum()
    per_mouse["mouse_library_umi"] = per_mouse.groupby(["group", "mouse_id"])[
        "umi"
    ].transform("sum")
    per_mouse["relative_abundance"] = per_mouse["umi"] / per_mouse["mouse_library_umi"]

    group_summary = per_mouse.groupby(
        ["ckey", "cdr3", "v_gene", "group"],
        as_index=False,
    ).agg(
        mice=("mouse_id", "nunique"),
        umi=("umi", "sum"),
        relative_abundance_sum=("relative_abundance", "sum"),
    )
    mouse_totals = data.groupby("group")["mouse_id"].nunique().to_dict()
    n_g1 = int(mouse_totals.get("g1", 0))
    n_g5 = int(mouse_totals.get("g5", 0))
    n_g6 = int(mouse_totals.get("g6", 0))
    n_allogeneic = n_g5 + n_g6
    if n_g1 < 1 or n_allogeneic < 2:
        raise ValueError(
            f"Approach 3 / {stratum}: at least one g1 mouse and two "
            "allogeneic mice are required."
        )

    mice = _group_pivot(group_summary, "mice").rename(
        columns={
            "g1": "g1_mice",
            "g5": "g5_mice",
            "g6": "g6_mice",
        }
    )
    umi = _group_pivot(group_summary, "umi").rename(
        columns={
            "g1": "g1_umi",
            "g5": "g5_umi",
            "g6": "g6_umi",
        }
    )
    relative = _group_pivot(
        group_summary,
        "relative_abundance_sum",
    ).rename(
        columns={
            "g1": "g1_relative_abundance_sum",
            "g5": "g5_relative_abundance_sum",
            "g6": "g6_relative_abundance_sum",
        }
    )

    statistics = mice.join(umi).join(relative).reset_index()
    statistics["allogeneic_mice"] = statistics["g5_mice"] + statistics["g6_mice"]
    statistics["allogeneic_umi"] = statistics["g5_umi"] + statistics["g6_umi"]
    statistics["g1_mouse_fraction"] = statistics["g1_mice"] / n_g1
    statistics["allogeneic_mouse_fraction"] = (
        statistics["allogeneic_mice"] / n_allogeneic
    )
    statistics["prevalence_delta_allogeneic_vs_g1"] = (
        statistics["allogeneic_mouse_fraction"] - statistics["g1_mouse_fraction"]
    )

    statistics["g1_mean_relative_abundance"] = (
        statistics["g1_relative_abundance_sum"] / n_g1
    )
    statistics["allogeneic_mean_relative_abundance"] = (
        statistics["g5_relative_abundance_sum"]
        + statistics["g6_relative_abundance_sum"]
    ) / n_allogeneic
    median_library = float(
        per_mouse[["group", "mouse_id", "mouse_library_umi"]]
        .drop_duplicates()["mouse_library_umi"]
        .median()
    )
    pseudocount = 0.5 / max(median_library, 1.0)
    statistics["log2_relative_abundance_fc_allogeneic_vs_g1"] = np.log2(
        (statistics["allogeneic_mean_relative_abundance"] + pseudocount)
        / (statistics["g1_mean_relative_abundance"] + pseudocount)
    )
    statistics["effect_g1_vs_allogeneic"] = -statistics[
        "log2_relative_abundance_fc_allogeneic_vs_g1"
    ]

    sequence_features = _sequence_features(statistics["cdr3"])
    statistics = pd.concat(
        [statistics, sequence_features],
        axis=1,
    )
    statistics["convergence_degree"] = _convergence_degree(statistics)

    statistics["score_prevalence"] = statistics[
        "prevalence_delta_allogeneic_vs_g1"
    ].rank(pct=True)
    statistics["score_abundance"] = statistics[
        "log2_relative_abundance_fc_allogeneic_vs_g1"
    ].rank(pct=True)
    statistics["score_convergence"] = statistics["convergence_degree"].rank(pct=True)
    statistics["integrated_score"] = (
        0.45 * statistics["score_prevalence"]
        + 0.35 * statistics["score_abundance"]
        + 0.20 * statistics["score_convergence"]
    )

    candidates = statistics[statistics["allogeneic_mice"] >= 2].copy()
    candidates = candidates.sort_values(
        [
            "integrated_score",
            "allogeneic_mice",
            "prevalence_delta_allogeneic_vs_g1",
        ],
        ascending=False,
    ).reset_index(drop=True)
    candidates["rank"] = np.arange(1, len(candidates) + 1)
    candidates["score"] = candidates["integrated_score"]
    candidates["stratum"] = stratum

    v_gene_ranking = (
        candidates.groupby("v_gene", as_index=False)
        .agg(
            candidate_count=("ckey", "size"),
            mean_integrated_score=("integrated_score", "mean"),
            max_integrated_score=("integrated_score", "max"),
            mean_prevalence_delta_allogeneic_vs_g1=(
                "prevalence_delta_allogeneic_vs_g1",
                "mean",
            ),
            mean_log2_relative_abundance_fc_allogeneic_vs_g1=(
                "log2_relative_abundance_fc_allogeneic_vs_g1",
                "mean",
            ),
            mean_effect_g1_vs_allogeneic=(
                "effect_g1_vs_allogeneic",
                "mean",
            ),
            mean_convergence_degree=("convergence_degree", "mean"),
        )
        .sort_values(
            [
                "candidate_count",
                "mean_integrated_score",
                "max_integrated_score",
            ],
            ascending=False,
        )
        .reset_index(drop=True)
    )
    v_gene_ranking["rank"] = np.arange(
        1,
        len(v_gene_ranking) + 1,
    )
    v_gene_ranking["score"] = v_gene_ranking["candidate_count"].astype(float)
    v_gene_ranking["stratum"] = stratum

    output = _output_dir()
    statistics.to_csv(
        output / f"{stratum}_all_clonotype_evidence.csv",
        index=False,
    )
    candidates.to_csv(
        output / f"{stratum}_candidate_clonotypes.csv",
        index=False,
    )
    v_gene_ranking.to_csv(
        output / f"{stratum}_v_gene_ranking.csv",
        index=False,
    )
    _plot_candidate_scores(candidates, stratum)

    top = candidates.head(20)
    top_genes = v_gene_ranking.head(5)["v_gene"].dropna().tolist()
    payload = {
        "stratum": stratum,
        "label": STRATUM_LABELS[stratum],
        "n_input_samples": int(data["sample_id"].nunique()),
        "n_analysis_units": int(data["analysis_unit"].nunique()),
        "n_mice_g1": n_g1,
        "n_mice_allogeneic": n_allogeneic,
        "n_candidate_clonotypes": len(candidates),
        "n_top20_absent_from_g1": int(top["g1_mice"].eq(0).sum()),
        "top_v_genes": top_genes,
        "top_candidate_clonotypes": top["ckey"].head(5).tolist(),
        "abundance_normalization": (
            "mean within-mouse relative abundance; mice weighted equally"
        ),
    }
    return save_conclusion(
        output,
        stratum,
        f"Approach 3 conclusion: {STRATUM_LABELS[stratum]}",
        payload,
        [
            (
                f"Approach 3 evaluated {n_g1} g1 mice and "
                f"{n_allogeneic} allogeneic mice using prevalence, "
                "mouse-normalized abundance, and local CDR3 convergence."
            ),
            (
                f"{len(candidates)} clonotypes were observed in at least "
                "two allogeneic mice; "
                f"{payload['n_top20_absent_from_g1']} of the top 20 were "
                "not observed in g1."
            ),
            (
                "The most strongly supported V segments were "
                f"{format_gene_list(top_genes)}."
            ),
            (
                "The integrated score is an evidence ranking and is not "
                "a calibrated probability of antigen specificity."
            ),
        ],
    )


def compile_summary(
    strata: list[str] | tuple[str, ...],
) -> pd.DataFrame:
    rows = []
    output = _output_dir()
    for stratum in strata:
        path = output / f"{stratum}_conclusion.json"
        if path.exists():
            rows.append(json.loads(path.read_text(encoding="utf-8")))

    summary = pd.DataFrame(rows)
    summary.to_csv(output / "stratum_summary.csv", index=False)
    return summary
