"""Approach 1: exact clonotype sets and count-based differential analysis."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import venn3
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from .figures import save_figure
from .reporting import format_gene_list, save_conclusion
from .strata import (
    STRATUM_LABELS,
    analysis_metadata,
    ensure_groups,
    repository_root,
    select_stratum,
)

APPROACH = "01_set_count"


def _output_dir() -> Path:
    path = repository_root() / "outputs" / "tables" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _count_matrix(df: pd.DataFrame) -> pd.DataFrame:
    matrix = df.pivot_table(
        index="ckey",
        columns="analysis_unit",
        values="umi",
        aggfunc="sum",
        fill_value=0,
    )
    return matrix.round().astype(int)


def _edgeR(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    """Fit the prespecified edgeR quasi-likelihood contrast: g1 / g5+g6."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    if set(groups.unique()) != {"allogeneic", "g1"}:
        raise ValueError(
            "edgeR requires g1 and allogeneic analysis units; "
            f"observed {sorted(groups.unique())}."
        )
    if (groups == "g1").sum() < 2 or (groups == "allogeneic").sum() < 2:
        raise ValueError(
            "edgeR requires at least two independent analysis units in both groups."
        )

    ordered_groups = pd.Categorical(
        groups,
        categories=["allogeneic", "g1"],
        ordered=True,
    )
    group_frame = pd.DataFrame(
        {"contrast_group": ordered_groups.astype(str)},
        index=counts.columns,
    )

    with (ro.default_converter + pandas2ri.converter).context():
        converter = ro.conversion.get_conversion()
        ro.globalenv["countData"] = converter.py2rpy(counts)
        ro.globalenv["coldata"] = converter.py2rpy(group_frame.reset_index(drop=True))

    ro.r(
        """
        suppressPackageStartupMessages(library(edgeR))
        countData <- as.matrix(countData)
        storage.mode(countData) <- "integer"
        group <- factor(
            coldata$contrast_group,
            levels=c("allogeneic", "g1")
        )
        y <- DGEList(counts=countData, group=group)
        design <- model.matrix(~group)
        keep <- filterByExpr(y, design=design)
        if (!any(keep)) {
            stop("No clonotypes passed edgeR expression filtering.")
        }
        y <- y[keep, , keep.lib.sizes=FALSE]
        y <- calcNormFactors(y)
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design)
        test <- glmQLFTest(fit, coef=2)
        edgeR_table <- topTags(test, n=Inf)$table
        """
    )

    with (ro.default_converter + pandas2ri.converter).context():
        output = ro.conversion.get_conversion().rpy2py(ro.globalenv["edgeR_table"])
    output = pd.DataFrame(output)
    output["feature_id"] = output.index.astype(str)
    return output.reset_index(drop=True)


def _fisher(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    g1_columns = groups.index[groups.eq("g1")]
    allogeneic_columns = groups.index[groups.eq("allogeneic")]
    g1 = counts[g1_columns].sum(axis=1).astype(float)
    allogeneic = counts[allogeneic_columns].sum(axis=1).astype(float)
    total_g1 = float(g1.sum())
    total_allogeneic = float(allogeneic.sum())

    p_values = []
    odds_ratios = []
    log2_fold_changes = []
    for g1_count, allogeneic_count in zip(g1.values, allogeneic.values):
        table = [
            [g1_count, max(total_g1 - g1_count, 0)],
            [
                allogeneic_count,
                max(total_allogeneic - allogeneic_count, 0),
            ],
        ]
        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
        g1_frequency = (g1_count + 0.5) / (total_g1 + 1.0)
        allogeneic_frequency = (allogeneic_count + 0.5) / (total_allogeneic + 1.0)
        odds_ratios.append(odds_ratio)
        p_values.append(p_value)
        log2_fold_changes.append(np.log2(g1_frequency / allogeneic_frequency))

    fdr = (
        multipletests(p_values, method="fdr_bh")[1]
        if p_values
        else np.array([], dtype=float)
    )
    return pd.DataFrame(
        {
            "feature_id": counts.index.astype(str),
            "fisher_log2fc_g1_vs_allogeneic": log2_fold_changes,
            "fisher_odds_ratio": odds_ratios,
            "fisher_p": p_values,
            "fisher_fdr": fdr,
        }
    )


def _v_gene_from_feature(feature: pd.Series) -> pd.Series:
    return feature.astype(str).str.rsplit("|", n=1).str[-1]


def _plot_set_overlap(sets: dict[str, set[str]], stratum: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    venn3(
        [sets["g1"], sets["g5"], sets["g6"]],
        set_labels=("g1", "g5", "g6"),
        ax=ax,
    )
    ax.set_title(f"Exact aaV clonotype overlap: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "exact_clonotype_overlap")


def _plot_volcano(features: pd.DataFrame, stratum: str) -> None:
    data = (
        features.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["logFC", "FDR"])
        .copy()
    )
    if data.empty:
        return

    data["minus_log10_fdr"] = -np.log10(data["FDR"].clip(lower=np.finfo(float).tiny))
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(
        data["logFC"],
        data["minus_log10_fdr"],
        s=9,
        alpha=0.55,
    )
    ax.axvline(0, linewidth=0.8)
    ax.axhline(-np.log10(0.05), linewidth=0.8, linestyle="--")
    ax.set_xlabel("edgeR log2 fold change: g1 / allogeneic")
    ax.set_ylabel("-log10 FDR")
    ax.set_title(f"Differential clonotype abundance: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "edger_volcano")


def _plot_ranking(ranking: pd.DataFrame, stratum: str) -> None:
    top = ranking.head(15).sort_values(
        ["n_edger_significant", "mean_effect_g1_vs_allogeneic"],
        ascending=True,
    )
    if top.empty:
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(top["v_gene"], top["n_edger_significant"])
    ax.set_xlabel("Significant aaV clonotypes: edgeR FDR < 0.05")
    ax.set_ylabel("V segment")
    ax.set_title(f"Count-based V-segment ranking: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "v_segment_ranking")


def run_stratum(repertoire: pd.DataFrame, stratum: str) -> dict:
    data = select_stratum(repertoire, stratum)
    ensure_groups(data, ("g1", "g5", "g6"), f"Approach 1 / {stratum}")

    sets = {
        group: set(data.loc[data["group"].eq(group), "ckey"].unique())
        for group in ("g1", "g5", "g6")
    }
    g1_exclusive = sets["g1"] - sets["g5"] - sets["g6"]
    set_table = (
        pd.DataFrame({"feature_id": sorted(g1_exclusive)})
        .assign(v_gene=lambda frame: _v_gene_from_feature(frame["feature_id"]))
        .groupby("v_gene", as_index=False)
        .size()
        .rename(columns={"size": "n_g1_exclusive"})
        .sort_values("n_g1_exclusive", ascending=False)
    )

    model_data = data[data["group"].isin(["g1", "g5", "g6"])].copy()
    counts = _count_matrix(model_data)
    metadata = analysis_metadata(model_data).loc[counts.columns]
    groups = metadata["group"].replace({"g5": "allogeneic", "g6": "allogeneic"})

    edger = _edgeR(counts, groups)
    fisher = _fisher(counts, groups)
    features = edger.merge(fisher, on="feature_id", how="outer")
    features["v_gene"] = _v_gene_from_feature(features["feature_id"])
    features["stratum"] = stratum
    features["effect_g1_vs_allogeneic"] = features["logFC"]
    features["consensus_significant"] = features["FDR"].lt(0.05) & features[
        "fisher_fdr"
    ].lt(0.05)

    significant = features[features["FDR"].lt(0.05)].copy()
    ranking = (
        significant.groupby("v_gene", as_index=False)
        .agg(
            n_edger_significant=("feature_id", "size"),
            n_consensus=("consensus_significant", "sum"),
            mean_effect_g1_vs_allogeneic=(
                "effect_g1_vs_allogeneic",
                "mean",
            ),
            max_effect_g1_vs_allogeneic=(
                "effect_g1_vs_allogeneic",
                "max",
            ),
            median_fdr=("FDR", "median"),
        )
        .merge(set_table, on="v_gene", how="outer")
        .fillna(
            {
                "n_edger_significant": 0,
                "n_consensus": 0,
                "n_g1_exclusive": 0,
            }
        )
    )
    ranking = ranking.sort_values(
        [
            "n_edger_significant",
            "mean_effect_g1_vs_allogeneic",
            "n_g1_exclusive",
        ],
        ascending=[False, False, False],
        na_position="last",
    ).reset_index(drop=True)
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    ranking["score"] = ranking["n_edger_significant"].astype(float)
    ranking["stratum"] = stratum

    output = _output_dir()
    features.to_csv(
        output / f"{stratum}_clonotype_statistics.csv",
        index=False,
    )
    ranking.to_csv(
        output / f"{stratum}_v_gene_ranking.csv",
        index=False,
    )
    set_table.to_csv(
        output / f"{stratum}_g1_exclusive_v_genes.csv",
        index=False,
    )
    _plot_set_overlap(sets, stratum)
    _plot_volcano(features, stratum)
    _plot_ranking(ranking, stratum)

    top_genes = (
        ranking[ranking["n_edger_significant"].gt(0)]
        .head(5)["v_gene"]
        .dropna()
        .tolist()
    )
    payload = {
        "stratum": stratum,
        "label": STRATUM_LABELS[stratum],
        "n_input_samples": int(model_data["sample_id"].nunique()),
        "n_analysis_units": int(model_data["analysis_unit"].nunique()),
        "n_g1_exclusive_clonotypes": len(g1_exclusive),
        "n_edger_significant_clonotypes": int(features["FDR"].lt(0.05).sum()),
        "n_fisher_edger_consensus": int(features["consensus_significant"].sum()),
        "top_v_genes": top_genes,
    }
    return save_conclusion(
        output,
        stratum,
        f"Approach 1 conclusion: {STRATUM_LABELS[stratum]}",
        payload,
        [
            (
                f"Approach 1 analyzed {payload['n_analysis_units']} independent "
                f"units and identified {payload['n_g1_exclusive_clonotypes']} "
                "exact aaV clonotypes observed in g1 and absent from g5 and g6."
            ),
            (
                f"edgeR detected {payload['n_edger_significant_clonotypes']} "
                "clonotypes at FDR < 0.05; "
                f"{payload['n_fisher_edger_consensus']} were also supported "
                "by Fisher's exact test."
            ),
            (f"The highest-ranked V segments were {format_gene_list(top_genes)}."),
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
