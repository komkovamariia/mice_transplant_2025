"""Approach 1: exact aaV sets and count-based differential analysis."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib_venn import venn3
from scipy.stats import fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests

from .count_statistics import historical_fisher_log2_ratio
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
FDR_THRESHOLD = 0.05


def _output_dir() -> Path:
    path = repository_root() / "results" / APPROACH
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
        y <- calcNormFactors(y)
        design <- model.matrix(~group)
        keep <- filterByExpr(y, design=design)
        if (!any(keep)) {
            stop("No clonotypes passed edgeR expression filtering.")
        }
        y <- y[keep, , keep.lib.sizes=FALSE]
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
            [allogeneic_count, max(total_allogeneic - allogeneic_count, 0)],
        ]
        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
        odds_ratios.append(odds_ratio)
        p_values.append(p_value)
        log2_fold_changes.append(
            historical_fisher_log2_ratio(g1_count, allogeneic_count)
        )

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


def _empty_panel(ax, message: str) -> None:
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        transform=ax.transAxes,
        wrap=True,
    )


def _plot_set_overlap(sets: dict[str, set[str]], stratum: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    venn3(
        [sets["g1"], sets["g5"], sets["g6"]],
        set_labels=("g1", "g5", "g6"),
        ax=ax,
    )
    ax.set_title(f"Exact aaV clonotype overlap: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "01_exact_clonotype_overlap")


def _plot_exclusive_trav(set_table: pd.DataFrame, stratum: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    top = set_table.head(20).sort_values("n_g1_exclusive")
    if top.empty:
        _empty_panel(ax, "No g1-exclusive aaV clonotypes were observed.")
    else:
        ax.barh(top["v_gene"], top["n_g1_exclusive"], color="#4F6D7A")
        ax.set_xlabel("g1-exclusive aaV clonotypes")
        ax.set_ylabel("TRAV segment")
    ax.set_title(f"TRAV usage among g1-exclusive clonotypes: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "02_g1_exclusive_trav_usage")


def _volcano_data(features: pd.DataFrame) -> pd.DataFrame:
    data = (
        features.replace([np.inf, -np.inf], np.nan)
        .dropna(subset=["logFC", "FDR"])
        .copy()
    )
    data["minus_log10_fdr"] = -np.log10(data["FDR"].clip(lower=np.finfo(float).tiny))
    data["significant"] = data["FDR"].lt(FDR_THRESHOLD)
    return data


def _draw_volcano(ax, data: pd.DataFrame) -> None:
    if data.empty:
        _empty_panel(ax, "No edgeR-tested aaV clonotypes were available.")
        return
    colors = np.where(data["significant"], "#C44E52", "#A7A7A7")
    ax.scatter(data["logFC"], data["minus_log10_fdr"], s=14, alpha=0.65, c=colors)
    ax.axvline(1, linewidth=0.8, linestyle="--", color="black")
    ax.axvline(-1, linewidth=0.8, linestyle="--", color="black")
    ax.axhline(
        -np.log10(FDR_THRESHOLD),
        linewidth=0.8,
        linestyle="--",
        color="black",
    )
    ax.set_xlabel("edgeR log2 fold change: g1 / g5+g6")
    ax.set_ylabel("-log10 FDR")


def _plot_volcano(features: pd.DataFrame, stratum: str, *, labels: bool) -> None:
    data = _volcano_data(features)
    fig, ax = plt.subplots(figsize=(8, 6))
    _draw_volcano(ax, data)
    if labels and not data.empty:
        top = data.sort_values(["FDR", "logFC"], ascending=[True, False]).head(10)
        for row in top.itertuples():
            ax.annotate(
                str(row.feature_id),
                (row.logFC, row.minus_log10_fdr),
                fontsize=6,
                xytext=(3, 3),
                textcoords="offset points",
            )
    qualifier = " with leading clonotypes" if labels else ""
    ax.set_title(
        f"edgeR differential clonotype abundance{qualifier}: {STRATUM_LABELS[stratum]}"
    )
    stem = "04_edger_volcano_labeled" if labels else "03_edger_volcano"
    save_figure(fig, APPROACH, stratum, stem)


def _draw_ma(ax, data: pd.DataFrame) -> None:
    data = data.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["logCPM", "logFC", "FDR"]
    )
    if data.empty:
        _empty_panel(ax, "No edgeR-tested aaV clonotypes were available.")
        return
    colors = np.where(data["FDR"].lt(FDR_THRESHOLD), "#C44E52", "#A7A7A7")
    ax.scatter(data["logCPM"], data["logFC"], s=14, alpha=0.65, c=colors)
    ax.axhline(1, linewidth=0.8, linestyle="--", color="black")
    ax.axhline(-1, linewidth=0.8, linestyle="--", color="black")
    ax.set_xlabel("Mean clonotype abundance, logCPM")
    ax.set_ylabel("edgeR log2 fold change: g1 / g5+g6")


def _plot_ma(features: pd.DataFrame, stratum: str, *, labels: bool) -> None:
    data = features.copy()
    fig, ax = plt.subplots(figsize=(8, 6))
    _draw_ma(ax, data)
    if labels:
        finite = data.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["logCPM", "logFC", "FDR"]
        )
        top = finite.sort_values(["FDR", "logFC"], ascending=[True, False]).head(10)
        for row in top.itertuples():
            ax.annotate(
                str(row.feature_id),
                (row.logCPM, row.logFC),
                fontsize=6,
                xytext=(3, 3),
                textcoords="offset points",
            )
    qualifier = " with leading clonotypes" if labels else ""
    ax.set_title(f"edgeR MA plot{qualifier}: {STRATUM_LABELS[stratum]}")
    stem = "06_edger_ma_labeled" if labels else "05_edger_ma"
    save_figure(fig, APPROACH, stratum, stem)


def _plot_top_clonotypes(features: pd.DataFrame, stratum: str) -> None:
    significant = features[
        features["FDR"].lt(FDR_THRESHOLD) & features["logFC"].gt(0)
    ].copy()
    top = significant.sort_values(["logFC", "FDR"], ascending=[False, True]).head(15)
    top = top.sort_values("logFC")
    fig, ax = plt.subplots(figsize=(9, 7))
    if top.empty:
        _empty_panel(ax, "No g1-enriched aaV clonotypes passed edgeR FDR < 0.05.")
    else:
        ax.barh(top["feature_id"], top["logFC"], color="#4C72B0")
        ax.set_xlabel("edgeR log2 fold change: g1 / g5+g6")
        ax.set_ylabel("aaV clonotype")
        ax.tick_params(axis="y", labelsize=6)
    ax.set_title(f"Leading g1-enriched aaV clonotypes: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "07_top_g1_enriched_clonotypes")


def _plot_trav_effect(ranking: pd.DataFrame, stratum: str) -> None:
    data = ranking[
        ranking["n_edger_significant"].gt(0)
        & ranking["mean_effect_g1_vs_allogeneic"].notna()
    ].copy()
    top = data.sort_values(
        ["mean_effect_g1_vs_allogeneic", "n_edger_significant"],
        ascending=[False, False],
    ).head(15)
    top = top.sort_values("mean_effect_g1_vs_allogeneic")
    fig, ax = plt.subplots(figsize=(8, 6))
    if top.empty:
        _empty_panel(
            ax, "No TRAV segment contained an edgeR-significant aaV clonotype."
        )
    else:
        ax.barh(top["v_gene"], top["mean_effect_g1_vs_allogeneic"], color="#55A868")
        ax.set_xlabel("Mean edgeR log2 fold change among significant aaV clonotypes")
        ax.set_ylabel("TRAV segment")
    ax.set_title(f"TRAV segments ranked by mean effect: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "08_trav_mean_effect")


def _plot_trav_count(ranking: pd.DataFrame, stratum: str) -> None:
    top = ranking[ranking["n_edger_significant"].gt(0)].head(15).copy()
    top = top.sort_values("n_edger_significant")
    fig, ax = plt.subplots(figsize=(8, 6))
    if top.empty:
        _empty_panel(
            ax, "No TRAV segment contained an edgeR-significant aaV clonotype."
        )
    else:
        ax.barh(top["v_gene"], top["n_edger_significant"], color="#8172B3")
        ax.set_xlabel("edgeR-significant aaV clonotypes")
        ax.set_ylabel("TRAV segment")
    ax.set_title(
        f"TRAV segments ranked by significant-clonotype count: "
        f"{STRATUM_LABELS[stratum]}"
    )
    save_figure(fig, APPROACH, stratum, "09_trav_significant_feature_count")


def _plot_method_concordance(features: pd.DataFrame, stratum: str) -> None:
    data = features.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["logFC", "fisher_log2fc_g1_vs_allogeneic"]
    )
    fig, ax = plt.subplots(figsize=(7, 6))
    if data.empty:
        _empty_panel(
            ax, "No aaV clonotypes were jointly estimable by edgeR and Fisher testing."
        )
    else:
        colors = np.where(data["consensus_significant"], "#C44E52", "#6E6E6E")
        ax.scatter(
            data["fisher_log2fc_g1_vs_allogeneic"],
            data["logFC"],
            s=14,
            alpha=0.6,
            c=colors,
        )
        lower = float(
            min(data["fisher_log2fc_g1_vs_allogeneic"].min(), data["logFC"].min())
        )
        upper = float(
            max(data["fisher_log2fc_g1_vs_allogeneic"].max(), data["logFC"].max())
        )
        ax.plot(
            [lower, upper], [lower, upper], linestyle="--", color="black", linewidth=0.8
        )
        rho = spearmanr(
            data["fisher_log2fc_g1_vs_allogeneic"],
            data["logFC"],
        ).statistic
        ax.text(
            0.03, 0.97, f"Spearman rho = {rho:.2f}", transform=ax.transAxes, va="top"
        )
        ax.set_xlabel("Fisher log2 pooled-count ratio: g1 / g5+g6")
        ax.set_ylabel("edgeR log2 fold change: g1 / g5+g6")
    ax.set_title(f"Fisher and edgeR effect concordance: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "10_fisher_edger_concordance")


def _build_ranking(features: pd.DataFrame, set_table: pd.DataFrame) -> pd.DataFrame:
    edger_significant = features[features["FDR"].lt(FDR_THRESHOLD)]
    fisher_significant = features[features["fisher_fdr"].lt(FDR_THRESHOLD)]
    consensus = features[features["consensus_significant"]]

    edger_by_v = edger_significant.groupby("v_gene", as_index=False).agg(
        n_edger_significant=("feature_id", "size"),
        mean_effect_g1_vs_allogeneic=("logFC", "mean"),
        max_effect_g1_vs_allogeneic=("logFC", "max"),
        median_edger_fdr=("FDR", "median"),
    )
    fisher_by_v = fisher_significant.groupby("v_gene", as_index=False).agg(
        n_fisher_significant=("feature_id", "size"),
        mean_fisher_effect_g1_vs_allogeneic=(
            "fisher_log2fc_g1_vs_allogeneic",
            "mean",
        ),
        median_fisher_fdr=("fisher_fdr", "median"),
    )
    consensus_by_v = (
        consensus.groupby("v_gene", as_index=False)
        .size()
        .rename(columns={"size": "n_consensus"})
    )

    ranking = edger_by_v.merge(fisher_by_v, on="v_gene", how="outer")
    ranking = ranking.merge(consensus_by_v, on="v_gene", how="outer")
    ranking = ranking.merge(set_table, on="v_gene", how="outer")
    for column in (
        "n_edger_significant",
        "n_fisher_significant",
        "n_consensus",
        "n_g1_exclusive",
    ):
        ranking[column] = ranking[column].fillna(0).astype(int)

    ranking["evidence_score"] = (
        3.0 * ranking["n_consensus"]
        + 2.0 * ranking["n_edger_significant"]
        + ranking["n_fisher_significant"]
        + np.log1p(ranking["n_g1_exclusive"])
    )
    ranking["effect_direction"] = "undetermined"
    ranking.loc[
        ranking["mean_effect_g1_vs_allogeneic"].gt(0),
        "effect_direction",
    ] = "g1_enriched"
    ranking.loc[
        ranking["mean_effect_g1_vs_allogeneic"].lt(0),
        "effect_direction",
    ] = "g5_g6_enriched"
    missing_edger = ranking["mean_effect_g1_vs_allogeneic"].isna()
    ranking.loc[
        missing_edger & ranking["mean_fisher_effect_g1_vs_allogeneic"].gt(0),
        "effect_direction",
    ] = "g1_enriched"
    ranking.loc[
        missing_edger & ranking["mean_fisher_effect_g1_vs_allogeneic"].lt(0),
        "effect_direction",
    ] = "g5_g6_enriched"

    ranking["evidence_class"] = "set_exclusive_only"
    ranking.loc[ranking["n_fisher_significant"].gt(0), "evidence_class"] = (
        "fisher_supported"
    )
    ranking.loc[ranking["n_edger_significant"].gt(0), "evidence_class"] = (
        "edger_supported"
    )
    ranking.loc[ranking["n_consensus"].gt(0), "evidence_class"] = (
        "edger_fisher_consensus"
    )

    ranking["absolute_effect"] = (
        ranking["mean_effect_g1_vs_allogeneic"]
        .abs()
        .fillna(ranking["mean_fisher_effect_g1_vs_allogeneic"].abs())
    )
    ranking = ranking.sort_values(
        ["evidence_score", "absolute_effect", "n_g1_exclusive", "v_gene"],
        ascending=[False, False, False, True],
        na_position="last",
    ).reset_index(drop=True)
    ranking["rank"] = np.arange(1, len(ranking) + 1)
    ranking["score"] = ranking["evidence_score"]
    return ranking.drop(columns="absolute_effect")


def _result_paragraphs(payload: dict, top_genes: list[str]) -> list[str]:
    paragraphs = [
        (
            f"This stratum included {payload['n_analysis_units']} independent analysis "
            f"units from {payload['n_input_samples']} samples. Exact set subtraction "
            f"identified {payload['n_g1_exclusive_clonotypes']} aaV clonotypes present "
            "in g1 and absent from both g5 and g6."
        )
    ]
    if payload["n_edger_significant_clonotypes"]:
        paragraphs.append(
            f"edgeR identified {payload['n_edger_significant_clonotypes']} differential "
            f"aaV clonotypes at FDR < {FDR_THRESHOLD:.2f}; "
            f"{payload['n_fisher_edger_consensus']} also passed Fisher FDR control."
        )
    elif payload["n_fisher_significant_clonotypes"]:
        paragraphs.append(
            "No aaV clonotype passed edgeR FDR control, while Fisher testing identified "
            f"{payload['n_fisher_significant_clonotypes']} candidates. These signals "
            "require cautious interpretation because Fisher testing pools counts across "
            "analysis units."
        )
    else:
        paragraphs.append(
            "Neither edgeR nor Fisher testing identified a differential aaV clonotype "
            f"at FDR < {FDR_THRESHOLD:.2f}."
        )
    paragraphs.append(
        f"The highest-ranked distinctive TRAV segments were {format_gene_list(top_genes)}."
    )
    return paragraphs


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
    features["consensus_significant"] = features["FDR"].lt(FDR_THRESHOLD) & features[
        "fisher_fdr"
    ].lt(FDR_THRESHOLD)

    ranking = _build_ranking(features, set_table)
    ranking["stratum"] = stratum

    output = _output_dir()
    counts.to_csv(output / f"{stratum}_count_matrix.csv")
    features.to_csv(output / f"{stratum}_clonotype_statistics.csv", index=False)
    ranking.to_csv(output / f"{stratum}_trav_ranking.csv", index=False)
    set_table.to_csv(output / f"{stratum}_g1_exclusive_trav.csv", index=False)

    _plot_set_overlap(sets, stratum)
    _plot_exclusive_trav(set_table, stratum)
    _plot_volcano(features, stratum, labels=False)
    _plot_volcano(features, stratum, labels=True)
    _plot_ma(features, stratum, labels=False)
    _plot_ma(features, stratum, labels=True)
    _plot_top_clonotypes(features, stratum)
    _plot_trav_effect(ranking, stratum)
    _plot_trav_count(ranking, stratum)
    _plot_method_concordance(features, stratum)

    distinctive = ranking[ranking["evidence_score"].gt(0)]
    top_genes = distinctive.head(5)["v_gene"].dropna().tolist()
    payload = {
        "stratum": stratum,
        "label": STRATUM_LABELS[stratum],
        "n_input_samples": int(model_data["sample_id"].nunique()),
        "n_analysis_units": int(model_data["analysis_unit"].nunique()),
        "n_g1_exclusive_clonotypes": len(g1_exclusive),
        "n_edger_significant_clonotypes": int(features["FDR"].lt(FDR_THRESHOLD).sum()),
        "n_fisher_significant_clonotypes": int(
            features["fisher_fdr"].lt(FDR_THRESHOLD).sum()
        ),
        "n_fisher_edger_consensus": int(features["consensus_significant"].sum()),
        "top_trav_segments": top_genes,
    }
    return save_conclusion(
        output,
        stratum,
        f"Approach 1 conclusion: {STRATUM_LABELS[stratum]}",
        payload,
        _result_paragraphs(payload, top_genes),
    )


def _plot_cross_stratum_trav(summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11, 8))
    if summary.empty:
        _empty_panel(
            ax, "No distinctive TRAV segments were identified across the strata."
        )
    else:
        recurrence = (
            summary.groupby("v_gene")["stratum"].nunique().sort_values(ascending=False)
        )
        genes = recurrence.head(20).index
        matrix = (
            summary[summary["v_gene"].isin(genes)]
            .pivot_table(
                index="v_gene",
                columns="stratum",
                values="evidence_score",
                aggfunc="max",
                fill_value=0,
            )
            .reindex(index=genes)
        )
        image = ax.imshow(matrix.to_numpy(float), aspect="auto", cmap="viridis")
        ax.set_xticks(
            np.arange(len(matrix.columns)), matrix.columns, rotation=45, ha="right"
        )
        ax.set_yticks(np.arange(len(matrix.index)), matrix.index)
        ax.set_xlabel("Biological stratum")
        ax.set_ylabel("TRAV segment")
        fig.colorbar(image, ax=ax, label="Integrated count-based evidence score")
    ax.set_title("Distinctive TRAV evidence across eight biological strata")
    save_figure(fig, APPROACH, None, "distinctive_trav_across_eight_strata")


def compile_summary(strata) -> pd.DataFrame:
    rows = []
    ranking_rows = []
    output = _output_dir()
    for stratum in strata:
        conclusion_path = output / f"{stratum}_conclusion.json"
        ranking_path = output / f"{stratum}_trav_ranking.csv"
        if conclusion_path.exists():
            rows.append(json.loads(conclusion_path.read_text(encoding="utf-8")))
        if ranking_path.exists():
            ranking = pd.read_csv(ranking_path)
            ranking = ranking[ranking["evidence_score"].gt(0)].copy()
            ranking["stratum_label"] = STRATUM_LABELS[stratum]
            ranking_rows.append(ranking)

    summary = pd.DataFrame(rows)
    csv_summary = summary.copy()
    if "top_trav_segments" in csv_summary:
        csv_summary["top_trav_segments"] = csv_summary["top_trav_segments"].apply(
            lambda values: "; ".join(values)
        )
    csv_summary.to_csv(output / "eight_stratum_summary.csv", index=False)

    if ranking_rows:
        distinctive = pd.concat(ranking_rows, ignore_index=True)
        recurrence = distinctive.groupby("v_gene").agg(
            n_strata_with_evidence=("stratum", "nunique"),
            best_rank=("rank", "min"),
            total_evidence_score=("evidence_score", "sum"),
        )
        distinctive = distinctive.merge(recurrence, on="v_gene", how="left")
        distinctive = distinctive.sort_values(
            ["n_strata_with_evidence", "total_evidence_score", "stratum", "rank"],
            ascending=[False, False, True, True],
        )
    else:
        distinctive = pd.DataFrame(
            columns=[
                "v_gene",
                "stratum",
                "stratum_label",
                "rank",
                "evidence_score",
                "evidence_class",
                "effect_direction",
                "n_strata_with_evidence",
            ]
        )
    distinctive.to_csv(output / "distinctive_trav_summary.csv", index=False)
    _plot_cross_stratum_trav(distinctive)

    if distinctive.empty:
        leading = []
    else:
        leading = (
            distinctive.sort_values(
                ["n_strata_with_evidence", "total_evidence_score"],
                ascending=[False, False],
            )["v_gene"]
            .drop_duplicates()
            .head(10)
            .tolist()
        )
    final_text = (
        "# Approach 1: cross-stratum conclusion\n\n"
        "Eight prespecified biological strata were analyzed independently. The TRAV "
        f"segments with the broadest count-based support were "
        f"{format_gene_list(leading, limit=10)}. Detailed evidence, effect direction, "
        "and stratum-specific rank are retained in `distinctive_trav_summary.csv`; "
        "the eight separate conclusion files should be used for biological "
        "interpretation of compartment-specific signals.\n"
    )
    (output / "cross_stratum_conclusion.md").write_text(final_text, encoding="utf-8")
    return summary
