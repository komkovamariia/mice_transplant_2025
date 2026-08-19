"""Approach 1: exact clonotype set operations and count-based differential analysis."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from .figures import save_figure
from .strata import STRATUM_LABELS, ensure_groups, repository_root, sample_metadata, select_stratum

APPROACH = "01_set_count"


def _output_dir() -> Path:
    path = repository_root() / "outputs" / "tables" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _count_matrix(df: pd.DataFrame) -> pd.DataFrame:
    matrix = df.pivot_table(index="ckey", columns="sample_id", values="umi", aggfunc="sum", fill_value=0)
    return matrix.round().astype(int)


def _edgeR(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    if set(groups.unique()) != {"ctrl", "g1"}:
        raise ValueError(f"edgeR requires ctrl and g1; observed {sorted(groups.unique())}.")
    if (groups == "g1").sum() < 2 or (groups == "ctrl").sum() < 2:
        raise ValueError("edgeR requires at least two samples in both g1 and ctrl for this stratum.")

    groups = pd.Categorical(groups, categories=["ctrl", "g1"], ordered=True)
    group_df = pd.DataFrame({"group_simple": groups.astype(str)}, index=counts.columns)
    with (ro.default_converter + pandas2ri.converter).context():
        conv = ro.conversion.get_conversion()
        ro.globalenv["countData"] = conv.py2rpy(counts)
        ro.globalenv["coldata"] = conv.py2rpy(group_df.reset_index(drop=True))

    ro.r("""
        suppressPackageStartupMessages(library(edgeR))
        countData <- as.matrix(countData)
        storage.mode(countData) <- "integer"
        group <- factor(coldata$group_simple, levels=c("ctrl", "g1"))
        y <- DGEList(counts=countData, group=group)
        design <- model.matrix(~group)
        keep <- filterByExpr(y, design=design)
        y <- y[keep, , keep.lib.sizes=FALSE]
        y <- calcNormFactors(y)
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design)
        test <- glmQLFTest(fit, coef=2)
        edgeR_table <- topTags(test, n=Inf)$table
    """)
    with (ro.default_converter + pandas2ri.converter).context():
        out = ro.conversion.get_conversion().rpy2py(ro.globalenv["edgeR_table"])
    out = pd.DataFrame(out)
    out["feature_id"] = out.index.astype(str)
    return out.reset_index(drop=True)


def _fisher(counts: pd.DataFrame, groups: pd.Series) -> pd.DataFrame:
    g1_cols = groups.index[groups.eq("g1")]
    ctrl_cols = groups.index[groups.eq("ctrl")]
    g1 = counts[g1_cols].sum(axis=1).astype(float)
    ctrl = counts[ctrl_cols].sum(axis=1).astype(float)
    total_g1 = float(g1.sum())
    total_ctrl = float(ctrl.sum())
    pvals, odds, logfc = [], [], []
    for a, b in zip(g1.values, ctrl.values):
        table = [[a, max(total_g1 - a, 0)], [b, max(total_ctrl - b, 0)]]
        od, p = fisher_exact(table, alternative="two-sided")
        odds.append(od); pvals.append(p)
        pa = (a + 0.5) / (total_g1 + 1.0)
        pb = (b + 0.5) / (total_ctrl + 1.0)
        logfc.append(np.log2(pa / pb))
    fdr = multipletests(pvals, method="fdr_bh")[1] if pvals else []
    return pd.DataFrame({"feature_id": counts.index.astype(str), "Fisher_log2FC": logfc,
                         "Fisher_odds_ratio": odds, "Fisher_p": pvals, "Fisher_FDR": fdr})


def _v_gene_from_feature(feature: pd.Series) -> pd.Series:
    return feature.astype(str).str.rsplit("|", n=1).str[-1]


def _plot_set_overlap(sets: dict[str, set[str]], stratum: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    venn3([sets["g1"], sets["g5"], sets["g6"]], set_labels=("g1", "g5", "g6"), ax=ax)
    ax.set_title(f"Exact aaV clonotype overlap: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "exact_clonotype_overlap")


def _plot_volcano(features: pd.DataFrame, stratum: str) -> None:
    d = features.replace([np.inf, -np.inf], np.nan).dropna(subset=["logFC", "FDR"]).copy()
    if d.empty: return
    d["minus_log10_fdr"] = -np.log10(d["FDR"].clip(lower=np.finfo(float).tiny))
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(d["logFC"], d["minus_log10_fdr"], s=9, alpha=0.55)
    ax.axvline(0, linewidth=0.8); ax.axhline(-np.log10(0.05), linewidth=0.8, linestyle="--")
    ax.set_xlabel("edgeR log2 fold change (g1 / g5+g6)"); ax.set_ylabel("-log10 FDR")
    ax.set_title(f"Differential clonotype abundance: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "edger_volcano")


def _plot_ranking(ranking: pd.DataFrame, stratum: str) -> None:
    top = ranking.head(15).sort_values(["n_edger_significant", "mean_log2fc"], ascending=True)
    if top.empty: return
    fig, ax = plt.subplots(figsize=(8, 6)); ax.barh(top["v_gene"], top["n_edger_significant"])
    ax.set_xlabel("Significant aaV clonotypes (edgeR FDR < 0.05)"); ax.set_ylabel("V segment")
    ax.set_title(f"Count-based V-segment ranking: {STRATUM_LABELS[stratum]}")
    save_figure(fig, APPROACH, stratum, "v_segment_ranking")


def run_stratum(repertoire: pd.DataFrame, stratum: str) -> dict:
    df = select_stratum(repertoire, stratum)
    ensure_groups(df, ("g1", "g5", "g6"), f"Approach 1 / {stratum}")
    sets = {g: set(df.loc[df["group"].eq(g), "ckey"].unique()) for g in ("g1", "g5", "g6")}
    g1_only = sets["g1"].difference(sets["g5"], sets["g6"])
    set_table = (pd.DataFrame({"feature_id": sorted(g1_only)})
                 .assign(v_gene=lambda x: _v_gene_from_feature(x["feature_id"]))
                 .groupby("v_gene", as_index=False).size()
                 .rename(columns={"size": "n_g1_exclusive"})
                 .sort_values("n_g1_exclusive", ascending=False))

    model_df = df[df["group"].isin(["g1", "g5", "g6"])].copy()
    counts = _count_matrix(model_df)
    meta = sample_metadata(model_df).loc[counts.columns]
    groups = meta["group"].replace({"g5": "ctrl", "g6": "ctrl"}).rename("group_simple")
    edger = _edgeR(counts, groups); fisher = _fisher(counts, groups)
    features = edger.merge(fisher, on="feature_id", how="outer")
    features["v_gene"] = _v_gene_from_feature(features["feature_id"]); features["stratum"] = stratum
    features["consensus_significant"] = features["FDR"].lt(0.05) & features["Fisher_FDR"].lt(0.05)

    sig = features[features["FDR"].lt(0.05)].copy()
    rank = (sig.groupby("v_gene", as_index=False)
            .agg(n_edger_significant=("feature_id", "size"), n_consensus=("consensus_significant", "sum"),
                 mean_log2fc=("logFC", "mean"), max_log2fc=("logFC", "max"), median_fdr=("FDR", "median"))
            .merge(set_table, on="v_gene", how="outer")
            .fillna({"n_edger_significant": 0, "n_consensus": 0, "n_g1_exclusive": 0}))
    rank = rank.sort_values(["n_edger_significant", "mean_log2fc", "n_g1_exclusive"],
                            ascending=[False, False, False], na_position="last").reset_index(drop=True)
    rank["rank"] = np.arange(1, len(rank) + 1); rank["score"] = rank["n_edger_significant"].astype(float); rank["stratum"] = stratum

    out = _output_dir(); features.to_csv(out / f"{stratum}_clonotype_statistics.csv", index=False)
    rank.to_csv(out / f"{stratum}_v_gene_ranking.csv", index=False)
    set_table.to_csv(out / f"{stratum}_g1_exclusive_v_genes.csv", index=False)
    _plot_set_overlap(sets, stratum); _plot_volcano(features, stratum); _plot_ranking(rank, stratum)

    conclusion = {"stratum": stratum, "label": STRATUM_LABELS[stratum],
                  "n_samples": int(model_df["sample_id"].nunique()),
                  "n_g1_exclusive_clonotypes": int(len(g1_only)),
                  "n_edger_significant_clonotypes": int(features["FDR"].lt(0.05).sum()),
                  "n_fisher_edger_consensus": int(features["consensus_significant"].sum()),
                  "top_v_genes": rank.head(5)["v_gene"].dropna().tolist()}
    (out / f"{stratum}_conclusion.json").write_text(json.dumps(conclusion, indent=2)); return conclusion


def compile_summary(strata) -> pd.DataFrame:
    rows = []; out = _output_dir()
    for stratum in strata:
        path = out / f"{stratum}_conclusion.json"
        if path.exists(): rows.append(json.loads(path.read_text()))
    summary = pd.DataFrame(rows); summary.to_csv(out / "stratum_summary.csv", index=False); return summary
