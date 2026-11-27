"""Mouse-level g1 differential abundance in Approach 1."""

from __future__ import annotations

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .count_statistics import pooled_fisher_rows
from .spike_in import aa_v_key, combine_count_tables


def count_result_path(directory, chain, name):
    return directory / (name if chain == "TRA" else f"trb_{name}")


def prepare_counts(tables, metadata, pool_within_mouse):
    counts = combine_count_tables({group: tables[group] for group in ("g1", "g2", "g5", "g6")})

    sample_metadata = (
        metadata[metadata["sample_id"].isin(counts.columns)]
        .drop_duplicates("sample_id")
        .set_index("sample_id")
        .loc[counts.columns]
        .copy()
    )
    sample_metadata["analysis_unit"] = sample_metadata.index.astype(str)
    if pool_within_mouse:
        if sample_metadata["mouse_no"].isna().any():
            raise ValueError("Mouse identifiers are required for pooled inference.")
        sample_metadata["analysis_unit"] = (
            sample_metadata["group_no"].astype(str) + "|" + sample_metadata["mouse_id"].astype(str)
        )
        unit_by_sample = sample_metadata["analysis_unit"].to_dict()
        counts.columns = [unit_by_sample[str(column)] for column in counts.columns]
        counts = counts.T.groupby(level=0, observed=True).sum().T.astype(int)
        sample_metadata = (
            sample_metadata.reset_index()
            .drop_duplicates("analysis_unit")
            .set_index("analysis_unit")
            .loc[counts.columns]
        )
    else:
        sample_metadata = sample_metadata.set_index("analysis_unit", drop=True)

    sample_metadata["contrast_group"] = sample_metadata["group_no"].replace(
        {"g5": "allogeneic", "g6": "allogeneic", "g1": "g1"}
    )
    return counts, sample_metadata


def fit_fisher(counts, sample_metadata, n_workers):
    fisher_samples_g1 = sample_metadata.index[sample_metadata["group_no"].eq("g1")].tolist()
    fisher_samples_ctrl = sample_metadata.index[
        sample_metadata["group_no"].isin(["g5", "g6"])
    ].tolist()

    g1_counts = counts[fisher_samples_g1].sum(axis=1).astype(float)
    control_counts = counts[fisher_samples_ctrl].sum(axis=1).astype(float)
    fisher_rows = pooled_fisher_rows(
        counts.index, g1_counts, control_counts, n_jobs=n_workers
    )

    fisher_results = pd.DataFrame(
        fisher_rows,
        columns=[
            "feature_id",
            "fisher_g1_count",
            "fisher_allogeneic_count",
            "Fisher_log2OR_g1",
            "Fisher_p_g1",
        ],
    )
    if len(fisher_results):
        fisher_results["Fisher_FDR_g1"] = multipletests(
            fisher_results["Fisher_p_g1"], method="fdr_bh"
        )[1]
    else:
        fisher_results["Fisher_FDR_g1"] = pd.Series(dtype=float)
    return fisher_results


def fit_edger(counts, sample_metadata):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri

    edgeR_columns = ["logFC", "logCPM", "F", "PValue", "FDR", "feature_id"]
    edgeR_res_g1 = pd.DataFrame(columns=edgeR_columns)
    status = "not fitted"

    edgeR_units = sample_metadata.index[
        sample_metadata["contrast_group"].isin(["g1", "allogeneic"])
    ].tolist()
    countData_g1 = counts[edgeR_units].copy()
    coldata_g1 = sample_metadata.loc[edgeR_units, ["contrast_group"]].copy()
    replicate_counts = coldata_g1["contrast_group"].value_counts()

    if replicate_counts.get("g1", 0) < 2 or replicate_counts.get("allogeneic", 0) < 2:
        status = (
            "edgeR was not fitted because at least two independent analysis units "
            "were required in both contrast groups."
        )
    else:
        feature_lookup = pd.DataFrame(
            {
                "feature_key": [f"feature_{index:07d}" for index in range(len(countData_g1))],
                "feature_id": list(countData_g1.index),
            }
        )
        countData_g1.index = feature_lookup["feature_key"]
        with (ro.default_converter + pandas2ri.converter).context():
            converter = ro.conversion.get_conversion()
            ro.globalenv["countData"] = converter.py2rpy(countData_g1)
            ro.globalenv["coldata"] = converter.py2rpy(coldata_g1.reset_index(drop=True))

        ro.r(
            """
            suppressPackageStartupMessages(library(edgeR))
            countData <- as.matrix(countData)
            storage.mode(countData) <- "integer"
            group <- factor(coldata$contrast_group, levels=c("allogeneic", "g1"))
            design <- model.matrix(~ group)
            y <- DGEList(counts=countData, group=group)
            y <- calcNormFactors(y)
            keep <- filterByExpr(y, design=design)
            if (any(keep)) {
                y <- y[keep, , keep.lib.sizes=FALSE]
                y <- estimateDisp(y, design)
                fit <- glmQLFit(y, design)
                test <- glmQLFTest(fit, coef=2)
                edgeR_table <- topTags(test, n=Inf)$table
            } else {
                edgeR_table <- data.frame()
            }
            """
        )
        with (ro.default_converter + pandas2ri.converter).context():
            converted = ro.conversion.get_conversion().rpy2py(ro.globalenv["edgeR_table"])
        edgeR_res_g1 = pd.DataFrame(converted)
        if not edgeR_res_g1.empty:
            edgeR_res_g1["feature_key"] = edgeR_res_g1.index.astype(str)
            edgeR_res_g1 = edgeR_res_g1.merge(feature_lookup, on="feature_key", how="left")
            edgeR_res_g1 = edgeR_res_g1.drop(columns="feature_key")
            status = "edgeR quasi-likelihood model fitted successfully."
        else:
            edgeR_res_g1 = pd.DataFrame(columns=edgeR_columns)
            status = "No aaV clonotype passed edgeR expression filtering."
    return edgeR_res_g1, status


def rank_v_segments(remaining, edger, fisher_results, stratum, stratum_label):
    edger_plot = edger.copy()
    edger_plot["v_gene"] = edger_plot["feature_id"].map(lambda value: aa_v_key(value)[1])

    if not edger_plot.empty:
        edger_v = edger_plot.groupby("v_gene", as_index=False).agg(
            n_edger_tested=("feature_id", "count"),
            min_edger_fdr=("FDR", "min"),
        )
        edgeR_significant = (
            edger_plot[edger_plot["FDR"].lt(0.05)]
            .groupby("v_gene", as_index=False)
            .agg(
                n_edger_significant=("feature_id", "count"),
                mean_edger_log2fc=("logFC", "mean"),
                max_edger_log2fc=("logFC", "max"),
            )
        )
        edger_v = edger_v.merge(edgeR_significant, on="v_gene", how="left")
    else:
        edger_v = pd.DataFrame(
            columns=[
                "v_gene", "n_edger_tested", "n_edger_significant",
                "mean_edger_log2fc", "max_edger_log2fc", "min_edger_fdr",
            ]
        )

    exclusive = pd.DataFrame(sorted(remaining), columns=["cdr3aa", "v_gene"])
    exclusive_v = (
        exclusive.groupby("v_gene", as_index=False)
        .size()
        .rename(columns={"size": "n_g1_exclusive_aav"})
    )

    fisher_parsed = fisher_results.copy()
    if not fisher_parsed.empty:
        fisher_parsed["v_gene"] = fisher_parsed["feature_id"].map(lambda value: aa_v_key(value)[1])
        fisher_v = fisher_parsed.groupby("v_gene", as_index=False).agg(
            min_fisher_fdr=("Fisher_FDR_g1", "min"),
        )
        fisher_significant = (
            fisher_parsed[fisher_parsed["Fisher_FDR_g1"].lt(0.05)]
            .groupby("v_gene", as_index=False)
            .agg(
                n_fisher_significant=("feature_id", "count"),
                mean_fisher_log2or=("Fisher_log2OR_g1", "mean"),
            )
        )
        fisher_v = fisher_v.merge(fisher_significant, on="v_gene", how="left")
    else:
        fisher_v = pd.DataFrame(
            columns=["v_gene", "n_fisher_significant", "mean_fisher_log2or", "min_fisher_fdr"]
        )

    ranking = exclusive_v.merge(edger_v, on="v_gene", how="outer").merge(
        fisher_v, on="v_gene", how="outer"
    )
    count_columns = [
        "n_g1_exclusive_aav", "n_edger_tested", "n_edger_significant", "n_fisher_significant"
    ]
    for column in count_columns:
        ranking[column] = pd.to_numeric(ranking[column], errors="coerce").fillna(0).astype(int)

    ranking["edgeR_g1_enriched"] = (
        ranking["n_edger_significant"].gt(0)
        & ranking["mean_edger_log2fc"].gt(0)
    )
    ranking["fisher_g1_enriched"] = (
        ranking["n_fisher_significant"].gt(0)
        & ranking["mean_fisher_log2or"].gt(0)
    )
    ranking["evidence_class"] = np.select(
        [
            ranking["edgeR_g1_enriched"] & ranking["fisher_g1_enriched"],
            ranking["edgeR_g1_enriched"],
            ranking["fisher_g1_enriched"],
            ranking["n_g1_exclusive_aav"].gt(0),
        ],
        [
            "edgeR and Fisher concordant",
            "edgeR-supported",
            "Fisher-supported",
            "exact-set descriptive",
        ],
        default="tested without positive enrichment evidence",
    )
    ranking["stratum"] = stratum
    ranking["stratum_label"] = stratum_label
    ranking = ranking.sort_values(
        ["n_g1_exclusive_aav", "v_gene"],
        ascending=[False, True],
    ).reset_index(drop=True)
    ranking["stratum_rank"] = np.arange(1, len(ranking) + 1)
    return ranking
