"""Exact aaV structural and overlap profiles for V-segment figures."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .spike_in import SUBTRACTION_GROUPS, aa_v_key


def v_region_proportions(tables):
    counts = tables["g1"].sum(axis=1)
    features = [aa_v_key(key) for key in counts.index]
    in_g5 = counts.index.isin(tables["g5"].index)
    in_g6 = counts.index.isin(tables["g6"].index)
    regions = np.select([in_g5 & in_g6, in_g5, in_g6], ["111", "110", "101"], default="100")
    burden = pd.DataFrame({"v": [v for _, v in features], "region": regions, "umi": counts.to_numpy()})
    wide = (
        burden.groupby(["v", "region"])["umi"].sum().unstack(fill_value=0)
        .reindex(columns=["100", "110", "101", "111"], fill_value=0)
    )
    return wide.div(wide.sum(axis=1), axis=0).fillna(0)


def structural_remainder(tables):
    return set(tables["g1"].index).difference(
        *(tables[group].index for group in SUBTRACTION_GROUPS)
    )


def structural_summary(tables):
    g1 = set(tables["g1"].index)
    remaining = structural_remainder(tables)
    n_total = len(g1)

    full = pd.Series([aa_v_key(feature)[1] for feature in g1], dtype="object").value_counts()
    retained = pd.Series(
        [aa_v_key(feature)[1] for feature in remaining],
        dtype="object",
    ).value_counts()

    summary = full.rename("n_full").to_frame()
    summary = summary.join(retained.rename("n_remaining"), how="left").fillna(0)
    summary["n_remaining"] = summary["n_remaining"].astype(int)
    summary["initial_share"] = summary["n_full"] / n_total
    summary["retention"] = summary["n_remaining"] / summary["n_full"]
    summary["score_s"] = (
        summary["initial_share"] * summary["retention"]
    ).round(12)

    summary["v"] = summary.index.astype(str)
    return summary.reset_index(drop=True)


def v_plot_profile(tables):
    summary = structural_summary(tables).rename(columns={
        "n_full": "cdr3_in_full_g1",
        "n_remaining": "cdr3_remaining_in_g1",
    })
    summary["frequency_in_full_g1_pct"] = 100 * summary["initial_share"]
    summary["retained_cdr3_share_pct"] = 100 * summary["retention"]
    summary["frequency_remaining_in_g1_pct_relative_to_original"] = (
        100 * summary["cdr3_remaining_in_g1"] / summary["cdr3_in_full_g1"].sum()
    )
    profile = summary.join(v_region_proportions(tables), on="v")
    probabilities = profile[["100", "110", "101", "111"]].replace(0, np.nan)
    profile["region_entropy"] = -(probabilities * np.log2(probabilities)).sum(axis=1) / np.log2(4)
    return profile


def subset_v_tables(tables, sample_ids):
    selected = {}
    for group, table in tables.items():
        subset = table.loc[:, table.columns.isin(sample_ids)]
        selected[group] = subset.loc[subset.sum(axis=1).gt(0)]
    return selected
