"""Exploratory cutoffs for the V-segment retention figures."""

from __future__ import annotations

from statistics import fmean
from typing import NamedTuple

import numpy as np


class PlotThresholds(NamedTuple):
    initial_share_pct: float
    retention_pct: float


FIXED_THRESHOLDS = PlotThresholds(0.6, 70.0)


def mean_plot_thresholds(fits):
    if not fits:
        return None
    return PlotThresholds(
        fmean(fit.initial_share_pct for fit in fits),
        fmean(fit.retention_pct for fit in fits),
    )


def _split(values, min_upper=2):
    # Use an exact split between observed values, with two segments on each side.
    ordered = np.sort(values)
    n = len(ordered)
    boundaries = np.arange(2, n - min_upper + 1)
    boundaries = boundaries[ordered[boundaries - 1] < ordered[boundaries]]
    if not len(boundaries):
        return None
    cumulative = np.cumsum(ordered)
    lower = cumulative[boundaries - 1]
    upper = cumulative[-1] - lower
    difference = upper / (n - boundaries) - lower / boundaries
    variance = boundaries * (n - boundaries) * difference**2
    split = boundaries[np.argmax(variance)]
    return float((ordered[split - 1] + ordered[split]) / 2)


def fit_plot_thresholds(initial_share_pct, retention_pct):
    """Fit the share split on a log scale, then split retention above it."""
    share = np.asarray(initial_share_pct, dtype=float)
    retention = np.asarray(retention_pct, dtype=float)
    observed = np.isfinite(share) & np.isfinite(retention) & (share > 0)
    share, retention = share[observed], retention[observed]
    log_share = np.log10(share)
    share_split = _split(log_share, min_upper=4)
    if share_split is None:
        return None
    retention_split = _split(retention[log_share >= share_split])
    if retention_split is None:
        return None
    return PlotThresholds(10**share_split, retention_split)


DROPOUT_GROUPS = ("g1", "g5", "g6")


def mouse_omissions(tables, metadata, *, cohort_tables=None):
    """Omit all libraries of one mouse from one group at a time."""
    mouse_by_sample = metadata.set_index("sample_id")["mouse_id"]
    cohort = tables if cohort_tables is None else cohort_tables
    for group in DROPOUT_GROUPS:
        table = tables[group]
        sample_mice = mouse_by_sample.loc[table.columns]
        cohort_mice = sorted(mouse_by_sample.loc[cohort[group].columns].unique())
        present = set(sample_mice)
        for mouse in cohort_mice:
            if mouse not in present:
                yield group, mouse, None, "no_panel_libraries"
                continue
            columns = sample_mice.index[sample_mice.ne(mouse)]
            if not len(columns):
                yield group, mouse, None, "skipped_empty_group"
                continue
            filtered = table.loc[:, columns]
            held_out = {**tables, group: filtered.loc[filtered.sum(axis=1).gt(0)]}
            yield group, mouse, held_out, "omitted"


def fit_mouse_dropout_thresholds(tables, metadata, *, cohort_tables=None):
    from .structural_v import structural_summary

    records = []

    def record(level, group, mouse, cutoffs, status, n_versions, n_fitted):
        records.append({
            "level": level, "omitted_group": group, "omitted_mouse": mouse, "status": status,
            "n_versions": n_versions, "n_fitted": n_fitted,
            "initial_share_pct": cutoffs.initial_share_pct if cutoffs else None,
            "retention_pct": cutoffs.retention_pct if cutoffs else None,
        })

    group_fits = {group: [] for group in DROPOUT_GROUPS}
    n_versions = dict.fromkeys(DROPOUT_GROUPS, 0)
    for group, mouse, training, status in mouse_omissions(tables, metadata, cohort_tables=cohort_tables):
        n_versions[group] += 1
        candidate = None
        if training is not None:
            if training["g1"].empty:
                status = "empty_target"
            else:
                summary = structural_summary(training)
                candidate = fit_plot_thresholds(100 * summary["initial_share"], 100 * summary["retention"])
                status = "fitted" if candidate else "no_split"
        if candidate:
            group_fits[group].append(candidate)
        record("mouse", group, mouse, candidate, status, 1, int(bool(candidate)))

    group_means = []
    for group, fits in group_fits.items():
        mean = mean_plot_thresholds(fits)
        if mean:
            group_means.append(mean)
        record("group", group, "", mean, "mean" if mean else "no_valid_fit", n_versions[group], len(fits))

    cutoffs = mean_plot_thresholds(group_means) if len(group_means) == len(DROPOUT_GROUPS) else None
    status = "mean_of_group_means" if cutoffs else "unavailable_missing_group_mean"
    record("final", "", "", cutoffs, status, sum(n_versions.values()), sum(map(len, group_fits.values())))
    return cutoffs, records
