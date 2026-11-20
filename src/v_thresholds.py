"""Empirical upper-right screening for V-segment retention."""

from __future__ import annotations

from statistics import fmean
from typing import NamedTuple

import numpy as np


class PlotThresholds(NamedTuple):
    initial_share_pct: float
    retention_pct: float

    def includes(self, share_pct, retention_pct):
        return (share_pct >= self.initial_share_pct) & (retention_pct >= self.retention_pct)


THRESHOLD_METHOD = "conditional_quantiles"
SHARE_QUANTILE = 0.50
RETENTION_QUANTILE = 0.75
FIXED_THRESHOLDS = PlotThresholds(0.6, 70.0)


def mean_plot_thresholds(fits):
    if not fits:
        return None
    return PlotThresholds(
        fmean(fit.initial_share_pct for fit in fits),
        fmean(fit.retention_pct for fit in fits),
    )


def fit_plot_thresholds(initial_share_pct, retention_pct):
    """Find the median share and the upper retention quartile above that share."""
    share = np.asarray(initial_share_pct, dtype=float)
    retention = np.asarray(retention_pct, dtype=float)
    observed = np.isfinite(share) & np.isfinite(retention) & (share > 0)
    share, retention = share[observed], retention[observed]
    if share.size < 4:
        return None
    share_cutoff = float(np.quantile(share, SHARE_QUANTILE, method="linear"))
    supported = retention[share >= share_cutoff]
    retained = supported[supported > 0]
    if retained.size == 0:
        return None
    retention_cutoff = float(np.quantile(supported, RETENTION_QUANTILE, method="linear"))
    if retention_cutoff <= 0:
        retention_cutoff = float(retained.min())
    return PlotThresholds(share_cutoff, retention_cutoff)


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
    """Fit the displayed cutoffs and audit their sensitivity to mouse omission."""
    from .structural_v import structural_summary

    records = []

    def record(level, group, mouse, cutoffs, status, n_versions, n_fitted, n_selected=None):
        records.append({
            "method": THRESHOLD_METHOD, "level": level,
            "omitted_group": group, "omitted_mouse": mouse, "status": status,
            "n_versions": n_versions, "n_fitted": n_fitted,
            "initial_share_pct": cutoffs.initial_share_pct if cutoffs else None,
            "retention_pct": cutoffs.retention_pct if cutoffs else None,
            "n_selected": n_selected,
        })

    def fit(profile_tables):
        summary = structural_summary(profile_tables)
        share = 100 * summary["initial_share"]
        retention = 100 * summary["retention"]
        cutoffs = fit_plot_thresholds(share, retention)
        return cutoffs, int(cutoffs.includes(share, retention).sum()) if cutoffs else 0

    group_fits = {group: [] for group in DROPOUT_GROUPS}
    n_versions = dict.fromkeys(DROPOUT_GROUPS, 0)
    for group, mouse, training, status in mouse_omissions(tables, metadata, cohort_tables=cohort_tables):
        n_versions[group] += 1
        candidate, n_selected = None, None
        if training is not None:
            candidate, n_selected = fit(training)
            status = "fitted" if candidate else "no_valid_cutoffs"
        if candidate:
            group_fits[group].append(candidate)
        record("mouse", group, mouse, candidate, status, 1, int(bool(candidate)), n_selected)

    for group, fits in group_fits.items():
        mean = mean_plot_thresholds(fits)
        status = "mean_of_omissions" if mean else "no_valid_fit"
        record("group", group, "", mean, status, n_versions[group], len(fits))

    cutoffs, n_selected = fit(tables)
    status = "fitted_on_displayed_profile" if cutoffs else "no_valid_cutoffs"
    record("final", "", "", cutoffs, status, 1, int(bool(cutoffs)), n_selected)
    return cutoffs, records
