"""Exploratory cutoffs for the V-segment retention figures."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


class PlotThresholds(NamedTuple):
    initial_share_pct: float
    retention_pct: float


FIXED_THRESHOLDS = PlotThresholds(0.6, 70.0)


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
