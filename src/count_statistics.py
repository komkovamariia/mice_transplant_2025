"""Dependency-light statistical definitions shared by count-based analyses."""

from __future__ import annotations

import math


def historical_fisher_log2_ratio(
    g1_count: float,
    allogeneic_count: float,
) -> float:
    """Return the pseudocount-adjusted pooled-count ratio used historically."""
    if g1_count < 0 or allogeneic_count < 0:
        raise ValueError("Fisher pooled counts must be non-negative.")
    return math.log2((float(g1_count) + 1.0) / (float(allogeneic_count) + 1.0))
