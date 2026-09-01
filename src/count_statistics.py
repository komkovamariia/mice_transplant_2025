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


def _fisher_batch(task):
    """Evaluate a batch without serializing a notebook's global namespace."""
    import numpy as np
    from scipy.stats import fisher_exact

    rows, total_g1, total_control = task
    results = []
    for feature_id, g1_count, control_count in rows:
        table = [
            [g1_count, control_count],
            [max(total_g1 - g1_count, 0), max(total_control - control_count, 0)],
        ]
        odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
        results.append({
            "feature_id": feature_id,
            "fisher_g1_count": g1_count,
            "fisher_allogeneic_count": control_count,
            "Fisher_log2FC_g1": np.log2((g1_count + 1) / (control_count + 1)),
            "Fisher_odds_ratio": odds_ratio,
            "Fisher_p_g1": p_value,
        })
    return results


def pooled_fisher_rows(feature_ids, g1_counts, control_counts, *, n_jobs=None, batch_size=512):
    """Exact historical tests in ordered batches; BH correction stays global."""
    from src.runtime import parallel_map

    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    rows = list(zip(feature_ids, g1_counts, control_counts, strict=True))
    total_g1 = float(g1_counts.sum())
    total_control = float(control_counts.sum())
    tasks = [(rows[i:i + batch_size], total_g1, total_control)
             for i in range(0, len(rows), batch_size)]
    batches = parallel_map(_fisher_batch, tasks, n_jobs=n_jobs, prefer="processes")
    return [row for batch in batches for row in batch]
