"""Count-based statistics used by Approach 1."""

from __future__ import annotations

import math


def historical_fisher_log2_ratio(
    g1_count: float,
    allogeneic_count: float,
) -> float:
    """Return the pooled-count log2 ratio with the historical pseudocount."""
    if g1_count < 0 or allogeneic_count < 0:
        raise ValueError("Fisher pooled counts must be non-negative.")
    return math.log2((float(g1_count) + 1.0) / (float(allogeneic_count) + 1.0))


def _fisher_batch(task):
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
        results.append(
            {
                "feature_id": feature_id,
                "fisher_g1_count": g1_count,
                "fisher_allogeneic_count": control_count,
                "Fisher_log2FC_g1": np.log2(
                    (g1_count + 1) / (control_count + 1)
                ),
                "Fisher_odds_ratio": odds_ratio,
                "Fisher_p_g1": p_value,
            }
        )
    return results


def pooled_fisher_rows(
    feature_ids,
    g1_counts,
    control_counts,
    *,
    n_jobs=None,
    batch_size=512,
):
    """Run Fisher tests in ordered batches; BH correction is applied later."""
    from src.runtime import parallel_map

    if batch_size < 1:
        raise ValueError("batch_size must be positive")

    pairs = list(zip(g1_counts, control_counts, strict=True))
    # With fixed library totals, identical pooled counts define the same 2×2 table.
    rows = [(pair, *pair) for pair in dict.fromkeys(pairs)]
    total_g1 = float(g1_counts.sum())
    total_control = float(control_counts.sum())
    tasks = [
        (rows[start : start + batch_size], total_g1, total_control)
        for start in range(0, len(rows), batch_size)
    ]
    batches = parallel_map(
        _fisher_batch,
        tasks,
        n_jobs=n_jobs,
        prefer="processes",
    )
    by_pair = {row["feature_id"]: row for batch in batches for row in batch}
    return [
        {**by_pair[pair], "feature_id": feature_id}
        for feature_id, pair in zip(feature_ids, pairs, strict=True)
    ]


def _two_sided_fisher_pvalues(a, b, total_g1, total_control):
    """Sum hypergeometric probabilities no larger than the observed probability."""
    import numpy as np
    from scipy.stats import hypergeom

    total = total_g1 + total_control
    draws = a + b
    mode = (draws + 1) * (total_g1 + 1) // (total + 2)
    # Swap columns for observations above the mode. This makes every search
    # start with the lower tail and keeps the same two-sided probability.
    swap = a > mode
    good = np.where(swap, total_control, total_g1)
    observed = np.where(swap, b, a)
    mode = (draws + 1) * (good + 1) // (total + 2)
    observed_mass = hypergeom.pmf(observed, total, good, draws)
    mode_mass = hypergeom.pmf(mode, total, good, draws)
    p = np.ones(len(draws))
    off_mode = (mode_mass - observed_mass) > 1e-14 * mode_mass
    k, g, n = observed[off_mode], good[off_mode], draws[off_mode]
    cutoff = observed_mass[off_mode] * (1 + 1e-14)
    left = mode[off_mode] + 1
    right = np.minimum(n, g) + 1
    while np.any(left < right):
        active = left < right
        middle = (left[active] + right[active]) // 2
        within = hypergeom.pmf(middle, total, g[active], n[active]) <= cutoff[active]
        left[active] = np.where(within, left[active], middle + 1)
        right[active] = np.where(within, middle, right[active])
    p[off_mode] = hypergeom.cdf(k, total, g, n) + hypergeom.sf(left - 1, total, g, n)
    return p.clip(max=1)


def pooled_fisher_arrays(g1_counts, control_counts):
    """Evaluate the two-sided Fisher test once per distinct pooled count pair."""
    import numpy as np

    total_g1, total_control = int(g1_counts.sum()), int(control_counts.sum())
    base = int(control_counts.max(initial=0)) + 1
    keys = g1_counts.astype(np.int64) * base + control_counts.astype(np.int64)
    unique_keys, inverse = np.unique(keys, return_inverse=True)
    a, b = np.divmod(unique_keys, base)
    odds = np.full(len(a), np.inf)
    denominator = b.astype(float) * (total_g1 - a)
    np.divide(a.astype(float) * (total_control - b), denominator, out=odds, where=denominator > 0)
    degenerate = (a + b == 0) | (a + b == total_g1 + total_control) | (total_g1 == 0) | (total_control == 0)
    odds[degenerate] = np.nan
    probabilities = np.ones(len(a))
    relevant = ~degenerate
    probabilities[relevant] = _two_sided_fisher_pvalues(a[relevant], b[relevant], total_g1, total_control)
    return odds[inverse], probabilities[inverse]
