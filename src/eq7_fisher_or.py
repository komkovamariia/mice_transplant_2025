"""Sensitivity of aaV-level Fisher enrichment to the direction estimator.

Fisher exact p-values and BH FDR are unchanged. Only the sign of enrichment
switches from raw pooled UMI ratio to the 2x2 table's log odds ratio.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .differential_abundance import fit_fisher
from .eq7_aav_fraction import summarize_positive_fraction
from .spike_in import aa_v_key

ODDS_FISHER = "Fisher_odds_positive_fraction"
RAW_FISHER = "Fisher_positive_fraction"


def add_log2_odds_ratio(fisher: pd.DataFrame) -> pd.DataFrame:
    """Haldane-Anscombe corrected odds ratio, accounting for group totals."""
    out = fisher.copy()
    a = out["fisher_g1_count"].to_numpy(dtype=float)
    b = out["fisher_allogeneic_count"].to_numpy(dtype=float)
    n_g1, n_allogeneic = a.sum(), b.sum()
    if n_g1 <= 0 or n_allogeneic <= 0:
        raise ValueError("Both groups must have positive total functional UMI counts.")
    if (a < 0).any() or (b < 0).any():
        raise ValueError("Negative UMI counts.")
    out["Fisher_log2OR_g1"] = np.log2(
        ((a + 0.5) * (n_allogeneic - b + 0.5))
        / ((b + 0.5) * (n_g1 - a + 0.5))
    )
    return out


def fisher_fraction_sensitivity(counts, units, excluded=(), *, workers=4):
    """Recalculate ONLY Fisher on the retained mice of an existing split.

    No edgeR, rarefaction, mouse split generation or subtraction is performed.
    """
    keys = pd.MultiIndex.from_frame(units[["group_no", "mouse_id"]])
    remaining = units.loc[~keys.isin(excluded)]
    contrast_cols = remaining.index[
        remaining["contrast_group"].isin(("g1", "allogeneic"))
    ]
    train = counts.loc[:, remaining.index]
    train = train.loc[train[contrast_cols].sum(axis=1).gt(0)]
    genes = pd.Index(sorted({aa_v_key(k)[1] for k in train.index}), name="v_gene")
    fisher = add_log2_odds_ratio(fit_fisher(train, remaining, n_workers=workers))
    raw_frac, raw_pos, n_tested = summarize_positive_fraction(
        fisher, "Fisher_FDR_g1", "Fisher_log2FC_g1", genes,
    )
    odds_frac, odds_pos, odds_tested = summarize_positive_fraction(
        fisher, "Fisher_FDR_g1", "Fisher_log2OR_g1", genes,
    )
    np.testing.assert_array_equal(n_tested, odds_tested)
    return pd.DataFrame({
        "v_gene": genes,
        RAW_FISHER: raw_frac,
        ODDS_FISHER: odds_frac,
        "n_Fisher_raw_positive": raw_pos,
        "n_Fisher_odds_positive": odds_pos,
        "n_Fisher_tested": n_tested,
    })
