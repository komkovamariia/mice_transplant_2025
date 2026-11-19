"""Fit exploratory Equation 7 weights from leave-one-mouse-out evidence."""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

from .differential_abundance import fit_edger, fit_fisher, rank_v_segments
from .spike_in import EDGER_SUCCESS, aa_v_key
from .structural_v import structural_remainder, v_region_proportions


FEATURES = ("edgeR_g1_enriched", "fisher_g1_enriched", "log1p_n_g1_exclusive_aav")


def shortlist(ranking, top):
    return ranking.sort_values(
        ["only_g1_share", "n_full", "v_gene"], ascending=[False, False, True],
    ).head(top)


def with_structural_fields(ranking, tables, *, remainder=None):
    if remainder is None:
        remainder = structural_remainder(tables)
    regions = v_region_proportions(tables)
    full = Counter(aa_v_key(key)[1] for key in tables["g1"].index)
    remaining = Counter(aa_v_key(key)[1] for key in remainder)
    rows = pd.DataFrame({"v_gene": list(full)}).merge(
        ranking[["v_gene", *FEATURES[:2]]], on="v_gene", how="left", validate="1:1",
    )
    for term in FEATURES[:2]:
        rows[term] = rows[term].eq(True).astype(int)
    rows["only_g1_share"] = rows["v_gene"].map(regions["100"])
    rows["n_full"] = rows["v_gene"].map(full)
    rows["n_g1_exclusive_aav"] = rows["v_gene"].map(remaining).fillna(0).astype(int)
    rows["log1p_n_g1_exclusive_aav"] = np.log1p(rows["n_g1_exclusive_aav"])
    return rows


def fold_evidence(training, held_mouse_table, group, mouse, counts, units, *, n_workers):
    """Compute V-level training evidence and the omitted mouse's response."""
    train_units = units.loc[~(units["mouse_id"].eq(mouse) & units["group_no"].eq(group))]
    train_counts = counts.loc[:, train_units.index]
    train_counts = train_counts.loc[train_counts.sum(axis=1).gt(0)]
    edger, status = fit_edger(train_counts, train_units)
    if status != EDGER_SUCCESS:
        return None, status
    fisher = fit_fisher(train_counts, train_units, n_workers)
    g1_only = set(training["g1"].index) - set(training["g5"].index) - set(training["g6"].index)
    remainder = structural_remainder(training)
    ranking = rank_v_segments(remainder, edger, fisher, "all_combined", "Mouse-omitted training set")
    rows = with_structural_fields(ranking, training, remainder=remainder)

    observed = set(held_mouse_table.index)
    available = Counter(aa_v_key(key)[1] for key in observed)
    recovered = Counter(aa_v_key(key)[1] for key in observed & g1_only)
    rows["heldout_n_v"] = rows["v_gene"].map(available).fillna(0).astype(int)
    rows["heldout_n_recovered"] = rows["v_gene"].map(recovered).fillna(0).astype(int)
    sign = 1 if group == "g1" else -1
    rows["response"] = sign * rows["heldout_n_recovered"] / rows["heldout_n_v"].replace(0, np.nan)
    rows["heldout_mouse"] = mouse
    rows["heldout_group"] = group
    return rows, "fitted"


def calibrate_weights(rows, top, ridge):
    """Estimate nonnegative weights using squared error on signed held-out recovery."""
    mouse_key = ["heldout_group", "heldout_mouse"]
    selected = pd.concat([shortlist(fold, top) for _, fold in rows.groupby(mouse_key)], ignore_index=True)
    selected["used_for_fit"] = selected["response"].notna()
    fitted = selected.loc[selected["used_for_fit"]].copy()
    if set(fitted["heldout_group"]) != {"g1", "g5", "g6"}:
        raise ValueError("Weight fitting requires usable omitted-mouse data from g1, g5 and g6.")
    n_mice = fitted.groupby("heldout_group")["heldout_mouse"].transform("nunique")
    n_segments = fitted.groupby(mouse_key)["v_gene"].transform("size")
    fitted["regression_weight"] = 1 / (3 * n_mice * n_segments)
    selected["regression_weight"] = fitted["regression_weight"].reindex(selected.index, fill_value=0.0)
    weights = fitted["regression_weight"].to_numpy()
    x = fitted[list(FEATURES)].to_numpy(dtype=float)
    y = fitted["response"].to_numpy(dtype=float)
    if np.ptp(y) == 0:
        raise ValueError("Equation 7 calibration has no variation in held-out recovery.")
    mean_x = np.average(x, axis=0, weights=weights)
    mean_y = np.average(y, weights=weights)
    std_x = np.sqrt(np.average((x - mean_x)**2, axis=0, weights=weights))
    active = std_x > 1e-12
    coefficients = np.zeros(len(FEATURES))
    if not active.any():
        raise ValueError("Equation 7 inputs are constant across all selected fold rows.")
    standardized = (x[:, active] - mean_x[active]) / std_x[active]
    design = np.vstack([standardized * np.sqrt(weights[:, None]), np.sqrt(ridge) * np.eye(active.sum())])
    response = np.concatenate([(y - mean_y) * np.sqrt(weights), np.zeros(active.sum())])
    solution = lsq_linear(design, response, bounds=(0, np.inf))
    if not solution.success:
        raise RuntimeError(f"Equation 7 weight fit did not converge: {solution.message}")
    coefficients[active] = solution.x / std_x[active]
    fit = pd.DataFrame({"term": FEATURES, "weight": coefficients, "varied_in_training": active})
    return fit, selected


def score_candidates(ranking, fit, top):
    rows = shortlist(ranking, top).copy()
    weights = fit.set_index("term")["weight"]
    contributions = [f"contribution_{term}" for term in FEATURES]
    for term, contribution in zip(FEATURES, contributions, strict=True):
        rows[contribution] = rows[term] * weights[term]
    rows["calibrated_score"] = rows[contributions].sum(axis=1)
    return rows.sort_values(["calibrated_score", "only_g1_share", "v_gene"], ascending=[False, False, True])
