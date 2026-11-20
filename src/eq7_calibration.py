"""Calibrate Equation 7 on training-fold V segments and omitted mice."""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from .differential_abundance import fit_edger, fit_fisher, rank_v_segments
from .spike_in import EDGER_SUCCESS, aa_v_key
from .structural_v import structural_remainder, v_region_proportions
from .v_thresholds import DROPOUT_GROUPS, THRESHOLD_METHOD, fit_plot_thresholds


FEATURES = ("edgeR_g1_enriched", "fisher_g1_enriched", "log1p_n_g1_exclusive_aav")
STANDARDIZED_FEATURES = tuple(f"z_{feature}" for feature in FEATURES)
CONTRIBUTIONS = tuple(f"contribution_{feature}" for feature in FEATURES)
RIDGE_GRID = (0.001, 0.01, 0.1, 1.0, 10.0)
DEFAULT_RIDGE = 0.1
RESPONSE_COLUMNS = ["heldout_n_v", "heldout_n_recovered", "response"]
MOUSE_COLUMNS = ["heldout_group", "heldout_mouse"]


class CalibrationUnavailable(ValueError):
    """The training folds do not supply enough evidence to fit Equation 7."""


def profile_coordinates(ranking):
    rows = ranking.copy()
    rows["initial_share_pct"] = 100 * (rows["n_full"] / rows["n_full"].sum())
    rows["retention_pct"] = 100 * (rows["n_g1_exclusive_aav"] / rows["n_full"])
    return rows


def shortlist(ranking, top=None, *, frozen_cutoffs=None):
    rows = profile_coordinates(ranking)
    cutoff = frozen_cutoffs
    if cutoff is None:
        cutoff = fit_plot_thresholds(rows["initial_share_pct"], rows["retention_pct"])
    rows["threshold_method"] = THRESHOLD_METHOD if frozen_cutoffs is None else "frozen_" + THRESHOLD_METHOD
    rows["training_share_cutoff_pct"] = cutoff.initial_share_pct if cutoff else np.nan
    rows["training_retention_cutoff_pct"] = cutoff.retention_pct if cutoff else np.nan
    if cutoff is None:
        return rows.iloc[:0]
    eligible = rows.loc[cutoff.includes(rows["initial_share_pct"], rows["retention_pct"])].sort_values(
        ["only_g1_share", "n_full", "v_gene"], ascending=[False, False, True],
    )
    return eligible if top is None else eligible.head(top)


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


def training_evidence(training, excluded, counts, units, *, n_workers):
    """Fit evidence with all libraries of each omitted mouse excluded."""
    unit_keys = pd.MultiIndex.from_frame(units[["group_no", "mouse_id"]])
    train_units = units.loc[~unit_keys.isin(excluded)]
    train_counts = counts.loc[:, train_units.index]
    contrast_units = train_units.index[train_units["contrast_group"].isin(["g1", "allogeneic"])]
    train_counts = train_counts.loc[train_counts.loc[:, contrast_units].sum(axis=1).gt(0)]
    edger, status = fit_edger(train_counts, train_units)
    if status != EDGER_SUCCESS:
        return None, status
    fisher = fit_fisher(train_counts, train_units, n_workers)
    remainder = structural_remainder(training)
    ranking = rank_v_segments(remainder, edger, fisher, "all_combined", "Mouse-omitted training set")
    return with_structural_fields(ranking, training, remainder=remainder), "fitted"


def attach_mouse_response(evidence, training, held_mouse_table, group, mouse):
    rows = evidence.copy()
    g1_only = set(training["g1"].index).difference(training["g5"].index, training["g6"].index)
    observed = set(held_mouse_table.index)
    available = Counter(aa_v_key(key)[1] for key in observed)
    recovered = Counter(aa_v_key(key)[1] for key in observed & g1_only)
    rows["heldout_n_v"] = rows["v_gene"].map(available).fillna(0).astype(int)
    rows["heldout_n_recovered"] = rows["v_gene"].map(recovered).fillna(0).astype(int)
    sign = 1 if group == "g1" else -1
    rows["response"] = sign * rows["heldout_n_recovered"] / rows["heldout_n_v"].replace(0, np.nan)
    rows["heldout_mouse"] = mouse
    rows["heldout_group"] = group
    return rows


def select_calibration_rows(rows, top):
    """Give groups, mice and eligible folds equal weight within their parent units."""
    if rows.empty:
        raise CalibrationUnavailable("No usable inner evidence fits.")
    fold_key = [*MOUSE_COLUMNS, "evidence_split"] if "evidence_split" in rows else MOUSE_COLUMNS
    selected = pd.concat([shortlist(fold, top) for _, fold in rows.groupby(fold_key)], ignore_index=True)
    selected["used_for_fit"] = np.isfinite(selected["response"])
    fitted = selected.loc[selected["used_for_fit"]].copy()
    missing = set(DROPOUT_GROUPS).difference(fitted["heldout_group"])
    if missing:
        raise CalibrationUnavailable("No eligible training responses for " + ", ".join(sorted(missing)) + ".")
    n_mice = fitted.groupby("heldout_group")["heldout_mouse"].transform("nunique")
    n_folds = fitted.groupby(MOUSE_COLUMNS)["evidence_split"].transform("nunique") if "evidence_split" in fitted else 1
    n_segments = fitted.groupby(fold_key)["v_gene"].transform("size")
    fitted["regression_weight"] = 1 / (len(DROPOUT_GROUPS) * n_mice * n_folds * n_segments)
    selected["regression_weight"] = fitted["regression_weight"].reindex(selected.index, fill_value=0.0)
    return selected


def fit_selected(selected, ridge):
    """Estimate nonnegative weighted ridge slopes on the standardized inputs."""
    selected = selected.copy()
    fitted = selected.loc[selected["used_for_fit"]]
    weights = fitted["regression_weight"].to_numpy()
    x = fitted[list(FEATURES)].to_numpy(dtype=float)
    y = fitted["response"].to_numpy(dtype=float)
    mean_x = np.average(x, axis=0, weights=weights)
    mean_y = np.average(y, weights=weights)
    std_x = np.sqrt(np.average((x - mean_x)**2, axis=0, weights=weights))
    active = std_x > 1e-12
    coefficients = np.zeros(len(FEATURES))
    standardized_weights = coefficients.copy()
    if active.any() and np.ptp(y) > 0:
        standardized = (x[:, active] - mean_x[active]) / std_x[active]
        design = np.vstack([standardized * np.sqrt(weights[:, None]), np.sqrt(ridge) * np.eye(active.sum())])
        response = np.concatenate([(y - mean_y) * np.sqrt(weights), np.zeros(active.sum())])
        standardized_weights[active], _ = nnls(design, response)
        coefficients[active] = standardized_weights[active] / std_x[active]
    intercept = mean_y - mean_x @ coefficients
    fit = pd.DataFrame({
        "term": FEATURES, "weight": coefficients, "varied_in_training": active,
        "training_mean": mean_x, "training_sd": std_x,
        "standardized_weight": standardized_weights, "intercept": intercept,
        "baseline_prediction": mean_y,
    })
    selected[list(STANDARDIZED_FEATURES)] = (
        (selected[list(FEATURES)].to_numpy() - mean_x) / np.where(active, std_x, 1)
    ) * active
    selected["E"] = selected[list(FEATURES)].to_numpy() @ coefficients
    selected["fitted_response"] = intercept + selected["E"]
    selected["residual"] = selected["response"] - selected["fitted_response"]
    return fit, selected


def calibrate_weights(rows, top, ridge):
    return fit_selected(select_calibration_rows(rows, top), ridge)


def score_candidates(ranking, fit, top=None, *, frozen_cutoffs=None):
    rows = shortlist(ranking, top, frozen_cutoffs=frozen_cutoffs).copy()
    weights = fit.set_index("term")["weight"]
    for term, contribution in zip(FEATURES, CONTRIBUTIONS, strict=True):
        rows[contribution] = rows[term] * weights.get(term, np.nan)
    rows["E"] = rows[list(CONTRIBUTIONS)].sum(axis=1, min_count=len(FEATURES))
    rows["calibrated_score"] = rows["E"]
    rows["predicted_response"] = rows["E"] + (fit["intercept"].iat[0] if not fit.empty else np.nan)
    return rows.sort_values(["E", "only_g1_share", "v_gene"], ascending=[False, False, True])
