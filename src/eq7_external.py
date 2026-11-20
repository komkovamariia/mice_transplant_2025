"""Transport the frozen TRB score to PRJNA954716 without refitting it."""

import numpy as np
import pandas as pd

from .differential_abundance import prepare_counts
from .eq7_calibration import (
    MOUSE_COLUMNS, RESPONSE_COLUMNS, CalibrationUnavailable, attach_mouse_response, profile_coordinates,
    score_candidates, training_evidence,
)
from .spike_in import SUBTRACTION_GROUPS
from .structural_sampling import downsample_group_tables
from .structural_v import subset_v_tables
from .v_thresholds import fit_plot_thresholds


MATCHED_GROUPS = ("g1", "g5")
EMPTY_BACKGROUND = tuple(group for group in SUBTRACTION_GROUPS if group not in MATCHED_GROUPS)


def matched_panel(raw, metadata, *, depth, seed, n_workers):
    panel = subset_v_tables({group: raw[group] for group in MATCHED_GROUPS}, metadata["sample_id"])
    n_mice = metadata.loc[metadata["sample_id"].isin(
        column for table in panel.values() for column in table.columns
    )].groupby("group_no")["mouse_id"].nunique()
    if n_mice.reindex(MATCHED_GROUPS, fill_value=0).min() < 2:
        raise CalibrationUnavailable("The matched panel requires at least two mice per group for edgeR.")
    structural, audit = downsample_group_tables(panel, depth=depth, seed=seed)
    qualified = metadata.loc[metadata["sample_id"].isin(
        column for table in structural.values() for column in table.columns
    )].groupby("group_no")["mouse_id"].nunique()
    if qualified.reindex(MATCHED_GROUPS, fill_value=0).min() < 2:
        raise CalibrationUnavailable("Fewer than two mice in a matched group reached the structural depth.")
    empty = panel["g1"].iloc[:0, :0]
    controls = dict.fromkeys(EMPTY_BACKGROUND, empty)
    structural.update(controls)
    counts, units = prepare_counts({**panel, **controls}, metadata, pool_within_mouse=True)
    evidence, status = training_evidence(structural, [], counts, units, n_workers=n_workers)
    if evidence is None:
        raise CalibrationUnavailable(status)
    return evidence, structural, audit


def panel_cutoffs(evidence):
    profile = profile_coordinates(evidence)
    return fit_plot_thresholds(profile["initial_share_pct"], profile["retention_pct"])


def mouse_predictions(reference, candidates, external, metadata, baseline_prediction):
    results = []
    for sample in metadata.itertuples(index=False):
        held = subset_v_tables({sample.group_no: external[sample.group_no]}, [sample.sample_id])
        table = held[sample.group_no]
        if not len(table.columns):
            continue
        observed = attach_mouse_response(candidates, reference, table, sample.group_no, sample.mouse_id)
        results.append(observed.assign(baseline_prediction=baseline_prediction, external_run=sample.sample_id))
    if not results:
        return candidates.iloc[:0].assign(
            **{column: pd.Series(dtype=float) for column in [*MOUSE_COLUMNS, *RESPONSE_COLUMNS]},
            baseline_prediction=pd.Series(dtype=float), external_run=pd.Series(dtype=str),
        )
    return pd.concat(results, ignore_index=True)


def mouse_prediction_losses(predictions, metadata):
    groups = predictions.assign(
        absolute_error=(predictions["response"] - predictions["predicted_response"]).abs(),
        baseline_abs_error=(predictions["response"] - predictions["baseline_prediction"]).abs(),
    ).groupby(MOUSE_COLUMNS)
    losses = groups.agg(
        n_observed=("response", "count"), mae=("absolute_error", "mean"), baseline_mae=("baseline_abs_error", "mean"),
    ).reset_index()
    counts = metadata[["group_no", "mouse_id", "sample_id"]].rename(columns={
        "group_no": "heldout_group", "mouse_id": "heldout_mouse", "sample_id": "external_run",
    })
    report = counts.merge(losses, on=MOUSE_COLUMNS, how="left", validate="1:1")
    report["n_observed"] = report["n_observed"].fillna(0).astype(int)
    return report.assign(n_frozen_candidates=predictions["v_gene"].nunique())


def compare_shortlists(internal, external, internal_evidence, external_evidence):
    """Compare E over the union, charging one unit for each exclusive candidate."""
    score = (
        internal.set_index("v_gene")["E"].rename("E_internal").to_frame()
        .join(external.set_index("v_gene")["E"].rename("E_external"), how="outer")
    )
    score["passes_internal"] = score.index.isin(internal["v_gene"])
    score["passes_external"] = score.index.isin(external["v_gene"])
    common = score["passes_internal"] & score["passes_external"]
    score["membership"] = np.select(
        [common, score["passes_internal"]], ["both", "internal_only"], default="external_only",
    )
    score["delta_E"] = score["E_external"] - score["E_internal"]
    score["absolute_delta_E"] = score["delta_E"].abs()
    largest = score[["E_internal", "E_external"]].max(axis=1)
    score["relative_delta_E"] = (score["absolute_delta_E"] / largest.replace(0, np.nan)).where(
        largest.ne(0) | ~common, 0.0,
    )
    score["union_loss"] = score["relative_delta_E"].where(common, 1.0)
    for prefix, evidence in (("internal", internal_evidence), ("external", external_evidence)):
        coordinates = profile_coordinates(evidence).set_index("v_gene")[["initial_share_pct", "retention_pct"]]
        score = score.join(coordinates.add_prefix(f"{prefix}_"), validate="1:1")
    order = score[["E_internal", "E_external"]].max(axis=1)
    return score.loc[order.sort_values(ascending=False).index].reset_index()


def summarize_transfer(mouse_losses, compared):
    within_group = mouse_losses.groupby("heldout_group")[["mae", "baseline_mae"]].mean()
    valid_groups = {"g1", "g5"}.issubset(within_group["mae"].dropna().index)
    common = compared.loc[compared["membership"].eq("both")]
    n_union = len(compared)
    return {
        "n_mice_scored": int(mouse_losses["mae"].notna().sum()),
        "balanced_mouse_mae": float(within_group["mae"].mean()) if valid_groups else None,
        "balanced_baseline_mae": float(within_group["baseline_mae"].mean()) if valid_groups else None,
        "n_internal": int(compared["passes_internal"].sum()),
        "n_external": int(compared["passes_external"].sum()),
        "n_common": len(common),
        "jaccard": len(common) / n_union if n_union else None,
        "mean_abs_delta_E_common": float(common["absolute_delta_E"].mean()) if len(common) else None,
        "mean_union_loss": float(compared["union_loss"].mean()) if n_union else None,
    }
