"""Plot the design matrix, fitted responses, and V-segment scores for Equation 7."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator

from .article_figures import _grid, _ticks, empty_axis
from .eq7_calibration import CONTRIBUTIONS, STANDARDIZED_FEATURES


INPUT_LABELS = ("edgeR verdict", "Fisher verdict", "log(1 + residual aaV)")
TERM_COLORS = tuple(plt.get_cmap("magma")(v) for v in (0.28, 0.53, 0.73))
GROUP_MARKERS = dict(zip(("g1", "g5", "g6"), zip(("o", "s", "^"), TERM_COLORS), strict=True))


def plot_eq7_inputs(rows, weights, chain, status):
    fig, ax = plt.subplots(figsize=(7.4, 5.8), layout="constrained")
    ax.set_title(f"{chain} Equation 7 regression inputs", fontsize=13, pad=10)
    _ticks(ax, 9)
    if weights.empty:
        empty_axis(ax, status)
        return fig

    fitted = rows.loc[rows["used_for_fit"]].sort_values(
        ["heldout_group", "heldout_mouse", "evidence_split", "v_gene"],
    )
    values = fitted[list(STANDARDIZED_FEATURES)].to_numpy(dtype=float)
    extent = max(1, np.abs(values).max())
    heatmap = ax.imshow(
        values, aspect="auto", interpolation="nearest", cmap="magma",
        norm=Normalize(-extent, extent),
    )
    labels = [
        label if varied else f"{label}\n(constant)"
        for label, varied in zip(INPUT_LABELS, weights["varied_in_training"], strict=True)
    ]
    ax.set_xticks(range(len(labels)), labels)
    boundaries = fitted.groupby(["heldout_group", "heldout_mouse"], sort=False).size()
    stops = boundaries.cumsum()
    ax.set_yticks(
        stops.to_numpy() - boundaries.to_numpy() / 2 - 0.5,
        [mouse for group, mouse in boundaries.index],
    )
    for stop in stops.iloc[:-1]:
        ax.axhline(stop - 0.5, color="white", linewidth=0.8)
    ax.set_ylabel("Omitted mice")
    fig.colorbar(heatmap, ax=ax, label="Standardized predictor (z)", pad=0.04)
    return fig


def plot_eq7_regression(rows, outer, outer_status, chain, status):
    """Show the fit response distribution and mouse-level internal predictions."""
    fig, (fit, heldout) = plt.subplots(1, 2, figsize=(7.1, 3.9), layout="constrained")
    for ax in (fit, heldout):
        _ticks(ax, 8.5)
        _grid(ax, "y")
    fit.set_title("A. Calibration")
    heldout.set_title("B. Internal holdouts")
    fit.set(ylabel="Signed recovery (%)", xlabel="Calibration group")
    heldout.set(xlabel="Mean predicted recovery (%)", ylabel="Mean observed recovery (%)")

    if {"used_for_fit", "response", "fitted_response"} <= set(rows):
        training = rows.loc[rows["used_for_fit"]].dropna(subset=["response", "fitted_response"])
    else:
        training = rows.iloc[:0]
    if training.empty:
        empty_axis(fit, status)
    else:
        grouped = training.groupby("heldout_group")
        groups = [group for group in GROUP_MARKERS if group in grouped.groups]
        positions = np.arange(len(groups))
        for offset, column, tone, label in (
            (-0.18, "response", TERM_COLORS[0], "Observed"),
            (0.18, "fitted_response", TERM_COLORS[2], "Fitted"),
        ):
            series = [100 * grouped.get_group(group)[column].to_numpy() for group in groups]
            fit.boxplot(series, positions=positions + offset, widths=0.29, patch_artist=True,
                        showfliers=False, boxprops={"facecolor": tone, "alpha": 0.65},
                        medianprops={"color": "black"}, manage_ticks=False)
            fit.plot([], [], color=tone, linewidth=5, label=label)
        fit.axhline(0, color="black", linewidth=0.55)
        fit.set_xticks(positions, groups)
        fit.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)

    if {"heldout_group", "response", "predicted_response"} <= set(outer):
        controls = outer.loc[outer["heldout_group"].isin(("g5", "g6"))]
        means = controls.groupby(["heldout_group", "heldout_mouse"])[
            ["response", "predicted_response"]
        ].mean().dropna() * 100
    else:
        means = outer.iloc[:0]
    if means.empty:
        empty_axis(heldout, "No eligible internal predictions")
        return fig

    extent = max(1.0, np.abs(means.to_numpy()).max()) * 1.12
    heldout.plot([-extent, extent], [-extent, extent], "--", color="black", lw=0.8, zorder=1)
    for group, (marker, color) in GROUP_MARKERS.items():
        if group in means.index.get_level_values(0):
            values = means.loc[group]
            heldout.scatter(values["predicted_response"], values["response"], facecolors="none",
                            edgecolors=[color], marker=marker, label=group, s=35, zorder=3)
    heldout.set(xlim=(-extent, extent), ylim=(-extent, extent))
    heldout.locator_params(axis="both", nbins=5)
    heldout.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    return fig


def plot_ridge_selection(curve, chain):
    fig, ax = plt.subplots(figsize=(6.5, 4.1), layout="constrained")
    ax.set_title(f"{chain} ridge penalty selection", fontsize=13, pad=8)
    if curve["mse"].isna().all():
        message = "Control mice did not provide comparable tuning estimates"
        empty_axis(ax, message)
        return fig
    _ticks(ax, 9)
    _grid(ax, "y")
    ax.set_xscale("log")
    ax.errorbar(curve["ridge"], curve["mse"], yerr=curve["mse_se"].fillna(0),
                color=TERM_COLORS[0], marker="o", capsize=3, label="Ridge regression")
    ax.plot(curve["ridge"], curve["baseline_mse"], "--", color=TERM_COLORS[2], label="Intercept-only baseline")
    chosen = curve.loc[curve["selected"], "ridge"].item()
    ax.axvline(chosen, lw=1, color=TERM_COLORS[1], label=f"Selected λ = {chosen:g}")
    ax.set_xticks(curve["ridge"], [f"{ridge:g}" for ridge in curve["ridge"]])
    ax.set_xlabel("Ridge penalty (λ)")
    ax.set_ylabel("Mouse-balanced mean squared error")
    ax.legend(frameon=False, fontsize=8, loc="upper center",
              bbox_to_anchor=(0.5, -0.17), ncol=2)
    return fig


def plot_eq7_scores(scored, chain, status):
    height = max(3.3, 1.75 + 0.245 * len(scored))
    fig, ax = plt.subplots(figsize=(8.3, height), layout="constrained")
    ax.set_title(f"{chain} V segments passing both dynamic cutoffs", fontsize=13, pad=8)
    _ticks(ax, 9)
    _grid(ax, "x")
    ax.set_xlabel("Equation 7 ranking score (E)")
    if scored.empty:
        empty_axis(ax, "No V segments passed both cutoffs.")
    elif scored["E"].isna().all():
        empty_axis(ax, status)
    else:
        y = np.arange(len(scored))
        left = np.zeros(len(scored))
        for term, label, color in zip(CONTRIBUTIONS, INPUT_LABELS, TERM_COLORS, strict=True):
            part = scored[term].to_numpy()
            ax.barh(y, part, left=left, height=0.74, label=label, color=color)
            left += part
        padding = left.max() * 0.24 if left.max() > 0 else 1
        ax.set_xlim(0, left.max() + padding)
        for pos, total in zip(y, left, strict=True):
            ax.annotate(f"{total:.3g}", (total, pos), xytext=(4, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=8)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.set_yticks(y, scored["v_gene"])
        ax.invert_yaxis()
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3,
                  frameon=False, fontsize=8)
    return fig


def plot_external_predictions(predictions, losses, summary):
    """Show external coverage and paired prediction errors when estimable."""
    fig, (coverage, errors) = plt.subplots(1, 2, figsize=(7.1, 3.6), layout="constrained")
    for ax in (coverage, errors):
        _ticks(ax, 8.5)
        _grid(ax, "x")
    coverage.set_title("A. Coverage")
    errors.set_title("B. Prediction error")
    positions = np.arange(len(losses))
    labels = losses["heldout_mouse"].astype(str)
    nominal = losses["n_frozen_candidates"].iat[0] if len(losses) else 0
    coverage.set(xlabel="Frozen candidates observed", yticks=positions, yticklabels=labels)
    coverage.set_xlim(-0.35, max(1, nominal, losses["n_observed"].max() if len(losses) else 0) + 0.35)
    coverage.xaxis.set_major_locator(MaxNLocator(integer=True, min_n_ticks=2))
    coverage.invert_yaxis()
    for group, (marker, color) in GROUP_MARKERS.items():
        mask = losses["heldout_group"].eq(group).to_numpy()
        if mask.any():
            coverage.scatter(losses.loc[mask, "n_observed"], positions[mask], marker=marker,
                             color=color, label=group, s=32, zorder=2)
    coverage.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=3)
    scored = losses.dropna(subset=["mae", "baseline_mae"])
    if scored.empty:
        empty_axis(errors, "MAE unavailable\nNo mouse had an observed\nfrozen candidate")
        return fig
    y = np.arange(len(scored))
    ridge, baseline = 100 * scored["mae"], 100 * scored["baseline_mae"]
    for row, a, b in zip(y, ridge, baseline, strict=True):
        errors.plot([a, b], [row - 0.13, row + 0.13], color="0.55", lw=0.8)
    errors.scatter(ridge, y - 0.13, color=TERM_COLORS[0], label="Frozen score", zorder=2)
    errors.scatter(baseline, y + 0.13, color=TERM_COLORS[2], marker="s", label="Intercept only", zorder=2)
    errors.set(xlabel="Mouse MAE (percentage points)", yticks=y, yticklabels=scored["heldout_mouse"])
    errors.invert_yaxis()
    errors.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2)
    return fig


def plot_external_comparison(compared, summary):
    """Compare the structural coordinates against the frozen cutoffs."""
    fig, axes = plt.subplots(1, 2, figsize=(7.1, max(3.1, 1.25 + 0.32 * len(compared))),
                             sharey=True, layout="constrained")
    if compared.empty:
        for ax in axes:
            empty_axis(ax, "No candidates passed the frozen cutoffs")
        return fig

    y = np.arange(len(compared))
    for ax, coordinate, title, prefix in zip(
        axes, ("initial_share_pct", "retention_pct"),
        ("Initial share (%)", "Retention after subtraction (%)"), "AB", strict=True,
    ):
        _ticks(ax, 8.5)
        _grid(ax, "x")
        ax.set_title(f"{prefix}. {title}")
        ax.set_xlabel(title)
        ax.axvline(summary["frozen_gate"][coordinate], ls="--", color="black", linewidth=0.8)
        for cohort, color, marker, label in (
            ("internal", TERM_COLORS[0], "o", "Internal reference"),
            ("external", TERM_COLORS[2], "s", "Independent cohort"),
        ):
            offset = -0.12 if cohort == "internal" else 0.12
            ax.scatter(compared[f"{cohort}_{coordinate}"], y + offset, color=color,
                       marker=marker, label=label, s=25, zorder=2)
        ax.set_yticks(y, compared["v_gene"])
        ax.xaxis.set_major_locator(MaxNLocator(5))
    axes[0].invert_yaxis()
    axes[0].legend(loc="upper left", bbox_to_anchor=(0, -0.23), frameon=False, ncol=2)
    return fig
