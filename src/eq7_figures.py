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
        [f"{group} · {mouse}" for group, mouse in boundaries.index],
    )
    for stop in stops.iloc[:-1]:
        ax.axhline(stop - 0.5, color="white", linewidth=0.8)
    ax.set_ylabel("Omitted mice · one row per split and V segment")
    fig.colorbar(heatmap, ax=ax, label="Weighted standardized predictor (z)", pad=0.04)
    return fig


def plot_eq7_regression(rows, outer, outer_status, chain, status):
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.7))
    fig.subplots_adjust(left=0.09, right=0.98, wspace=0.29, top=0.81, bottom=0.28)
    fig.suptitle(f"{chain} Equation 7 fit and internal tuning", fontsize=13)
    panels = (
        (rows, "E", "Equation 7 score (E)", "A. Calibration observations", status),
        (outer, "predicted_response", "Predicted fraction (intercept + E)", "B. Omitted control mice",
         "No eligible predictions for control mice."),
    )
    values = [
        d[["response", predicted]].dropna().to_numpy()
        for d, predicted in ((rows, "fitted_response"), (outer, "predicted_response")) if predicted in d
    ]
    limit = max([0.04, *(np.abs(v).max() for v in values if v.size)]) * 1.08
    handles = {}
    for ax, (frame, x, xlabel, title, empty_message) in zip(axes, panels, strict=True):
        _ticks(ax, 9)
        _grid(ax, "both")
        ax.set_title(title, fontsize=10, pad=8)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel("Observed signed recovery fraction", fontsize=9)
        if x not in frame or frame.empty:
            empty_axis(ax, empty_message)
            continue
        frame = frame.dropna(subset=[x, "response"])
        if x != "E":
            frame = frame.loc[frame["heldout_group"].isin(("g5", "g6"))]
        if frame.empty:
            empty_axis(ax, empty_message)
            continue
        ax.set_ylim(-limit, limit)
        if x == "E":
            intercept = (frame["fitted_response"] - frame["E"]).iat[0]
            lower, upper = frame["E"].min(), frame["E"].max()
            padding = max(0.001, 0.14 * (upper - lower))
            ax.set_xlim(lower - padding, upper + padding)
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.ticklabel_format(axis="x", useOffset=False)
            if lower == upper:
                ax.axhline(intercept + lower, lw=1.1, color="black", zorder=1,
                           label="One fitted value")
            else:
                ax.plot([lower, upper], [intercept + lower, intercept + upper], lw=1.1, color="black",
                        zorder=1, label="Fitted line: intercept + E")
            ax.legend(loc="lower right", frameon=False, fontsize=7.6)
        else:
            ax.set_xlim(-limit, limit)
            ax.plot([-limit, limit], [-limit, limit], "--", lw=0.8, color="black", zorder=1,
                    label="Equal observed and predicted")
            ax.legend(loc="lower right", frameon=False, fontsize=7.6)
        for group, (marker, color) in GROUP_MARKERS.items():
            subset = frame.loc[frame["heldout_group"].eq(group)]
            if not subset.empty:
                points = subset.groupby([x, "response"], as_index=False).size()
                sizes = 18 + 27 * np.sqrt(points["size"] / points["size"].max())
                handles[group] = ax.scatter(
                    points[x], points["response"], marker=marker, facecolors="none", edgecolors=[color],
                    s=sizes, alpha=0.65 if x == "E" else 0.85, linewidths=0.7, zorder=2,
                )
        n_mice = len(frame[["heldout_group", "heldout_mouse"]].drop_duplicates())
        ax.text(0.04, 0.96, f"{n_mice} mice · {len(frame)} observations\nOverlapping points are grouped",
                transform=ax.transAxes,
                ha="left", va="top", fontsize=7.8)
    if handles:
        groups = [group for group in GROUP_MARKERS if group in handles]
        fig.legend([handles[group] for group in groups], groups, loc="lower center",
                   bbox_to_anchor=(0.5, 0.025), ncol=len(groups), frameon=False, fontsize=9)
    if not outer_status.empty:
        missing = {"g5", "g6"}.difference(outer_status.loc[outer_status["status"].eq("scored"), "outer_group"])
        if missing:
            fig.text(0.73, 0.028, f"No eligible tuning scores for {', '.join(sorted(missing))}.",
                     ha="center", fontsize=7.8)
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
    ax.set_ylabel("Mean squared error, balanced across g5 and g6")
    ax.legend(frameon=False, fontsize=8)
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
        padding = max(1e-3, left.max() * 0.15)
        ax.set_xlim(0, left.max() + padding)
        for pos, total in zip(y, left, strict=True):
            ax.annotate(f"{total:.3g}", (total, pos), xytext=(4, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=8)
        ax.set_yticks(y, scored["v_gene"])
        ax.invert_yaxis()
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3,
                  frameon=False, fontsize=8)
    return fig


def plot_external_predictions(predictions, losses, summary):
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.5), layout="constrained", sharex=True, sharey=True)
    fig.suptitle("TRB recovery in independent mice · PRJNA954716", fontsize=13)
    present = predictions.dropna(subset=["response", "predicted_response"])
    extent = max(0.05, np.abs(present[["response", "predicted_response"]].to_numpy()).max(initial=0)) * 1.15
    panels = {"g1": "Intact BALB/c", "g5": "BALB/c with C57BL/6 skin graft"}
    for ax, (group, title) in zip(axes, panels.items(), strict=True):
        ax.set_title(title, fontsize=10, pad=8)
        subset = present.loc[present["heldout_group"].eq(group)]
        if subset.empty:
            empty_axis(ax, "No selected V segment observed in this group")
            continue
        _ticks(ax, 8.5)
        _grid(ax, "both")
        ax.plot([-extent, extent], [-extent, extent], "--", lw=0.8, color="black", zorder=1)
        mice = subset.groupby("heldout_mouse")
        colors = plt.get_cmap("magma")(np.linspace(0.22, 0.78, mice.ngroups))
        for (mouse, frame), color in zip(mice, colors, strict=True):
            ax.scatter(frame["predicted_response"], frame["response"], label=mouse, c=[color],
                       alpha=0.82, s=36, linewidths=0.4, edgecolors="black", zorder=2)
        group_losses = losses.loc[losses["heldout_group"].eq(group)]
        ax.text(0.03, 0.97, f"{mice.ngroups}/{len(group_losses)} mice · MAE {group_losses['mae'].mean():.3f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=8)
        ax.set(xlim=(-extent, extent), ylim=(-extent, extent), aspect="equal")
        legend_loc = "lower right" if group == "g1" else "upper right"
        ax.legend(loc=legend_loc, frameon=False, fontsize=8, markerscale=0.75)
    fig.supxlabel("Predicted fraction from frozen internal reference (intercept + E)", fontsize=9)
    fig.supylabel("Observed signed recovery fraction", fontsize=9)
    if summary["balanced_mouse_mae"] is not None:
        fig.suptitle(
            f"TRB recovery in independent mice · PRJNA954716\n"
            f"Balanced MAE {summary['balanced_mouse_mae']:.3f} · "
            f"intercept-only {summary['balanced_baseline_mae']:.3f}", fontsize=12,
        )
    return fig


def plot_external_comparison(compared, summary):
    fig, ax = plt.subplots(figsize=(8.0, max(3.2, 1.65 + 0.3 * len(compared))), layout="constrained")
    ax.set_title("TRB candidates in the internal and independent cohorts", fontsize=12, pad=8)
    if compared.empty:
        empty_axis(ax, "No V segments passed the frozen cutoffs in either cohort")
        return fig
    _ticks(ax, 8.5)
    _grid(ax, "x")
    y = np.arange(len(compared))
    common = compared["membership"].eq("both")
    ax.hlines(y[common], compared.loc[common, "E_internal"], compared.loc[common, "E_external"],
              colors=TERM_COLORS[1], linewidth=1.2, alpha=0.65, zorder=1)
    for column, color, marker, label in (
        ("E_internal", TERM_COLORS[0], "o", "Internal · bone marrow"),
        ("E_external", TERM_COLORS[2], "s", "Independent · skin"),
    ):
        exists = compared[column].notna()
        ax.scatter(compared.loc[exists, column], y[exists], c=[color], marker=marker, label=label,
                   s=35, linewidths=0.45, edgecolors="black", zorder=2, clip_on=False)
    cutoff_note = f"{summary['n_common']}/{len(compared)} shared · union loss {summary['mean_union_loss']:.3f}"
    ax.set(xlabel=f"Equation 7 score with frozen coefficients\n{cutoff_note}", yticks=y, yticklabels=compared["v_gene"])
    ax.invert_yaxis()
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.31), ncol=2, frameon=False, fontsize=8)
    return fig
