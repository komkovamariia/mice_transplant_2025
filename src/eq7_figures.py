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

    fitted = rows.loc[rows["used_for_fit"]].sort_values(["heldout_group", "heldout_mouse", "v_gene"])
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
    ax.set_yticks(stops.to_numpy() - boundaries.to_numpy() / 2 - 0.5,
                  boundaries.index.get_level_values("heldout_mouse"))
    for stop in stops.iloc[:-1]:
        ax.axhline(stop - 0.5, color="white", linewidth=0.8)
    ax.set_ylabel("One row per omitted mouse × V segment")
    fig.colorbar(heatmap, ax=ax, label="Weighted standardized predictor (z)", pad=0.04)
    return fig


def plot_eq7_regression(rows, outer, outer_status, chain, status):
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.7))
    fig.subplots_adjust(left=0.09, right=0.98, wspace=0.28, top=0.83, bottom=0.22)
    fig.suptitle(f"{chain} Equation 7 regression and held-out checks", fontsize=13)
    panels = (
        (rows, "E", "Equation 7 score (E)", "A. Points used to learn weights", status),
        (outer, "predicted_response", "Predicted response (intercept + E)", "B. Outer mouse holdouts",
         "No eligible outer mouse predictions."),
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
        ax.set_ylabel("Observed signed aaV overlap", fontsize=9)
        if x not in frame or frame.empty:
            empty_axis(ax, empty_message)
            continue
        frame = frame.dropna(subset=[x, "response"])
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
                ax.axhline(intercept + lower, lw=1.1, color="#454545", zorder=1,
                           label="One fitted value")
            else:
                ax.plot([lower, upper], [intercept + lower, intercept + upper], lw=1.1, color="#454545",
                        zorder=1, label="Fitted line: intercept + E")
            ax.legend(loc="lower right", frameon=False, fontsize=7.6)
        else:
            ax.set_xlim(-limit, limit)
            ax.plot([-limit, limit], [-limit, limit], "--", lw=0.8, color="#999999", zorder=1,
                    label="Equal observed and predicted")
            ax.legend(loc="lower right", frameon=False, fontsize=7.6)
        for group, (marker, color) in GROUP_MARKERS.items():
            subset = frame.loc[frame["heldout_group"].eq(group)]
            if not subset.empty:
                handles[group] = ax.scatter(
                    subset[x], subset["response"], marker=marker, c=[color], alpha=0.65,
                    s=27, linewidths=0.35, edgecolors="#303030", zorder=2,
                )
        n_mice = len(frame[["heldout_group", "heldout_mouse"]].drop_duplicates())
        ax.text(0.04, 0.96, f"{n_mice} mice · {len(frame)} V–mouse observations", transform=ax.transAxes,
                ha="left", va="top", fontsize=7.8)
    if handles:
        groups = [group for group in GROUP_MARKERS if group in handles]
        fig.legend([handles[group] for group in groups], groups, loc="lower center",
                   bbox_to_anchor=(0.5, 0.055), ncol=len(groups), frameon=False, fontsize=9)
    if not outer_status.empty:
        missing = set(GROUP_MARKERS).difference(outer_status.loc[outer_status["status"].eq("scored"), "outer_group"])
        if missing:
            fig.text(0.73, 0.028, f"No eligible outer scores for {', '.join(sorted(missing))}.",
                     ha="center", fontsize=7.8)
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
