#!/usr/bin/env python3
"""Render the V-segment benchmark and weighted recovery comparison."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.eq7_weighted import MOLECULE_CONTRAST
from src.publication_figures import save_publication_figure


PURPLE = "#6c257c"
ORANGE = "#e7794f"
PALE_PURPLE = "#9466a5"
LIGHT = "#e6dde5"
SELECTION = {
    "structural_screen": "Adaptive shortlist",
    "S_ungated": r"$S(v)$ without cutoffs",
    "initial_share_ungated": "Initial share",
    "edgeR_binary_ungated": "Binary aaV edgeR",
    "Fisher_binary_ungated": "Binary aaV Fisher",
    "old_E_heuristic_ungated": "Historical heuristic",
    "V_usage_edgeR_logFC": "V-usage edgeR effect",
    "V_usage_edgeR_significant": "V-usage edgeR with FDR",
}
RANKERS = {
    "legacy_score_refit": r"Original $E(v)$ refit",
    "structural_score": r"Structural score $S(v)$",
    "count_contrast_ridge": "Count contrast ridge",
    "weighted_contrast_ridge": "Weighted contrast ridge",
}


def style(ax, direction="x"):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis=direction, color=LIGHT, linewidth=0.55)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5)


def selection_figure(benchmark, usage):
    summary = pd.read_csv(benchmark / "internal_1249_selection_summary.csv")
    if usage:
        summary = pd.concat([summary, pd.read_csv(usage / "internal_v_usage_edger_summary.csv")])
    order = {method: label for method, label in SELECTION.items() if method in set(summary["method"])}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), sharey=True, layout="constrained")
    y = np.arange(len(order))
    for ax, chain in zip(axes, ("TRA", "TRB"), strict=True):
        values = 100 * summary.loc[summary["chain"].eq(chain)].set_index("method").loc[list(order), "g1_minus_mean_g5g6"]
        for idx, (method, value) in enumerate(values.items()):
            color = ORANGE if method == "structural_screen" else PURPLE
            ax.hlines(idx, 0, value, color=color, lw=1, alpha=0.42)
            ax.scatter(value, idx, color=color, s=30, zorder=3)
            ax.annotate(f"{value:.2f}", (value, idx), xytext=(4 if value >= 0 else -4, 0),
                        textcoords="offset points", ha="left" if value >= 0 else "right", va="center")
        ax.set(xlabel=f"{chain} contrast (percentage points)", yticks=y, yticklabels=list(order.values()))
        if values.min() >= 0:
            ax.set_xlim(-0.6, 11.3)
        else:
            ax.set_xlim(-1.85, 0.25)
        ax.axvline(0, color="black", lw=0.65)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        style(ax)
    axes[0].set_ylim(len(order) - 0.4, -0.6)
    return fig


def weighted_figure(results):
    ranks = pd.read_csv(results / "internal_ranking_summary.csv")
    internal = pd.read_csv(results / "internal_error_summary.csv")
    external = pd.read_csv(results / "external_error_summary.csv")
    predictions = pd.read_csv(results / "external_predictions.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), layout="constrained")
    for ax, chain in zip(axes[0], ("TRA", "TRB"), strict=True):
        frame = ranks.loc[ranks["chain"].eq(chain)].pivot(index="model", columns="top_k", values=MOLECULE_CONTRAST)
        positions = np.arange(len(RANKERS))
        for idx, (k, color, marker) in enumerate(zip(frame.columns, (PURPLE, ORANGE), ("o", "s"), strict=True)):
            ax.scatter(frame.loc[list(RANKERS), k], positions + 0.26 * (idx - 0.5),
                       c=color, marker=marker, s=31, label=f"Top {k}", zorder=3)
        ax.set(yticks=positions, yticklabels=list(RANKERS.values()), ylim=(len(RANKERS) - 0.4, -0.6),
               xlabel=f"{chain} contrast (weighted UMIs per million)")
        if chain == "TRB":
            ax.tick_params(labelleft=False)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.legend(frameon=False, loc="best")
        style(ax)
    ratios = pd.concat([
        internal.assign(label=internal["chain"].map({"TRA": "Internal TRA", "TRB": "Internal TRB"})),
        external.assign(label="External TRB"),
    ])
    ratios = ratios.loc[ratios["model"].eq("weighted_contrast_ridge")]
    ax = axes[1, 0]
    rows = ("Internal TRA", "Internal TRB", "External TRB")
    ax.axvline(1, color="black", linestyle=":", lw=0.8)
    for idx, (penalty, color, marker) in enumerate(zip((1, 10), (PURPLE, ORANGE), ("o", "s"), strict=True)):
        item = ratios.loc[ratios["penalty"].eq(penalty)].set_index("label").loc[list(rows)]
        ax.scatter(item["mae"] / item["constant_mae"], np.arange(len(rows)) + 0.24 * (idx - 0.5),
                   c=color, marker=marker, label=rf"$\lambda={penalty}$", s=30)
    ax.set(xlabel="MAE relative to Constant baseline", yticks=range(len(rows)), yticklabels=rows, ylim=(2.65, -0.55))
    ax.legend(frameon=False, loc="upper right")
    style(ax)
    ax = axes[1, 1]
    sample = predictions.loc[
        predictions["study"].eq("PRJNA954716") & predictions["model"].eq("weighted_contrast_ridge")
        & predictions["penalty"].eq(1)
    ].sort_values("v_gene", key=lambda v: v.str.extract(r"(\d+)")[0].astype(int))
    loc = np.arange(len(sample))
    ax.scatter(sample["observed"], loc - 0.19, s=26, c=PURPLE, zorder=3)
    ax.scatter(sample["prediction"], loc + 0.19, s=26, c=ORANGE, marker="s", zorder=3)
    ax.axvline(sample["constant_baseline"].iat[0], color="black", linestyle=":", lw=0.8)
    handles = [
        Line2D([], [], linestyle="", marker="o", color=PURPLE, label="Observed"),
        Line2D([], [], linestyle="", marker="s", color=ORANGE, label="Predicted"),
    ]
    ax.legend(handles=handles, frameon=False, loc="upper right")
    ax.set(xlabel="External contrast (weighted reads per million)", yticks=loc, yticklabels=sample["v_gene"],
           ylim=(len(sample) - 0.4, -0.85))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    style(ax)
    return fig


def ablation_figure(benchmark):
    summary = pd.read_csv(benchmark / "external_1249_frozen_summary.csv")
    order = {
        "calibration_mean": "Constant baseline",
        "published_frozen_full_E": r"Ridge intercept + $E(v)$",
        "drop_edgeR_no_refit": "Without edgeR",
        "drop_Fisher_no_refit": "Without Fisher",
        "drop_structural_no_refit": "Without residual count",
        "full_intercept_only": "Ridge intercept",
    }
    data = summary.loc[summary["study"].eq("PRJNA650133") & summary["chain"].eq("TRA")]
    values = data.set_index("method").loc[list(order), "group_balanced_mae"] * 100
    fig, ax = plt.subplots(figsize=(6.3, 2.65), layout="constrained")
    y = np.arange(len(order))
    ax.scatter(values, y, c=[ORANGE, PURPLE, PALE_PURPLE, PALE_PURPLE, PALE_PURPLE, ORANGE], s=30, zorder=3)
    ax.axvline(values.iat[0], color="black", linestyle=":", lw=0.75)
    ax.set(xlabel="TRA MAE in marrow-graft recipients (percentage points)", xlim=(8.55, 10.65),
           yticks=y, yticklabels=list(order.values()), ylim=(len(order) - 0.5, -0.5))
    for position, value in zip(y, values, strict=True):
        ax.annotate(f"{value:.3f}", (value, position), xytext=(6, 0), textcoords="offset points", va="center")
    style(ax)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--weighted", type=Path, required=True)
    parser.add_argument("--v-usage", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for name, figure in (
        ("eq7_benchmark_selection", selection_figure(args.benchmark, args.v_usage)),
        ("eq7_weighted_recovery", weighted_figure(args.weighted)),
        ("eq7_benchmark_frozen_ablation", ablation_figure(args.benchmark)),
    ):
        save_publication_figure(figure, args.output / f"{name}.png", panel_letters=False)
        plt.close(figure)


if __name__ == "__main__":
    main()
