#!/usr/bin/env python3
"""Plot the paired binary and continuous evidence comparison."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.plot_eq7_comparison import ORANGE, PURPLE, selection_figure, style
from src.eq7_continuous import BINARY_MODEL, CONTINUOUS_MODEL, JOINT_CONTRAST
from src.publication_figures import save_publication_figure


MODELS = {
    BINARY_MODEL: ("Binary-evidence ridge", PURPLE, "o"),
    CONTINUOUS_MODEL: ("Continuous-evidence ridge", ORANGE, "s"),
    "structural_score": ("Structural score", "#383838", "D"),
}


def spread(ax, labels):
    positions = np.arange(len(labels))
    ax.set(yticks=positions, yticklabels=labels, ylim=(len(labels) - 0.4, -0.6))
    ax.xaxis.set_major_locator(MaxNLocator(5))
    style(ax)
    return positions


def draw_markers(ax, position, values, methods=MODELS):
    for i, (model, (label, color, marker)) in enumerate(methods.items()):
        ax.scatter(values[model], position + 0.23 * (i - (len(methods) - 1) / 2),
                   s=28, color=color, marker=marker, label=label, zorder=3)


def legend(ax):
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), borderaxespad=0,
              frameon=False, fontsize=7, handletextpad=0.3, labelspacing=0.33)


def plots(output, results):
    ranks = pd.read_csv(results / "internal_ranking_summary.csv")
    errors = pd.read_csv(results / "internal_error_summary.csv").assign(panel=lambda f: "Internal " + f["chain"])
    external_file = results / "external_error_summary.csv"
    if external_file.is_file():
        errors = pd.concat([errors, pd.read_csv(external_file).assign(panel="External TRB")])
    for panel, chain in zip("ab", ("TRA", "TRB"), strict=True):
        fig, ax = plt.subplots(figsize=(3.45, 2.2), layout="constrained")
        frame = ranks.loc[ranks["chain"].eq(chain)].pivot(index="top_k", columns="model", values=JOINT_CONTRAST)
        y = spread(ax, [f"Top {k}" for k in frame.index])
        draw_markers(ax, y, frame)
        ax.set_xlabel(f"{chain} joint recovery contrast")
        legend(ax)
        save_publication_figure(fig, output / f"eq7_joint_{panel}.png", panel_letters=False, width_in=3.45)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.45, 2.2), layout="constrained")
    methods = {name: MODELS[name] for name in (BINARY_MODEL, CONTINUOUS_MODEL)}
    frame = errors.loc[errors["model"].isin(methods) & errors["penalty"].eq(1)].copy()
    frame["relative_mae"] = frame["mae"] / frame["constant_mae"]
    frame = frame.pivot(index="panel", columns="model", values="relative_mae").reindex(
        [panel for panel in ("Internal TRA", "Internal TRB", "External TRB") if panel in set(frame["panel"])]
    )
    y = spread(ax, frame.index)
    draw_markers(ax, y, frame, methods)
    ax.axvline(1, color="#383838", lw=0.75, ls=":")
    ax.set_xlabel("MAE relative to Constant baseline")
    legend(ax)
    save_publication_figure(fig, output / "eq7_joint_c.png", panel_letters=False, width_in=3.45)
    plt.close(fig)

    predictions_file = results / "external_predictions.csv"
    if predictions_file.is_file():
        data = pd.read_csv(predictions_file)
        data = data.loc[data["study"].eq("PRJNA954716") & data["penalty"].eq(1)]
    else:
        data = pd.DataFrame()
    fig, ax = plt.subplots(figsize=(3.45, 2.2), layout="constrained")
    if not data.empty:
        frame = data.pivot(index="v_gene", columns="model", values="prediction")
        frame["observed"] = data.drop_duplicates("v_gene").set_index("v_gene")["observed"]
        y = spread(ax, frame.index)
        draw_markers(ax, y, frame, {**methods, "observed": ("Observed", "#383838", "x")})
        ax.set_xlabel("External TRB joint recovery contrast")
    else:
        data = pd.read_csv(results / "internal_predictions.csv")
        data = data.loc[data["chain"].eq("TRB") & data["penalty"].eq(1) & data["model"].isin(methods)]
        for model, (label, color, marker) in methods.items():
            points = data.loc[data["model"].eq(model)]
            ax.scatter(points["observed"], points["prediction"], s=8, c=color, marker=marker,
                       alpha=0.45, label=label)
        ax.axline((0, 0), slope=1, color="#383838", ls=":", lw=0.75)
        ax.set(xlabel="Observed TRB contrast", ylabel="Predicted TRB contrast")
        style(ax)
    legend(ax)
    save_publication_figure(fig, output / "eq7_joint_d.png", panel_letters=False, width_in=3.45)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--benchmark", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plots(args.output, args.results)
    if args.benchmark:
        usage = args.benchmark if (args.benchmark / "internal_v_usage_edger_summary.csv").is_file() else None
        figure = selection_figure(args.benchmark, usage)
        save_publication_figure(figure, args.output / "eq7_benchmark_selection.png", panel_letters=False)
        plt.close(figure)


if __name__ == "__main__":
    main()
