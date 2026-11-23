#!/usr/bin/env python3
"""Draw the phenotype-specific Eq. 7 figures and export the manuscript error table."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_splenic import MODELS
from src.eq7_continuous import JOINT_CONTRAST

SUBSETS = {"cd4": "CD4+", "cd8": "CD8+", "pooled": "Pooled"}
COLORS = dict(zip(MODELS, ("#744b85", "#146682", "#118171"), strict=True))
NAMES = dict(zip(MODELS, ("aaV-indicator ridge", "V-segment-effect ridge", "aaV-fraction ridge"), strict=True))
RANKERS = {"structural_score": "Structural rank", "continuous_edger": "V-segment edgeR",
           "continuous_fisher": "V-segment Fisher"}
SUBSET_COLORS = {"cd4": "#744b85", "cd8": "#146682", "pooled": "#be622b"}
HIGHLIGHTS = ("TRBV2", "TRBV14", "TRBV17")
INK = "#252735"


def style():
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Liberation Serif", "DejaVu Serif"],
                         "mathtext.fontset": "stix", "font.size": 11, "axes.labelcolor": INK,
                         "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.spines.left": True, "axes.spines.bottom": True, "pdf.fonttype": 42,
                         "axes.edgecolor": "#9996a3", "xtick.color": INK, "ytick.color": INK})


def save(fig, output, name):
    fig.savefig(output / f"{name}.pdf", bbox_inches="tight", pad_inches=0.14)
    plt.close(fig)


def panel_label(ax, label):
    ax.text(-0.05, 1.04, f"({label})", transform=ax.transAxes, fontweight="bold", va="bottom")


def transfer_plot(runs, output):
    directory = runs / "TRB_cd4"
    internal = pd.read_csv(directory / "full_internal_structure.csv").set_index("v_gene")
    external = pd.read_csv(directory / "external_structure.csv").set_index("v_gene")
    candidates = pd.read_csv(directory / "frozen_candidates.csv")
    cutoffs = candidates.iloc[0][["training_share_cutoff_pct", "training_retention_cutoff_pct"]]
    genes = sorted(set(internal.index[internal["selected"]]) | set(external.index[external["selected"]]) | set(HIGHLIGHTS))
    rows = [(internal, "Internal marrow study", "#744b85", "o"),
            (external, "External skin-graft study", "#be622b", "s")]
    fig, axes = plt.subplots(1, 2, figsize=(8.7, max(3.1, 1.3 + len(genes) * .26)), sharey=True)
    for ax, field, cutoff, title, letter in zip(
        axes, ("initial_share_pct", "retention_pct"), cutoffs,
        ("Initial share of aaV (%)", "aaV retention after subtraction (%)"), "ab", strict=True,
    ):
        ax.axvline(cutoff, color="#797780", ls="--", lw=1.2, zorder=1)
        for frame, label, color, marker in rows:
            data = frame.reindex(genes)
            y = np.arange(len(genes)) + (0.12 if marker == "s" else -0.12)
            ax.scatter(data[field], y, c=color, marker=marker, s=31, label=label, zorder=3)
        ax.set(yticks=np.arange(len(genes)), yticklabels=genes, xlabel=title)
        ax.set_title(f"Internal cutoff {cutoff:.2f}%", loc="right", fontsize=10)
        ax.grid(axis="y", color="#e7e5ec", lw=.65)
        panel_label(ax, letter)
    axes[0].invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=2, loc="lower center", bbox_to_anchor=(.5, -.02), frameon=False)
    fig.subplots_adjust(wspace=.19, bottom=.23)
    save(fig, output, "figure_13")


def ranking_panel(ax, runs, chain, letter):
    methods = {**NAMES, **RANKERS}
    rows, labels = [], []
    for subset, label in SUBSETS.items():
        scores = pd.read_csv(runs / f"{chain}_{subset}" / "internal_ranking_summary.csv")
        for top_k, group in scores.groupby("top_k"):
            rows.append(group.set_index("model")[JOINT_CONTRAST].reindex(methods))
            labels.append(f"{label} · top {top_k}")
    coordinates = np.arange(len(rows))
    for i, (model, label) in enumerate(methods.items()):
        values = [row[model] for row in rows]
        marker = ("o", "s", "^", "x", "+", "d")[i]
        shade = COLORS.get(model, ("#474752", "#707582", "#b29575")[max(0, i - len(MODELS))])
        ax.scatter(values, coordinates + (i - 2.5) * .084, marker=marker, c=shade,
                   s=23, linewidths=.85, label=label, zorder=3)
    ax.set(yticks=coordinates, yticklabels=labels, xlabel="Observed joint recovery contrast", title=chain)
    ax.invert_yaxis()
    ax.axvline(0, color="#a9a8ae", lw=.8)
    ax.grid(axis="y", color="#e7e5ec", lw=.7)
    panel_label(ax, letter)


def error_table(runs):
    records = []
    for chain in ("TRA", "TRB"):
        for subset in SUBSETS:
            directory = runs / f"{chain}_{subset}"
            info = pd.read_json(directory / "run.json", typ="series")
            internal = pd.read_csv(directory / "internal_error_summary.csv")
            records.append(internal.assign(subset=subset, evaluation="Internal",
                                           splits=np.prod(list(info["mice"].values()))))
    directory = runs / "TRB_cd4"
    external = pd.read_csv(directory / "external_error_summary.csv")
    records.append(external.assign(subset="cd4", evaluation="External skin graft", splits=np.nan))
    results = pd.concat(records, ignore_index=True)
    results["relative_mae"] = results["mae"] / results["constant_mae"].replace(0, np.nan)
    results["is_best"] = False
    main = results["model"].isin(MODELS)
    for _, rows in results.loc[main].groupby(["evaluation", "chain", "subset", "penalty"], sort=False):
        results.loc[rows.index, "is_best"] = rows["mae"].le(rows[["mae", "constant_mae"]].min().min() + 1e-12)
    return results


def error_panel(ax, errors):
    primary = errors.loc[errors["penalty"].eq(1) & errors["model"].isin(MODELS)]
    record, labels = [], []
    for (evaluation, chain, subset), group in primary.groupby(["evaluation", "chain", "subset"], sort=False):
        record.append(group.set_index("model").reindex(MODELS))
        labels.append(f"{chain} · {SUBSETS[subset]}" if evaluation == "Internal" else "TRB · external CD4+")
    matrix = np.stack([np.r_[row["relative_mae"], 1.0] for row in record])
    if np.isfinite(matrix).any():
        upper = max(1.03, np.nanmax(matrix))
        lower = min(.97, np.nanmin(matrix))
    else:
        lower, upper = .97, 1.03
    ax.imshow(matrix, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=1, vmin=lower, vmax=upper), aspect="auto")
    for y, row in enumerate(record):
        winners = [*row["is_best"], not row["relative_mae"].lt(1 - 1e-12).any()]
        for x, (value, winner) in enumerate(zip(matrix[y], winners, strict=True)):
            dark_cell = np.isfinite(value) and (value - lower < .2 * (1 - lower) or
                                                upper - value < .2 * (upper - 1))
            ax.text(x, y, f"{value:.2f}" if np.isfinite(value) else "NA", ha="center", va="center",
                    color="white" if dark_cell else INK, fontweight="bold" if winner else "normal")
    ax.set(yticks=np.arange(len(labels)), yticklabels=labels,
           xticks=np.arange(len(MODELS) + 1),
           xticklabels=["aaV indicator", "V effect", "aaV fraction", "Constant"],
           title="MAE divided by training-mean MAE")
    ax.tick_params(axis="x", labelsize=9)
    ax.tick_params(length=0)
    ax.axhline(len(labels) - 1.5, color="white", lw=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_label(ax, "c")


def eligibility_panel(ax, structure):
    for i, (subset, label) in enumerate(SUBSETS.items()):
        genes = structure.loc[structure["subset"].eq(subset)].set_index("v_gene").reindex(HIGHLIGHTS)
        margin = 100 * (np.minimum(genes["initial_share_pct"] / genes["share_cutoff_pct"],
                                    genes["retention_pct"] / genes["retention_cutoff_pct"]) - 1)
        ax.scatter(margin, np.arange(len(HIGHLIGHTS)) + (i - 1) * .18, s=37,
                   marker=("o", "s", "D")[i], color=SUBSET_COLORS[subset], label=label, zorder=3)
        if subset == "cd4":
            segment = genes.loc["TRBV14"]
            detail = f"{segment['retention_pct']:.2f}% vs {segment['retention_cutoff_pct']:.2f}% retention"
            ax.annotate(detail, (margin["TRBV14"], .82), xytext=(.16, .58), textcoords="axes fraction",
                        fontsize=9, color=SUBSET_COLORS[subset],
                        arrowprops={"arrowstyle": "-", "color": "#aba5b3", "linewidth": .8})
    ax.axvline(0, color="#797780", ls="--", lw=1.2)
    ax.set(yticks=np.arange(len(HIGHLIGHTS)), yticklabels=HIGHLIGHTS, ylim=(2.46, -.46),
           xlabel="Relative margin to limiting cutoff (%)", title="TRB structural profile · both tissues")
    ax.grid(axis="y", color="#e7e5ec", lw=.7)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.22), ncol=3, frameon=False, fontsize=9)
    ax.text(-0.10, 1.04, "(d)", transform=ax.transAxes, fontweight="bold", va="bottom")


def comparison_plot(runs, output, errors, structure):
    fig, ((a, b), (c, d)) = plt.subplots(2, 2, figsize=(9.2, 7.1),
                                       gridspec_kw={"height_ratios": [1, .95]})
    ranking_panel(a, runs, "TRA", "a")
    ranking_panel(b, runs, "TRB", "b")
    handles, labels = a.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.52, .992),
               ncol=3, frameon=False, fontsize=9.5)
    error_panel(c, errors)
    eligibility_panel(d, structure)
    fig.subplots_adjust(left=.16, right=.99, bottom=.11, top=.87, wspace=.40, hspace=.60)
    save(fig, output, "figure_14")


def external_plot(runs, output):
    fig = plt.figure(figsize=(9, 5.5), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])
    a, b, c = fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])
    values = pd.read_csv(runs / "TRB_cd4" / "external_predictions.csv")
    values = values.loc[values["penalty"].eq(1)]
    observed = values.groupby("v_gene")["observed"].first()
    genes = observed.sort_values(ascending=False).index
    spots = np.arange(len(genes))
    for i, (model, name) in enumerate(NAMES.items()):
        scores = values.loc[values["model"].eq(model)].set_index("v_gene").reindex(genes)
        a.scatter(spots + (i - 1) * .10, scores["prediction"], color=COLORS[model],
                  marker=("o", "s", "^")[i], s=40, label=name, zorder=3)
    a.scatter(spots, observed.reindex(genes), marker="x", color=INK, s=48, label="Observed", zorder=4)
    a.axhline(values["constant_baseline"].iat[0], color="#797780", ls="--", label="Training mean")
    a.axhline(0, color="#d3d0d8", lw=.8)
    a.set(xticks=spots, xticklabels=genes, ylabel="Joint recovery contrast",
          title="Skin graft · splenic CD4+ TRB")
    a.legend(loc="upper right", ncol=2, frameon=False, fontsize=9)
    panel_label(a, "a")

    for chain, ax, letter in zip(("TRA", "TRB"), (b, c), "bc", strict=True):
        rows = pd.read_csv(runs / f"{chain}_pooled" / "external_mouse_recovery.csv")
        order = rows.groupby("v_gene")["joint_recovery"].mean().sort_values(ascending=False).index
        for i, (mouse, group) in enumerate(rows.groupby("mouse")):
            ax.scatter(np.arange(len(order)) + (i - 1) * .17,
                       group.set_index("v_gene").loc[order, "joint_recovery"], s=25,
                       color=("#744b85", "#146682", "#be622b")[i], marker=("o", "s", "^")[i], label=mouse)
        ax.set(xticks=np.arange(len(order)), xticklabels=order, ylabel="Joint recovery index",
               title=f"Marrow recipients · {chain}")
        ax.tick_params(axis="x", labelsize=8.5)
        plt.setp(ax.get_xticklabels(), rotation=60 if chain == "TRA" else 45,
                 horizontalalignment="right", rotation_mode="anchor")
        ax.grid(axis="y", color="#e7e5ec", lw=.7)
        panel_label(ax, letter)
    b.legend(loc="upper right", ncol=3, frameon=False, fontsize=9)
    save(fig, output, "figure_15")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--structural", required=True, type=Path, help="Archived all-tissue TRB profiles and cutoffs.")
    args = parser.parse_args()
    style()
    args.output.mkdir(parents=True, exist_ok=True)
    errors = error_table(args.runs)
    errors.to_csv(args.output / "model_errors.csv", index=False)
    transfer_plot(args.runs, args.output)
    comparison_plot(args.runs, args.output, errors, pd.read_csv(args.structural))
    external_plot(args.runs, args.output)


if __name__ == "__main__":
    main()
