"""Render the splenic study figures using the manuscript's Palatino typography."""

import argparse
import sys
from pathlib import Path
from shutil import copyfile

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from pypdf import PdfReader, PdfWriter, Transformation

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.plot_eq7_splenic import NAMES
from scripts.run_splenic_figures import DISPLAY_SEED, METHODS
from src.article_figures import (COUNT_COLUMN, RETENTION_COLUMN, SHARE_COLUMN, eligible_segments,
                                 plot_dumbbell, plot_fisher, plot_localization, plot_subsets)
from src.repertoire_depth import plot_depth_check, plot_mds
from src.v_stability import plot_recovery, plot_regional_stability
from src.v_thresholds import PlotThresholds

PURPLE, ORANGE, GREEN = "#793e86", "#d16c48", "#347d70"
LIGHT = "#e8dfe8"
COLORS = dict(zip(NAMES, (PURPLE, ORANGE, GREEN), strict=True))
LABELS = {"cd4": "CD4+", "cd8": "CD8+", "pooled": "CD4+ and CD8+"}
CHAINS = ("TRA", "TRB")
METHOD_LABELS = {**NAMES, "structural_score": r"Structural rank $S(v)$",
                 "continuous_edger": "V-segment edgeR", "continuous_fisher": "V-segment Fisher"}


def style():
    mpl.rcParams.update({"font.family": ["P052", "serif"], "mathtext.fontset": "stix",
                         "font.size": 9, "xtick.labelsize": 8.4, "ytick.labelsize": 8.4,
                         "pdf.fonttype": 3, "axes.linewidth": .65, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlesize": 9.5, "axes.labelsize": 9,
                         "legend.fontsize": 8, "axes.unicode_minus": False})


def save(fig, path):
    fig.savefig(path, bbox_inches="tight", pad_inches=.08, facecolor="white")
    plt.close(fig)


def tidy(ax, axis="x"):
    ax.grid(axis=axis, color=LIGHT, linewidth=.6)
    ax.set_axisbelow(True)
    ax.tick_params(length=2.5, width=.6)


def tag(ax, letter):
    ax.text(0, 1.055, f"({letter})", transform=ax.transAxes, ha="left", va="bottom")


def mount_panels(panels, destination, columns=1):
    pages = [PdfReader(p).pages[0] for p in panels]
    total_width, gap, pad = 522, 12, 5
    column_width = (total_width - 2*pad - (columns-1)*gap) / columns
    heights = [float(p.mediabox.height) * column_width / float(p.mediabox.width) for p in pages]
    row_heights = [max(heights[i:i+columns]) for i in range(0, len(pages), columns)]
    scale_all = min(1, (480 - gap*(len(row_heights)-1)) / sum(row_heights))
    writer = PdfWriter()
    height = (sum(row_heights)*scale_all) + gap * (len(row_heights)-1) + 2*pad
    mounted = writer.add_blank_page(total_width, height)
    y = height - pad
    start_x = (total_width - columns*column_width*scale_all - (columns-1)*gap) / 2
    for index, (page, page_height) in enumerate(zip(pages, heights, strict=True)):
        row, column = divmod(index, columns)
        scale = scale_all * column_width / float(page.mediabox.width)
        x = start_x + column * (column_width*scale_all + gap)
        mounted.merge_transformed_page(page, Transformation().scale(scale).translate(x, y-page_height*scale_all))
        if column == columns-1:
            y -= row_heights[row]*scale_all + gap
    with destination.open("wb") as handle:
        writer.write(handle)


def cutoff(data, chain, subset):
    rows = pd.read_csv(data / chain / "thresholds.csv", float_precision="round_trip")
    row = rows.loc[rows["level"].eq("final") & rows["subset"].eq(subset)].iloc[0]
    return PlotThresholds(row["initial_share_pct"], row["retention_pct"])


def profiles(data):
    return {(chain, subset): pd.read_csv(data / chain / f"profile_{subset}.csv", float_precision="round_trip")
            for chain in CHAINS for subset in LABELS}


def retention_figure(data, sets):
    fig, axes = plt.subplots(2, 2, figsize=(7.65, 6.65))
    fig.subplots_adjust(left=.095, right=.86, top=.92, bottom=.17, wspace=.63, hspace=.53)
    norm = mpl.colors.Normalize(0, 1)
    cap = max(sets[(chain, s)][COUNT_COLUMN].max() for chain in CHAINS for s in ("cd4", "cd8"))
    for row, chain in enumerate(CHAINS):
        x_max = max(sets[(chain, s)][SHARE_COLUMN].max() for s in ("cd4", "cd8")) * 1.09
        for col, subset in enumerate(("cd4", "cd8")):
            ax = axes[row, col]
            limits = cutoff(data, chain, subset)
            frame = sets[(chain, subset)]
            chosen = frame[frame["selected"]].sort_values(["score_s", "v"], ascending=[False, True]).head(10)
            ax.scatter(frame[SHARE_COLUMN], frame[RETENTION_COLUMN], c=frame["region_entropy"],
                       cmap="magma", norm=norm, s=12 + 150 * frame[COUNT_COLUMN] / cap,
                       linewidths=0, alpha=.83, zorder=3)
            ax.axvline(limits.initial_share_pct, color=PURPLE, ls=":", lw=.8)
            ax.axhline(limits.retention_pct, color=PURPLE, ls=":", lw=.8)
            ax.set(xlim=(0, x_max), ylim=(0, 105), title=f"{chain}  {LABELS[subset]}")
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.set_xlabel("Initial aaV share (%)")
            if col == 0:
                ax.set_ylabel("Retention after subtraction (%)")
            if col:
                ax.tick_params(labelleft=False)
            ax.text(.01, .035, f"Cutoffs  {limits.initial_share_pct:.2f}%, {limits.retention_pct:.2f}%",
                    transform=ax.transAxes, fontsize=7.4, color=PURPLE)
            for location, rec in zip(np.linspace(.92, .20, len(chosen)),
                                     chosen.sort_values(RETENTION_COLUMN, ascending=False).to_dict("records"),
                                     strict=True):
                ax.annotate(rec["v"], (rec[SHARE_COLUMN], rec[RETENTION_COLUMN]),
                            (1.035, location), textcoords="axes fraction", ha="left", va="center",
                            fontsize=7.8, annotation_clip=False, clip_on=False,
                            arrowprops={"arrowstyle": "-", "color": ".45", "lw": .4, "shrinkB": 3})
            tag(ax, "abcd"[row * 2 + col])
    bar_ax = fig.add_axes((.29, .065, .40, .022))
    bar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="magma"), cax=bar_ax, orientation="horizontal")
    bar.set_label("Normalized Shannon entropy of overlap regions")
    bar.ax.tick_params(labelsize=8, width=.6, length=2.5)
    return fig


def threshold_sensitivity(data):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.3), layout="constrained")
    for r, chain in enumerate(CHAINS):
        raw = pd.read_csv(data / chain / "thresholds.csv", float_precision="round_trip")
        rows = raw[raw["level"].eq("mouse") & raw["subset"].eq("pooled")]
        limits = cutoff(data, chain, "pooled")
        for c, (field, name, fixed) in enumerate(zip(
            ("initial_share_pct", "retention_pct"), ("Initial share cutoff (%)", "Retention cutoff (%)"),
            limits, strict=True,
        )):
            ax = axes[r, c]
            ax.scatter(rows[field], np.arange(len(rows)), c=PURPLE, s=22)
            ax.axvline(fixed, color=".35", lw=.7, ls=":")
            ax.set(xlabel=name, title=chain, yticks=np.arange(len(rows)),
                   yticklabels=rows["omitted_mouse"].str.replace("_mm", " M", regex=False))
            ax.invert_yaxis()
            if c:
                ax.tick_params(labelleft=False)
            tag(ax, "abcd"[r * 2 + c])
            tidy(ax)
    return fig


def abundance_plot(data):
    fig, axes = plt.subplots(2, 1, figsize=(7.15, 6.65), layout="constrained")
    for ax, chain, letter in zip(axes, CHAINS, "ab", strict=True):
        table = pd.read_csv(data / chain / "abundance.csv")
        values = table.loc[table["n_edger_significant"].gt(0), ["v_gene", "mean_edger_log2fc"]]
        effect = "mean_edger_log2fc"
        columns = pd.concat([values.nlargest(8, effect), values.nsmallest(8, effect)]).drop_duplicates()
        columns = columns.sort_values(effect, ascending=True)
        ax.barh(columns["v_gene"], columns[effect], color=np.where(columns[effect].gt(0), ORANGE, PURPLE))
        ax.axvline(0, color="black", lw=.75)
        ax.set(xlabel=r"Mean signed $log_2$ fold change among significant aaV", title=chain)
        tidy(ax)
        tag(ax, letter)
    return fig


def selection_plot(data):
    fig, axes = plt.subplots(1, 2, figsize=(7.25, 3.40), sharey=True, layout="constrained")
    y = np.arange(len(METHODS))
    for ax, chain, letter in zip(axes, CHAINS, "ab", strict=True):
        frame = pd.read_csv(data / chain / "selection.csv").set_index("method").loc[list(METHODS)]
        v = frame["contrast_pp"].to_numpy()
        shade = [ORANGE, *([PURPLE] * (len(y) - 1))]
        ax.hlines(y, 0, v, color=shade, lw=1, alpha=.45)
        ax.scatter(v, y, color=shade, s=28)
        for xpos, ypos in zip(v, y):
            ax.annotate(f"{xpos:.2f}", (xpos, ypos), xytext=(4 if xpos >= 0 else -4, 0),
                        textcoords="offset points", va="center", ha="left" if xpos >= 0 else "right")
        extent = max(abs(v)) * .25
        ax.set(xlim=(min(0, min(v)-extent), max(0, max(v)+extent)), yticks=y, yticklabels=METHODS.values(),
               ylim=(len(y)-.45, -.6), title=chain, xlabel="Recovery contrast (percentage points)")
        ax.axvline(0, lw=.65, color=".35")
        ax.xaxis.set_major_locator(MaxNLocator(4))
        tag(ax, letter)
        tidy(ax)
    return fig


def external_structure(spleen):
    root = spleen / "TRB_cd4"
    it = pd.read_csv(root / "full_internal_structure.csv").set_index("v_gene")
    et = pd.read_csv(root / "external_structure.csv").set_index("v_gene")
    cand = pd.read_csv(root / "frozen_candidates.csv").iloc[0]
    union = it.index[it["selected"]].union(et.index[et["selected"]])
    labels = et.reindex(union)["initial_share_pct"].sort_values(ascending=False).index
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.0), sharey=True, layout="constrained")
    y = np.arange(len(labels))
    for ax, metric, threshold, label, letter in zip(
        axes, ("initial_share_pct", "retention_pct"), cand[["training_share_cutoff_pct", "training_retention_cutoff_pct"]],
        ("Initial aaV share (%)", "Retention after subtraction (%)"), "ab", strict=True,
    ):
        ax.axvline(threshold, ls=":", color="black", lw=.7)
        for source, offset, color, symbol in ((it, -.14, PURPLE, "o"), (et, .14, ORANGE, "s")):
            ax.scatter(source.reindex(labels)[metric], y+offset, color=color, marker=symbol, s=25, zorder=3)
        ax.set(xlabel=label, title=f"Internal cutoff {threshold:.2f}%", yticks=y, yticklabels=labels,
               ylim=(len(y)-.45, -.55))
        tidy(ax, "y")
        tag(ax, letter)
    handles = [Line2D([], [], marker=m, ls="", color=c, label=lab, ms=5)
               for m, c, lab in (("o", PURPLE, "Internal spleen CD4+ TRB"), ("s", ORANGE, "External spleen CD4+ TRB"))]
    axes[0].legend(handles=handles, loc="lower left", bbox_to_anchor=(0, -.33), ncol=2, frameon=False)
    return fig


def ranking_plot(spleen):
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 3.65))
    fig.subplots_adjust(left=.14, right=.98, bottom=.25, top=.87, wspace=.44)
    markers = ("o", "s", "^", "x", "+", "D")
    colors = (PURPLE, ORANGE, GREEN, "#575059", "#836c82", "#a48a75")
    for ax, chain, letter in zip(axes, CHAINS, "ab", strict=True):
        frames = []
        for subset in LABELS:
            t = pd.read_csv(spleen / f"{chain}_{subset}" / "internal_ranking_summary.csv")
            frames.append(t.assign(panel=t["top_k"].map(lambda n: f"{LABELS[subset]}  top {n}")))
        rows = pd.concat(frames).pivot(index="panel", columns="model", values="joint_recovery_contrast")
        labels = pd.concat(frames)["panel"].drop_duplicates()
        rows = rows.reindex(labels)
        y = np.arange(len(rows))
        for j, (model, name) in enumerate(METHOD_LABELS.items()):
            ax.scatter(rows[model], y+(j-2.5)*.067, color=colors[j], marker=markers[j],
                       label=name, s=23, linewidth=.75, zorder=3)
        maximum = np.nanmax(rows.to_numpy())
        minimum = min(0, np.nanmin(rows.to_numpy()))
        ax.set(xlim=(minimum, maximum*1.10), ylim=(len(y)-.43, -.5), yticks=y,
               yticklabels=[s.replace("CD4+ and CD8+", "Pooled") for s in rows.index], title=chain)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        tag(ax, letter)
        tidy(ax, "y")
    fig.supxlabel("Observed joint recovery contrast", y=.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, frameon=False, loc="lower center", bbox_to_anchor=(.52, .005))
    return fig


def error_plot(spleen):
    records = []
    for chain in CHAINS:
        for subset in LABELS:
            table = pd.read_csv(spleen / f"{chain}_{subset}" / "internal_error_summary.csv")
            records.append(table.assign(panel=f"{chain}  {'Pooled' if subset == 'pooled' else LABELS[subset]}"))
    external = pd.read_csv(spleen / "TRB_cd4" / "external_error_summary.csv")
    records.append(external.assign(panel="External  TRB CD4+"))
    table = pd.concat(records, ignore_index=True)
    table = table.loc[table["model"].isin(NAMES)].copy()
    table["ratio"] = table["mae"] / table["constant_mae"]
    labels = table["panel"].drop_duplicates()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.15), sharey=True)
    fig.subplots_adjust(left=.16, right=.98, bottom=.24, top=.87, wspace=.18)
    y = np.arange(len(labels))
    lower = min(.51, table["ratio"].min() - .06)
    for ax, penalty, letter in zip(axes, (1, 10), "ab", strict=True):
        frame = table.loc[table["penalty"].eq(penalty)].pivot(index="panel", columns="model", values="ratio").reindex(labels)
        best = frame.min(axis=1)
        ax.axvline(1, color="black", ls=":", lw=.75)
        ax.axhline(len(y)-1.5, color=LIGHT, lw=.8)
        for j, (model, name) in enumerate(NAMES.items()):
            offsets = y + (j - 1) * .18
            vals = frame[model].to_numpy()
            ax.scatter(vals, offsets, c=COLORS[model], marker=("o", "s", "^")[j], s=28, label=name, zorder=3)
            for value, ypos, winner in zip(vals, offsets, frame[model].eq(best) & best.lt(1), strict=True):
                if winner:
                    dx = 6 if value < .62 else -6
                    ax.annotate(f"{value:.2f}", (value, ypos), xytext=(dx, 0), textcoords="offset points",
                                ha="left" if dx > 0 else "right", va="center", fontsize=7.7)
        ax.set(yticks=y, yticklabels=labels, ylim=(len(y)-.42, -.55),
               xlim=(lower, 1.045), title=rf"$\lambda={penalty}$")
        ax.xaxis.set_major_locator(MaxNLocator(5))
        tidy(ax, "y")
        tag(ax, letter)
    fig.supxlabel("MAE divided by training-mean MAE", y=.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, frameon=False, loc="lower center", bbox_to_anchor=(.52, .025))
    return fig


def eligibility_plot(data, sets):
    frame = []
    for subset in LABELS:
        profile = sets[("TRB", subset)].copy()
        t = cutoff(data, "TRB", subset)
        profile["share_margin"] = profile[SHARE_COLUMN] - t.initial_share_pct
        profile["retention_margin"] = profile[RETENTION_COLUMN] - t.retention_pct
        frame.append(profile.assign(subset=subset).set_index("v"))
    genes = sorted(set.union(*(set(x.index[x["selected"]]) for x in frame)) | {"TRBV2", "TRBV14", "TRBV17"})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.75), layout="constrained", sharey=True)
    y = np.arange(len(genes))
    for ax, field, title, letter in zip(axes, ("share_margin", "retention_margin"),
                                        ("Initial share minus cutoff (pp)", "Retention minus cutoff (pp)"),
                                        "ab", strict=True):
        ax.axvline(0, color="black", ls=":", lw=.75)
        for j, (part, color, marker) in enumerate(zip(frame, (PURPLE, ORANGE, GREEN), ("o", "s", "^"), strict=True)):
            panel = part.reindex(genes)
            passes = panel["selected"].fillna(False).to_numpy()
            ax.scatter(panel[field], y+(j-1)*.22, marker=marker, s=25,
                       facecolors=[color if selected else "white" for selected in passes],
                       edgecolors=color, linewidth=.75, zorder=3)
        ax.set(xlabel=title, yticks=y, yticklabels=genes, ylim=(len(y)-.5, -.55))
        tidy(ax, "y")
        tag(ax, letter)
    handles = [Line2D([], [], marker=m, ls="", color=c, ms=5, label=label)
               for m, c, label in zip(("o", "s", "^"), (PURPLE, ORANGE, GREEN), LABELS.values(), strict=True)]
    handles.append(Line2D([], [], marker="o", ls="", markerfacecolor="white", color=".4", ms=5,
                          label="Fails one or both cutoffs"))
    axes[0].legend(handles=handles, frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0, -.36))
    return fig


def prediction_plot(spleen):
    fig = plt.figure(figsize=(7.35, 5.45), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])
    ax = fig.add_subplot(grid[0, :])
    values = pd.read_csv(spleen / "TRB_cd4" / "external_predictions.csv").query("penalty == 1")
    observed = values.groupby("v_gene")["observed"].first().sort_values(ascending=False)
    genes = observed.index
    x = np.arange(len(genes))
    for j, (model, label) in enumerate(NAMES.items()):
        rows = values.loc[values["model"].eq(model)].set_index("v_gene").reindex(genes)
        ax.scatter(x+(j-1)*.10, rows["prediction"], color=COLORS[model], s=28,
                   marker=("o", "s", "^")[j], label=label, zorder=3)
    ax.scatter(x, observed, color="black", s=37, marker="x", label="Observed contrast", zorder=4)
    ax.axhline(values["constant_baseline"].iat[0], color=".35", ls=":", lw=.8, label="Training mean")
    ax.set(xticks=x, xticklabels=genes, ylim=(0, observed.max()*1.1), ylabel="Joint recovery contrast",
           title="Skin transplantation  |  CD4+ TRB")
    ax.legend(ncol=2, frameon=False, loc="upper right")
    tidy(ax, "y")
    tag(ax, "a")
    for i, (chain, label) in enumerate(zip(CHAINS, "bc", strict=True)):
        ax = fig.add_subplot(grid[1, i])
        frame = pd.read_csv(spleen / f"{chain}_pooled" / "external_mouse_recovery.csv")
        order = frame.groupby("v_gene")["joint_recovery"].mean().sort_values(ascending=False).index
        x = np.arange(len(order))
        for j, (mouse, group) in enumerate(frame.groupby("mouse")):
            ax.scatter(x+(j-1)*.16, group.set_index("v_gene").loc[order, "joint_recovery"],
                       marker=("o", "s", "^")[j], c=(PURPLE, ORANGE, GREEN)[j], s=20, label=mouse)
        ax.set(xticks=x, xticklabels=order, ylabel="Joint recovery index",
               title=f"Marrow-graft recipients  |  {chain}")
        plt.setp(ax.get_xticklabels(), rotation=65 if chain == "TRA" else 40,
                 ha="right", rotation_mode="anchor", fontsize=7.5)
        ax.set_ylim(bottom=-.18)
        if i == 0:
            ax.legend(loc="upper right", ncol=3, frameon=False, fontsize=7.3)
        tidy(ax, "y")
        tag(ax, label)
    return fig


def plot_figures(data, spleen, out):
    sets = profiles(data)
    staging = out.parent / "panels"
    staging.mkdir(exist_ok=True)
    for chain in CHAINS:
        root = data / chain
        plot = plot_depth_check(pd.read_csv(root / "depth.csv"), chain, 15000)
        plot._suptitle.set_visible(False)
        for ax, letter in zip(plot.axes, "ab", strict=True):
            ax.set_title("")
            tag(ax, letter)
        save(plot, out / ("figure_2.pdf" if chain == "TRA" else "figure_S1.pdf"))
        mds = plot_mds(pd.read_csv(root / "mds.csv"), chain, 15000)
        mds._suptitle.set_visible(False)
        save(mds, out / ("figure_3.pdf" if chain == "TRA" else "figure_S2.pdf"))
        pooled = sets[(chain, "pooled")]
        threshold = cutoff(data, chain, "pooled")
        genes = eligible_segments(pooled, threshold).sort_values("score_s", ascending=False)["v"].tolist()
        for name, graphic in (
            ("dumbbell", plot_dumbbell(pooled, chain, thresholds=threshold)),
            ("localization", plot_localization(pooled, genes, chain, thresholds=threshold)),
        ):
            if graphic._suptitle is not None:
                graphic._suptitle.set_visible(False)
            graphic.set_size_inches(7.25, 3.20)
            graphic.axes[0].set_title(chain)
            tag(graphic.axes[0], "a" if chain == "TRA" else "b")
            save(graphic, staging / f"{chain}_{name}.pdf")
        comparison = pd.read_csv(root / "edgeR_fisher.csv")
        fisher = plot_fisher(comparison, chain)
        fisher.axes[0].set_title("")
        save(fisher, out / ("figure_9.pdf" if chain == "TRA" else "figure_S5.pdf"))
        subsets = plot_subsets({LABELS[s]: sets[(chain, s)] for s in ("cd4", "cd8")}, chain,
                               thresholds_by_subset={LABELS[s]: cutoff(data, chain, s) for s in ("cd4", "cd8")})
        subsets._suptitle.set_visible(False)
        save(subsets, out / ("figure_10.pdf" if chain == "TRA" else "figure_S6.pdf"))
        recovery = plot_recovery(pd.read_csv(root / "rarefaction.csv"), chain, "", selected_seed=DISPLAY_SEED)
        matrix = pd.read_csv(root / "rarefaction_regions.csv", index_col=0)
        summary = pd.read_csv(root / "rarefaction_regions_summary.csv", index_col=0)
        stability = plot_regional_stability(matrix, summary, DISPLAY_SEED, chain)
        recovery.set_size_inches(3.55, 2.65)
        stability.set_size_inches(3.55, 2.65)
        recovery.axes[0].set_title(f"{chain}  shortlist recovery")
        stability.axes[0].set_title("")
        legend = stability.axes[0].legend_
        if legend is not None:
            legend.remove()
        tag(recovery.axes[0], "a" if chain == "TRA" else "c")
        tag(stability.axes[0], "b" if chain == "TRA" else "d")
        save(recovery, staging / f"{chain}_stability_a.pdf")
        save(stability, staging / f"{chain}_stability_b.pdf")
        copyfile(root / "overlap.pdf", out / ("figure_4.pdf" if chain == "TRA" else "figure_S3.pdf"))
    for n, name in ((6, "dumbbell"), (7, "localization")):
        mount_panels([staging / f"{c}_{name}.pdf" for c in CHAINS], out / f"figure_{n}.pdf")
    mount_panels([staging / f"{c}_stability_{s}.pdf" for c in CHAINS for s in "ab"], out / "figure_11.pdf", columns=2)
    for n, graphic in (
        (5, retention_figure(data, sets)), (8, abundance_plot(data)), (12, selection_plot(data)),
        (13, external_structure(spleen)), (14, ranking_plot(spleen)), (15, error_plot(spleen)),
        (16, eligibility_plot(data, sets)), (17, prediction_plot(spleen)), ("S4", threshold_sensitivity(data)),
    ):
        save(graphic, out / f"figure_{n}.pdf")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--spleen-models", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    style()
    plot_figures(args.data, args.spleen_models, args.output)


if __name__ == "__main__":
    main()
