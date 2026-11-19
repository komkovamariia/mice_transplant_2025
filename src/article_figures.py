"""Approach 1 figures for the manuscript."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from .v_thresholds import FIXED_THRESHOLDS
from .figures import CORAL, MAGENTA, PURPLE, THRESHOLD_BLUE


# Canvas sizes follow the manuscript's PNG exports.
CANVAS_PIXELS = {
    "depth": (6138, 2349),
    "mds": (5594, 3687),
    "overlap": (6141, 4191),
    "dumbbell": (6138, 5100),
    "retention": (5040, 3900),
    "localization": (5847, 3888),
    "edger": (6138, 2913),
    "fisher": (3762, 3762),
    "spike_in": (6438, 4230),
    "cd4_cd8": (6870, 5725),
}
ARTICLE_FIGURES = {
    "depth": "depth_check",
    "mds": "mds",
    "overlap": "clonotype_overlap",
    "dumbbell": "v_segment_representation",
    "retention": "v_segment_retention",
    "localization": "v_segment_localization",
    "edger": "edger_summary",
    "fisher": "fisher_edger_concordance",
    "spike_in": "spike_in_sensitivity",
    "cd4_cd8": "cd4_cd8_structure",
}
ARTICLE_FONTS = ["Arial", "Liberation Sans", "DejaVu Sans"]
REGIONS = {
    "100": "g1 only",
    "110": "g1 ∩ g5\nwithout g6",
    "101": "g1 ∩ g6\nwithout g5",
    "111": "g1 ∩ g5 ∩ g6",
}
SHARE_COLUMN = "frequency_in_full_g1_pct"
RETENTION_COLUMN = "retained_cdr3_share_pct"
COUNT_COLUMN = "cdr3_in_full_g1"
SCORE_COLUMNS = ["score_s", SHARE_COLUMN, RETENTION_COLUMN, "v"]


def article_name(chain, panel):
    return f"{chain}_g1_{ARTICLE_FIGURES[panel]}.png"


def article_canvas(panel):
    dpi = 300 if panel == "spike_in" else 600
    return plt.figure(figsize=tuple(np.nextafter(pixels / dpi, np.inf) for pixels in CANVAS_PIXELS[panel]))


def panel_bbox(ax):
    fig = ax.figure
    fig.canvas.draw()
    bbox = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    return bbox.padded(0.06)


def article_dpi(filename):
    return 300 if filename.endswith("_spike_in_sensitivity.png") else 600


def empty_axis(ax, message):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(0.5, 0.5, message, ha="center", va="center", transform=ax.transAxes, wrap=True)


def panel_letter(ax, label, x=-0.10, y=1.10, fontsize=13):
    ax.text(x, y, label, fontsize=fontsize, fontweight="bold", transform=ax.transAxes)


def _ticks(ax, size, *, box=False):
    ax.tick_params(axis="both", labelsize=size, length=3, width=0.65, pad=3)
    for spine in ax.spines.values():
        spine.set_visible(box)
        spine.set_linewidth(0.65)
        spine.set_color("#333333")
    if not box:
        for side in ("left", "bottom"):
            ax.spines[side].set_visible(True)
    ax.grid(False)


def _grid(ax, axis):
    ax.grid(axis=axis, color="#e6e6e6", linewidth=0.6, alpha=0.38)
    ax.set_axisbelow(True)


def _title(ax, title, size=13, pad=5):
    ax.set_title(title, fontsize=size, fontweight="normal", pad=pad)


def eligible_segments(profile, thresholds=FIXED_THRESHOLDS):
    return profile.loc[
        profile[SHARE_COLUMN].ge(thresholds.initial_share_pct)
        & profile[RETENTION_COLUMN].ge(thresholds.retention_pct)
    ]


def plot_dumbbell(profile, chain, suffix="", *, thresholds=FIXED_THRESHOLDS):
    rows = eligible_segments(profile, thresholds).sort_values(
        SCORE_COLUMNS, ascending=[False, False, False, True],
    ).head(20)
    fig = article_canvas("dumbbell")
    ax = fig.add_axes((0.169, 0.071, 0.817, 0.879))
    _ticks(ax, 10)
    _grid(ax, "x")
    fig.suptitle(f"{chain} g1: V-segment representation before and after background subtraction{suffix}",
                 fontsize=13, y=0.988)
    ax.set_xlabel("Share of aaV clonotypes with the V segment in g1 (%)", fontsize=10)
    ax.set_ylabel(f"{chain} V segment", fontsize=10)
    ax.yaxis.set_label_coords(-0.172, 0.5)
    if rows.empty:
        empty_axis(ax, "No V segments met both selection thresholds.")
        return fig
    y = np.arange(len(rows))
    retained = rows["frequency_remaining_in_g1_pct_relative_to_original"]
    ax.hlines(y, retained, rows[SHARE_COLUMN], color="#d47dab", linewidth=1.25)
    ax.scatter(retained, y, s=42, color=PURPLE, label="After background subtraction", zorder=3)
    ax.scatter(rows[SHARE_COLUMN], y, s=42, color=CORAL, label="Complete g1 repertoire", zorder=3)
    ax.set_yticks(y, rows["v"])
    ax.invert_yaxis()
    ax.legend(frameon=False, loc="upper left", fontsize=9)
    return fig


def bubble_candidates(profile, limit, thresholds=FIXED_THRESHOLDS):
    plotted = profile.dropna(subset=SCORE_COLUMNS + [COUNT_COLUMN, "region_entropy"])
    plotted = plotted.loc[plotted[COUNT_COLUMN].gt(0) & plotted[SHARE_COLUMN].gt(0)]
    eligible = eligible_segments(plotted, thresholds)
    selected = eligible.sort_values(SCORE_COLUMNS, ascending=[False, False, False, True]).head(limit)
    return plotted, selected


def _retention_panel(fig, ax, color_axis, profile, title, *, compact=False, thresholds=FIXED_THRESHOLDS):
    from adjustText import adjust_text

    rows, selected = bubble_candidates(profile, 10 if compact else 11, thresholds)
    _ticks(ax, 7 if compact else 10, box=not compact)
    if compact:
        _ticks(color_axis, 6, box=True)
    _title(ax, title, size=9 if compact else 13, pad=4)
    if compact:
        ax.title.set_fontweight("bold")
    size = 7 if compact else 12
    ax.set_xlabel(
        "V-segment share in the complete group g1, %" if compact else
        "Initial share of unique aaV clonotypes in g1 (%)",
        fontsize=size, labelpad=6,
    )
    ax.set_ylabel("Retention after subtraction (%)", fontsize=size, labelpad=3)
    ax.axvline(thresholds.initial_share_pct, color=THRESHOLD_BLUE, linestyle=":", linewidth=0.70, alpha=0.87)
    ax.axhline(thresholds.retention_pct, color=THRESHOLD_BLUE, linestyle=":", linewidth=0.70, alpha=0.87)
    if rows.empty:
        empty_axis(ax, "No V segments were available.")
        color_axis.set_axis_off()
        return selected["v"].tolist()

    points = ax.scatter(
        rows[SHARE_COLUMN], rows[RETENTION_COLUMN], c=rows["region_entropy"],
        cmap="magma", vmin=0, vmax=max(0.84 if compact else 0.62, rows["region_entropy"].max()),
        s=(22 + 302 * rows[COUNT_COLUMN] / rows[COUNT_COLUMN].max()) * (0.60 if compact else 1),
        alpha=0.74, edgecolors="face", linewidths=0.35, zorder=2,
    )
    ax.set_xlim(0, max(thresholds.initial_share_pct * 1.2, rows[SHARE_COLUMN].max() * 1.12))
    ax.set_ylim(max(0, min(13, thresholds.retention_pct, rows[RETENTION_COLUMN].min() - 4)), 108)
    colorbar = fig.colorbar(points, cax=color_axis)
    colorbar.solids.set_alpha(1)
    colorbar.outline.set_linewidth(0.65)
    colorbar.set_label(
        "Normalized regional entropy" if not compact else "Normalized Shannon entropy (region)",
        fontsize=size, labelpad=7,
    )
    colorbar.ax.tick_params(labelsize=6 if compact else 10, width=0.65, length=3)

    texts = []
    for row in selected.itertuples(index=False):
        x, y = getattr(row, SHARE_COLUMN), getattr(row, RETENTION_COLUMN)
        texts.append(ax.text(x, y + 1.8, row.v, ha="center", va="bottom", fontsize=6 if compact else 9.5))
    if texts:
        adjust_text(
            texts, x=rows[SHARE_COLUMN].to_numpy(), y=rows[RETENTION_COLUMN].to_numpy(),
            target_x=selected[SHARE_COLUMN].to_numpy(), target_y=selected[RETENTION_COLUMN].to_numpy(),
            expand=(1.12, 1.25), force_static=(0.3, 0.35), force_text=(0.55, 0.7),
            iter_lim=300, ensure_inside_axes=True, ax=ax,
        )
        for label, row in zip(texts, selected.itertuples(index=False), strict=True):
            diameter = np.sqrt((22 + 302 * getattr(row, COUNT_COLUMN) / rows[COUNT_COLUMN].max())
                               * (0.60 if compact else 1))
            ax.annotate(
                "", xy=(getattr(row, SHARE_COLUMN), getattr(row, RETENTION_COLUMN)),
                xytext=label.get_position(), zorder=2.5,
                arrowprops={"arrowstyle": "-", "color": "#333333", "lw": 0.65,
                            "patchA": label, "shrinkA": 0.8, "shrinkB": diameter / 2},
            )
    return selected["v"].tolist()


def plot_retention(profile, chain, stratum_title, *, thresholds=FIXED_THRESHOLDS):
    fig = article_canvas("retention")
    ax = fig.add_axes((0.074, 0.110, 0.757, 0.844))
    color_axis = fig.add_axes((0.883, 0.110, 0.032, 0.844))
    genes = _retention_panel(
        fig, ax, color_axis, profile, f"{chain} g1: {stratum_title} V-segment retention",
        thresholds=thresholds,
    )
    return fig, genes


def _localization_panel(fig, ax, color_axis, profile, genes, chain, title, *, compact=False, thresholds=FIXED_THRESHOLDS):
    size = 7 if compact else 10
    _ticks(ax, size)
    for spine in ax.spines.values():
        spine.set_visible(False)
    _title(ax, title, size=9 if compact else 13, pad=7)
    if compact:
        ax.title.set_fontweight("bold")
    ax.set_xlabel("g1-containing overlap region", fontsize=size, labelpad=6)
    ax.set_ylabel(f"{chain} V segment", fontsize=size, labelpad=7)
    if not compact:
        ax.yaxis.set_label_coords(-0.191, 0.5)
    rows = eligible_segments(profile.loc[profile["v"].isin(genes)], thresholds)
    rows = rows.reindex(columns=["v", "score_s", *REGIONS])
    rows = rows.fillna({column: 0.0 for column in REGIONS})
    if rows.empty:
        empty_axis(ax, "No V segments passed both retention thresholds.")
        color_axis.set_axis_off()
        return
    rows = rows.sort_values(["100", "score_s", "v"], ascending=[False, False, True])
    values = rows[list(REGIONS)].to_numpy(dtype=float)
    image = ax.imshow(values, aspect="auto", interpolation="nearest", cmap="magma", vmin=0, vmax=1)
    region_labels = list(REGIONS.values())
    if compact:
        region_labels = ["g1 only", "g1 ∩ g5", "g1 ∩ g6", "g1 ∩ g5 ∩ g6"]
    ax.set_xticks(range(len(REGIONS)), region_labels)
    ax.set_yticks(range(len(rows)), rows["v"])
    ax.tick_params(length=0, pad=4)
    if not compact:
        ax.tick_params(axis="x", pad=12)
        ax.xaxis.labelpad = 15
    if compact:
        plt.setp(ax.get_xticklabels(), rotation=37, ha="right", rotation_mode="anchor")
    for (row, col), value in np.ndenumerate(values):
        rgb = image.cmap(image.norm(value))[:3]
        luminance = np.dot(rgb, (0.2126, 0.7152, 0.0722))
        ax.text(col, row, f"{value:.2f}", ha="center", va="center",
                color="white" if luminance < 0.52 else "black", fontsize=6.8 if compact else 9)
    colorbar = fig.colorbar(image, cax=color_axis, ticks=np.linspace(0, 1, 6))
    colorbar.outline.set_linewidth(0.65)
    colorbar.ax.tick_params(labelsize=size - 1, length=3, width=0.65)
    if not compact:
        colorbar.set_label(
            "Fraction of each V segment's g1 UMI in the region",
            fontsize=size, labelpad=8,
        )


def plot_localization(profile, genes, chain, suffix="", *, thresholds=FIXED_THRESHOLDS):
    fig = article_canvas("localization")
    ax = fig.add_axes((0.156, 0.148, 0.738, 0.794))
    color_axis = fig.add_axes((0.916, 0.148, 0.027, 0.794))
    title = f"{chain} g1: V-segment localization across g1-containing overlap regions{suffix}"
    _localization_panel(fig, ax, color_axis, profile, genes, chain, title, thresholds=thresholds)
    return fig


def plot_subsets(profiles, chain, *, thresholds_by_subset=None):
    fig = article_canvas("cd4_cd8")
    fig.suptitle(f"{chain} g1: separate CD8+ and CD4+ structural analyses",
                 fontsize=15, fontweight="bold", y=0.987)
    for row, subset in enumerate(("CD8+", "CD4+")):
        bottom = 0.546 - row * 0.473
        bubble = fig.add_axes((0.041, bottom - 0.036, 0.403, 0.390))
        bubble_color = fig.add_axes((0.450, bottom - 0.036, 0.016, 0.390))
        heatmap = fig.add_axes((0.605, bottom, 0.319, 0.354))
        heatmap_color = fig.add_axes((0.939, bottom, 0.016, 0.354))
        profile = profiles[subset]
        thresholds = (thresholds_by_subset or {}).get(subset, FIXED_THRESHOLDS)
        genes = _retention_panel(fig, bubble, bubble_color, profile,
                                 f"{subset} V-segment retention", compact=True, thresholds=thresholds)
        _localization_panel(fig, heatmap, heatmap_color, profile, genes, chain,
                            f"{subset} overlap-region composition", compact=True, thresholds=thresholds)
        panel_letter(bubble, "AC"[row], x=-0.067, y=1.09, fontsize=15)
        panel_letter(heatmap, "BD"[row], x=-0.275, y=1.10, fontsize=15)
    return fig


def plot_edger(ranking, chain, suffix=""):
    fig = article_canvas("edger")
    fig.suptitle(f"{chain} g1: edgeR V-segment differential-abundance summaries{suffix}",
                 fontsize=12, y=0.979)
    effect = "mean_edger_log2fc"
    count = "n_edger_significant"
    tested = ranking.loc[ranking[count].gt(0)]
    for i, (column, title, xlabel, color) in enumerate((
        (effect, "A. V segments enriched in g1", "Mean log2 fold change among significant aaV clonotypes", CORAL),
        (count, "B. V segments with the largest significant aaV support", "Number of significant aaV clonotypes", PURPLE),
    )):
        ax = fig.add_axes((0.084 + 0.497 * i, 0.085, 0.407, 0.787))
        _ticks(ax, 8)
        _grid(ax, "x")
        _title(ax, title, size=10, pad=10)
        ax.set_xlabel(xlabel, fontsize=9)
        if i == 0:
            ax.set_ylabel(f"{chain} V segment", fontsize=9)
        rows = tested.dropna(subset=[column]).sort_values([column, "v_gene"], ascending=[False, True]).head(15)
        if rows.empty:
            empty_axis(ax, "No edgeR-significant aaV clonotypes.")
        else:
            ax.barh(rows["v_gene"], rows[column], color=color, height=0.78)
            ax.invert_yaxis()
            if i == 1:
                ax.xaxis.set_major_locator(MaxNLocator(9, integer=True))
    return fig


def plot_fisher(comparison, chain, suffix=""):
    fig = article_canvas("fisher")
    ax = fig.add_axes((0.107, 0.100, 0.878, 0.826))
    _ticks(ax, 10)
    _title(ax, f"{chain} g1: Fisher and edgeR effect concordance{suffix}", size=14, pad=10)
    ax.set_xlabel("Fisher pooled-count log₂ ratio: g1 relative to g5+g6", fontsize=11)
    ax.set_ylabel("edgeR log₂ fold change: g1 relative to g5+g6", fontsize=11)
    columns = ["Fisher_log2FC_g1", "logFC"]
    rows = comparison[columns].replace([np.inf, -np.inf], np.nan).dropna()
    if rows.empty:
        empty_axis(ax, "No aaV clonotype was jointly available from Fisher and edgeR.")
        return fig
    x, y = rows.to_numpy().T
    lower, upper = rows.min().min(), rows.max().max()
    ax.plot([lower, upper], [lower, upper], "--", color="black", linewidth=0.85, zorder=1)
    ax.scatter(x, y, s=14, alpha=0.88, color=MAGENTA, linewidth=0, zorder=2)
    anchor = np.argmin((x - 0.8)**2 + (y - 2.8)**2)
    ax.annotate(
        "Each point: one aaV clonotype\n"
        "Dashed diagonal: equal Fisher and edgeR effect\n"
        "Distance from the line reflects method-specific scaling",
        xy=(x[anchor], y[anchor]), xytext=(0.035, 0.955),
        textcoords="axes fraction", ha="left", va="top", fontsize=8.7, linespacing=1.4,
        bbox={"facecolor": "white", "edgecolor": "black", "linewidth": 0.65,
              "boxstyle": "square,pad=0.55"},
        arrowprops={"arrowstyle": "-", "color": "black", "lw": 0.7, "relpos": (0.64, 0)},
    )
    return fig


def plot_spike_in(frame, chain):
    fig = article_canvas("spike_in")
    gene = frame.iloc[0]["v_gene"]
    family_size = int(frame.iloc[0]["n_selected"])
    fig.suptitle(f"{chain} g1: {gene} spike-in sensitivity", fontsize=26, y=0.995)
    doses = frame["fraction"].to_numpy()
    positions = np.arange(len(frame))
    tick_positions = [0, *positions[1::2]]
    if tick_positions[-1] != positions[-1]:
        tick_positions.append(positions[-1])
    labels = ["baseline" if dose == 0 else f"{100 * dose:g}%" for dose in doses]
    exact = frame["all_selected_in_remainder"]
    significant = frame["all_edger_positive_fdr"]
    positive_doses = frame["fraction"].gt(0)
    exact_passes = frame.index[exact & positive_doses]
    edger_passes = frame.index[significant & positive_doses]
    both = frame.index[exact & significant & positive_doses]
    common = int(both[0]) if len(both) else None
    panels = (
        ("recovery_fraction", "Exact aaV recovery", "Recovered aaV fraction"),
        ("edger_positive_recovery_fraction", "edgeR recovery", "Significant aaV fraction"),
        ("v_umi_share_g1", f"{gene} UMI share", "UMI share in g1 (%)"),
        ("v_final_rank", "Integrated V-segment rank", "Integrated rank"),
    )
    arrow = {"arrowstyle": "-", "color": "#222222", "lw": 0.85}
    for i, (column, title, ylabel) in enumerate(panels):
        ax = fig.add_axes((0.053 + 0.502 * (i % 2), 0.596 - 0.449 * (i // 2), 0.438, 0.344))
        _ticks(ax, 14)
        _title(ax, title, size=20, pad=4)
        panel_letter(ax, "ABCD"[i], x=-0.10, y=1.015, fontsize=22)
        y = frame[column] * (100 if i == 2 else 1)
        ax.plot(positions, y, "o-", color=CORAL, linewidth=1.9, markersize=8)
        ax.set_xticks(tick_positions, [labels[tick] for tick in tick_positions], rotation=45, ha="right")
        ax.set_xlabel("Added family UMI relative to the original g1 UMI library", fontsize=15)
        ax.set_ylabel(ylabel, fontsize=15)
        if common is not None:
            ax.axvline(common, color=PURPLE, linestyle=(0, (3, 3)), linewidth=1.25)
            if pd.notna(y.iloc[common]):
                ax.scatter(common, y.iloc[common], s=170, facecolors="none",
                           edgecolors=PURPLE, linewidths=1.4, zorder=4)
                if i < 2:
                    ax.scatter(common, y.iloc[common], s=61, color=PURPLE, zorder=3)
        if i < 2:
            ax.set_ylim(-0.05, 1.05)
            passed = exact_passes if i == 0 else edger_passes
            if i == 1 and len(exact_passes) and (not len(passed) or exact_passes[0] < passed[0]):
                first_exact = int(exact_passes[0])
                ax.scatter(first_exact, y.iloc[first_exact], s=61, color=PURPLE, zorder=3)
            if len(passed):
                first = int(passed[0])
                ax.scatter(first, y.iloc[first], s=61, color=PURPLE, zorder=3)
                ax.annotate(f"all {family_size} {'recovered' if i == 0 else 'significant'}\nfrom {labels[first]}",
                            (first, y.iloc[first]), xytext=(0.58, 0.86), textcoords="axes fraction",
                            fontsize=16, va="top", arrowprops=arrow)
        elif i == 2 and common is not None:
            ax.annotate(f"At {labels[common]}:\nUMI share is {y.iloc[common]:.3g}%",
                        (common, y.iloc[common]), xytext=(0.59, 0.74), textcoords="axes fraction",
                        fontsize=16, arrowprops=arrow)
        elif i == 3:
            ax.invert_yaxis()
            changes = y.notna() & y.ne(y.iloc[0])
            if pd.notna(y.iloc[0]) and changes.any():
                first = int(np.flatnonzero(changes)[0])
                ax.annotate(
                    f"Rank remains {y.iloc[0]:g} through {labels[first - 1]};\n"
                    f"changes to {y.iloc[first]:g} at {labels[first]}",
                    (first, y.iloc[first]), xytext=(0.23, 0.75),
                    textcoords="axes fraction", ha="center", fontsize=16, arrowprops=arrow,
                )
    if common is not None:
        note = f"Common detection threshold: {labels[common]} added family UMI; all {family_size} aaV recovered and significant"
    else:
        note = f"No tested dose achieved recovery and significance for all {family_size} aaV."
    fig.text(0.5, 0.032, note, ha="center", fontsize=10.5)
    return fig
