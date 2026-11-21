"""F2 and classical MDS for the g1/g2 depth check."""

from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from scipy.sparse import csc_matrix, csr_matrix
from scipy.spatial import procrustes
from scipy.stats import spearmanr

from .article_figures import article_canvas, empty_axis, _grid, _ticks, _title
from .figures import CORAL, PURPLE


DEPTHS = (10_000, 13_000, 15_000, 20_000, 30_000)
DEPTH_SEEDS = tuple(range(1201, 1251))
MDS_COLORS = {
    ("g1", "cd4"): CORAL,
    ("g1", "cd8"): "#bf2471",
    ("g2", "cd4"): PURPLE,
    ("g2", "cd8"): "#350a6d",
}
MDS_MARKERS = {"spleen": "o", "thymus": "^"}


def f2_matrix(counts):
    """Return the Bhattacharyya coefficient of the per-library aaV frequencies."""
    rooted = csc_matrix(np.sqrt(counts.to_numpy(dtype=float)))
    totals = counts.sum(axis=0).to_numpy()
    weights = np.divide(1, np.sqrt(totals), out=np.zeros(len(totals)), where=totals > 0)
    rooted = rooted.multiply(weights)
    f2 = np.clip((rooted.T @ rooted).toarray(), 0, 1)
    np.fill_diagonal(f2, 1)
    return f2


def classical_mds(f2):
    """Use the historical -log10(F2) distance and two classical MDS axes."""
    if not len(f2):
        return np.empty((0, 2))
    positive = f2[(f2 > 0) & (f2 < 1)]
    floor = positive.min() / 2 if len(positive) else 1e-12
    squared = np.square(-np.log10(np.maximum(f2, floor)))
    gram = -0.5 * (squared - squared.mean(axis=0) - squared.mean(axis=1)[:, None] + squared.mean())
    values, vectors = np.linalg.eigh(gram)
    coords = vectors[:, -2:][:, ::-1] * np.sqrt(np.maximum(values[-2:][::-1], 0))
    return np.pad(coords, ((0, 0), (0, 2 - coords.shape[1])))


def depth_check(counts, sample_seed, *, depths=DEPTHS, seeds=DEPTH_SEEDS):
    """Compare subsamples on a common set of libraries with the full-depth data."""
    counts = counts.loc[:, counts.sum(axis=0).ge(max(depths))]
    columns = ["depth", "seed", "n_libraries", "f2_rank_correlation", "mds_shape_discrepancy"]
    if counts.shape[1] < 3:
        return pd.DataFrame(columns=columns)
    reference = f2_matrix(counts)
    reference_mds = classical_mds(reference)
    upper = np.triu_indices_from(reference, k=1)
    reference_ranks = reference[upper]
    libraries = []
    for sample_id, counts_by_clone in counts.items():
        positive = np.flatnonzero(counts_by_clone.to_numpy())
        libraries.append((sample_id, positive, counts_by_clone.to_numpy()[positive]))
    rows = []
    for depth in depths:
        for seed in seeds:
            values, indices, offsets = [], [], [0]
            for sample_id, positive, counts_by_clone in libraries:
                rng = np.random.default_rng(sample_seed(seed, sample_id))
                drawn = rng.multivariate_hypergeometric(counts_by_clone, depth)
                retained = drawn > 0
                values.append(np.sqrt(drawn[retained] / depth))
                indices.append(positive[retained])
                offsets.append(offsets[-1] + int(retained.sum()))
            rooted = csr_matrix(
                (np.concatenate(values), np.concatenate(indices), offsets),
                shape=(len(libraries), len(counts)),
            )
            f2 = np.clip((rooted @ rooted.T).toarray(), 0, 1)
            np.fill_diagonal(f2, 1)
            ranks = f2[upper]
            correlation = np.nan
            if np.ptp(ranks) and np.ptp(reference_ranks):
                correlation = float(spearmanr(ranks, reference_ranks).statistic)
            coords = classical_mds(f2)
            discrepancy = np.nan
            if np.linalg.norm(reference_mds) and np.linalg.norm(coords):
                discrepancy = procrustes(reference_mds, coords)[2]
            rows.append((depth, seed, counts.shape[1], correlation, discrepancy))
    return pd.DataFrame(rows, columns=columns)


def plot_depth_check(frame, chain, target_depth):
    fig = article_canvas("depth")
    fig.suptitle(f"{chain} g1: depth check relative to g2 controls", fontsize=11.5, y=0.981)
    runs = frame["seed"].nunique() or len(DEPTH_SEEDS)
    fig.text(0.5, 0.915, f"{runs} independent subsamples at each depth; medians and 5th–95th percentiles",
             ha="center", fontsize=9)
    for i, (column, title, ylabel) in enumerate((
        ("f2_rank_correlation", "A. Reproducibility of F2 ranks", "F2 rank correlation"),
        ("mds_shape_discrepancy", "B. Stability of MDS", "MDS shape discrepancy"),
    )):
        ax = fig.add_axes((0.075 + 0.495 * i, 0.114, 0.423, 0.71))
        _ticks(ax, 8)
        _grid(ax, "both")
        _title(ax, title, size=10)
        ax.set_xlabel("Functional UMI per library (thousands)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        rows = frame.dropna(subset=[column])
        target = target_depth / 1000
        ax.axvline(target, color=PURPLE, linestyle=":", linewidth=0.85)
        if rows.empty:
            empty_axis(ax, "Insufficient comparable libraries at these depths.")
            continue
        grouped = rows.groupby("depth")[column]
        median = grouped.median()
        x = median.index.to_numpy() / 1000
        ax.errorbar(x, median, yerr=[median - grouped.quantile(0.05), grouped.quantile(0.95) - median],
                    color=CORAL, marker="o", markersize=4, linewidth=1.25, ecolor=PURPLE,
                    elinewidth=0.75, capsize=2)
        ax.set_xticks(x)
        if target_depth in median:
            value = median.loc[target_depth]
            ax.scatter(target, value, s=26, color=PURPLE, zorder=4)
            ax.annotate(f"{target_depth:,} UMI", (target, value), xytext=(25, -29 if i == 0 else 23),
                        textcoords="offset points", ha="left", fontsize=7,
                        arrowprops={"arrowstyle": "-", "color": "black", "lw": 0.6})
    return fig


def mds_points(counts, metadata):
    md = metadata.set_index("sample_id").reindex(counts.columns).copy()
    md["MDS1"], md["MDS2"] = classical_mds(f2_matrix(counts)).T
    # Orient the axes by the g1 groups, independent of eigensolver signs.
    g1_cd4 = md["group_no"].eq("g1") & md["_cell_subset"].eq("cd4")
    if g1_cd4.any() and not g1_cd4.all():
        if md.loc[g1_cd4, "MDS1"].mean() > md.loc[~g1_cd4, "MDS1"].mean():
            md["MDS1"] *= -1
    thymus = md["_tissue"].eq("thymus")
    if (g1_cd4 & thymus).any() and (g1_cd4 & ~thymus).any():
        if md.loc[g1_cd4 & thymus, "MDS2"].mean() < md.loc[g1_cd4 & ~thymus, "MDS2"].mean():
            md["MDS2"] *= -1
    return md


def plot_mds(points, chain, target_depth):
    fig = article_canvas("mds")
    fig.suptitle(f"{chain} g1: MDS relative to g2 control repertoires", fontsize=12, y=0.985)
    fig.text(0.5, 0.940, f"{target_depth:,} functional UMI per library", ha="center", fontsize=9.5)
    ax = fig.add_axes((0.072, 0.139, 0.750, 0.748))
    _ticks(ax, 9, box=True)
    fig.text(0.45, 0.055, "MDS1", ha="center", fontsize=11)
    fig.text(0.016, 0.51, "MDS2", va="center", rotation=90, fontsize=11)
    if len(points) < 2:
        empty_axis(ax, "Insufficient g1/g2 libraries for MDS.")
        return fig

    x = np.sort(points["MDS1"].to_numpy())
    gaps = np.diff(x)
    largest = gaps.argmax()
    span = x[-1] - x[0]
    margin = max(0.05 * span, 0.015)
    split = len(x) >= 4 and span > 0 and gaps[largest] > 0.34 * span and 0 < largest < len(x) - 2
    if split:
        left_end, right_start = x[largest:largest + 2]
        margin = min(margin, gaps[largest] / 4)
        left_span = left_end - x[0] + 2 * margin
        right_span = x[-1] - right_start + 2 * margin
        width = 0.73 * left_span / (left_span + right_span)
        ax.set_position((0.072, 0.139, width, 0.748))
        right = fig.add_axes((0.092 + width, 0.139, 0.73 - width, 0.748), sharey=ax)
        _ticks(right, 9, box=True)
        ax.set_xlim(x[0] - margin, left_end + margin)
        right.set_xlim(right_start - margin, x[-1] + margin)
        ax.spines["right"].set_visible(False)
        right.spines["left"].set_visible(False)
        right.tick_params(axis="y", left=False, labelleft=False)
        for axis, side in ((ax, 1), (right, 0)):
            axis.plot([side, side], [0, 1], linestyle="none", marker=[(-1, -0.6), (1, 0.6)],
                      markersize=9, markeredgewidth=0.85, color="black",
                      transform=axis.transAxes, clip_on=False)
        fig.text(0.45, 0.022, "Broken X axis; horizontal distances across the gap are not to scale.",
                 ha="center", fontsize=9)
        axes = (ax, right)
    else:
        ax.set_xlim(x[0] - margin, x[-1] + margin)
        axes = (ax,)

    handles = []
    for axis in axes:
        _grid(axis, "both")
        axis.xaxis.set_major_locator(MaxNLocator(max(3, int(9 * axis.get_position().width / 0.75))))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    for (group, subset), color in MDS_COLORS.items():
        for tissue, marker in MDS_MARKERS.items():
            rows = points.loc[points["group_no"].eq(group) & points["_cell_subset"].eq(subset)
                              & points["_tissue"].eq(tissue)]
            if rows.empty:
                continue
            handles.append(Line2D([], [], linestyle="none", marker=marker, markersize=5, color=color,
                                  label=f"{group} / {subset.upper()} / {tissue}"))
            for axis in axes:
                lower, upper = axis.get_xlim()
                visible = rows.loc[rows["MDS1"].between(lower, upper)]
                axis.scatter(visible["MDS1"], visible["MDS2"], marker=marker, s=27, color=color,
                             linewidth=0, alpha=0.9)
    ax.margins(y=0.09)
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.830, 0.885),
               fontsize=7.5, frameon=True, borderpad=0.4, handletextpad=0.45)
    return fig
