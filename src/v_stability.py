"""Sampling stability for the complete adaptive V-segment shortlist."""

import numpy as np
import pandas as pd

from .spike_in import SUBTRACTION_GROUPS
from .structural_sampling import downsample_group_tables
from .structural_v import structural_summary, v_region_proportions
from .v_thresholds import fit_plot_thresholds


RAREFACTION_SEEDS = tuple(range(1201, 1301))
BALANCE_SEEDS = tuple(range(2201, 2251))


def adaptive_profile(tables, *, include_regions=False):
    profile = structural_summary(tables)
    share, retention = 100 * profile["initial_share"], 100 * profile["retention"]
    cutoffs = fit_plot_thresholds(share, retention)
    profile["eligible"] = cutoffs.includes(share, retention) if cutoffs else False
    if include_regions:
        profile = profile.join(v_region_proportions(tables)["100"].rename("only_g1_share"), on="v")
    return profile.sort_values(["score_s", "v"], ascending=[False, True]), cutoffs


def rarefaction_profiles(tables, depth, seeds=RAREFACTION_SEEDS):
    profiles, shortlists = [], {}
    for seed in seeds:
        sampled, _ = downsample_group_tables(tables, depth=depth, seed=seed)
        profile, cutoffs = adaptive_profile(sampled, include_regions=True)
        profiles.append(profile.assign(seed=seed))
        shortlists[seed] = (profile.loc[profile["eligible"], "v"].tolist(), cutoffs)
    return pd.concat(profiles, ignore_index=True), shortlists


def representative_seed(shortlists):
    valid = {seed: set(genes) for seed, (genes, cutoffs) in shortlists.items() if cutoffs is not None}
    if not valid:
        return None, pd.DataFrame(columns=["seed", "mean_jaccard_distance", "selected_representative_seed"])
    distances = [
        {
            "seed": seed,
            "mean_jaccard_distance": np.mean([
                1 - len(genes & other) / len(genes | other) if genes | other else 0.0
                for other in valid.values()
            ]),
        }
        for seed, genes in valid.items()
    ]
    summary = pd.DataFrame(distances).sort_values(["mean_jaccard_distance", "seed"])
    selected = int(summary["seed"].iat[0])
    return selected, summary.assign(selected_representative_seed=summary["seed"].eq(selected))


def compare_shortlists(reference, shortlists):
    reference_genes, reference_cutoffs = reference
    anchor = set(reference_genes)
    rows = []
    for seed, (genes, cutoffs) in shortlists.items():
        observed = set(genes)
        overlap = len(anchor & observed)
        union = len(anchor | observed)
        comparable = reference_cutoffs is not None and cutoffs is not None
        rows.append({
            "seed": seed,
            "status": "compared" if comparable else "no_valid_cutoffs",
            "initial_share_pct": cutoffs.initial_share_pct if cutoffs else np.nan,
            "retention_pct": cutoffs.retention_pct if cutoffs else np.nan,
            "n_reference": len(anchor) if reference_cutoffs else np.nan,
            "n_replicate": len(observed) if cutoffs else np.nan,
            "n_recovered": overlap if comparable else np.nan,
            "reference_recovery": overlap / len(anchor) if comparable and anchor else np.nan,
            "jaccard": overlap / union if comparable and union else np.nan,
            "reference_segments": ";".join(reference_genes),
            "replicate_segments": ";".join(genes),
        })
    return pd.DataFrame(rows)


def balanced_shortlists(tables, metadata, seeds=BALANCE_SEEDS):
    groups = ("g1", *SUBTRACTION_GROUPS)
    mouse_by_sample = metadata.drop_duplicates("sample_id").set_index("sample_id")["mouse_id"]
    sample_mice = {group: mouse_by_sample.loc[tables[group].columns] for group in groups}
    n_mice = min(mice.nunique() for mice in sample_mice.values())
    if n_mice < 2:
        return {}, n_mice
    shortlists = {}
    for seed in seeds:
        rng = np.random.default_rng(seed)
        balanced = dict(tables)
        for group, mice in sample_mice.items():
            selection = rng.choice(sorted(mice.unique()), size=n_mice, replace=False)
            selected = tables[group].loc[:, mice.isin(selection).to_numpy()]
            balanced[group] = selected.loc[selected.sum(axis=1).gt(0)]
        profile, cutoffs = adaptive_profile(balanced)
        shortlists[seed] = (profile.loc[profile["eligible"], "v"].tolist(), cutoffs)
    return shortlists, n_mice


def plot_recovery(frame, chain, title, *, selected_seed=None):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    from .figures import CORAL, PURPLE

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    recovered = frame["n_recovered"].dropna().astype(int)
    if recovered.empty:
        ax.text(0.5, 0.5, "No comparable shortlists", ha="center", transform=ax.transAxes)
    else:
        size = int(frame["n_reference"].dropna().iat[0])
        counts = recovered.value_counts().reindex(range(size + 1), fill_value=0).sort_index()
        ax.bar(counts.index, counts, width=0.72, color=CORAL, alpha=0.88)
        median = recovered.median()
        ax.axvline(median, color=PURPLE, ls=":", lw=1.0, label=f"Median {median:g}/{size}")
        selected = frame.loc[frame["seed"].eq(selected_seed), "n_recovered"].dropna()
        if len(selected):
            recovered_at_seed = int(selected.iat[0])
            ax.scatter(recovered_at_seed, counts.loc[recovered_at_seed], s=26, facecolors="none",
                       edgecolors=PURPLE, linewidths=1.0, label=f"Seed {selected_seed}")
        ax.set_xlim(-0.6, size + 0.6)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(frameon=False, loc="upper left")
    ax.set(xlabel="Reference V segments recovered", ylabel="Number of draws", title=f"{chain} g1: {title}")
    ax.grid(False)
    fig.tight_layout()
    return fig


def regional_tables(profiles, genes, seed):
    matrix = profiles.pivot(index="seed", columns="v", values="only_g1_share").reindex(columns=genes) * 100
    summary = matrix.describe().T if genes else pd.DataFrame(columns=["mean"])
    summary["representative_seed_share_pct"] = matrix.loc[seed] if seed in matrix.index else np.nan
    return matrix, summary


def plot_regional_stability(matrix, summary, seed, chain):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from .figures import CORAL, PURPLE

    fig, ax = plt.subplots(figsize=(7.1, max(3.15, 1.4 + 0.26 * len(matrix.columns))), layout="constrained")
    if len(matrix.columns):
        genes = matrix.columns.tolist()
        positions = np.arange(1, len(genes) + 1)
        ax.boxplot([matrix[gene].dropna() for gene in genes], positions=positions, vert=False,
                   widths=0.50, patch_artist=True, boxprops={"facecolor": "#f5cdbe", "linewidth": 0.65},
                   medianprops={"color": PURPLE, "linewidth": 1.0},
                   whiskerprops={"linewidth": 0.65}, capprops={"linewidth": 0.65},
                   flierprops={"markersize": 2})
        mean = ax.scatter(summary.reindex(genes)["mean"], positions, s=22, color=CORAL, marker="D",
                          label="Mean", zorder=3)
        handles = [Patch(facecolor="#f5cdbe", edgecolor="black", label="IQR and median"), mean]
        chosen = summary.reindex(genes)["representative_seed_share_pct"]
        if chosen.notna().any():
            point = ax.scatter(chosen, positions, s=29, facecolors="none", edgecolors=PURPLE,
                               linewidths=1.0, label=f"Seed {seed}", zorder=4)
            handles.append(point)
        ax.set_yticks(positions, genes)
        ax.invert_yaxis()
        ax.legend(handles=handles, frameon=False, loc="lower center",
                  bbox_to_anchor=(0.5, 1.01), ncol=len(handles), columnspacing=1.2)
        ax.set_title(f"{chain} g1: V-segment overlap across rarefaction draws", pad=27)
    else:
        ax.text(0.5, 0.5, "No V segments passed both cutoffs", ha="center", transform=ax.transAxes)
        ax.set_title(f"{chain} g1: V-segment overlap across rarefaction draws")
    ax.set(xlabel="g1-only share of V-segment UMI (%)", ylabel="V segment")
    ax.grid(axis="x", color="#eadbed", linewidth=0.55)
    ax.set_axisbelow(True)
    return fig


def write_stability(raw, displayed, metadata, depth, seed, chain, results, save_figure, *, reruns=None):
    profiles, shortlists = rarefaction_profiles(raw, depth) if reruns is None else reruns
    reference_profile, cutoffs = adaptive_profile(displayed)
    genes = reference_profile.loc[reference_profile["eligible"], "v"].tolist()
    reference = genes, cutoffs
    prefix = chain.lower()
    profiles.to_csv(results / f"{prefix}_rarefaction_v_profiles.csv", index=False)

    rarefaction = compare_shortlists(reference, shortlists)
    rarefaction.to_csv(results / f"{prefix}_rarefaction_adaptive_stability.csv", index=False)
    save_figure(
        plot_recovery(rarefaction, chain, f"{depth:,}-UMI rarefaction stability", selected_seed=seed),
        f"{chain}_g1_rarefaction_stability.png",
    )

    balanced, n_mice = balanced_shortlists(displayed, metadata)
    balance = compare_shortlists(reference, balanced or dict.fromkeys(BALANCE_SEEDS, ([], None)))
    balance = balance.assign(mice_per_group=n_mice)
    if not balanced:
        balance["status"] = "insufficient_mice_per_group"
    balance.to_csv(results / f"{prefix}_balanced_mouse_adaptive_stability.csv", index=False)
    save_figure(
        plot_recovery(balance, chain, "balanced-mouse subtraction stability"),
        f"{chain}_g1_balanced_mouse_stability.png",
    )

    matrix, summary = regional_tables(profiles, genes, seed)
    matrix.to_csv(results / f"{prefix}_adaptive_only_g1_by_seed.csv")
    summary.to_csv(results / f"{prefix}_adaptive_only_g1_summary.csv")
    save_figure(plot_regional_stability(matrix, summary, seed, chain),
                f"{chain}_g1_adaptive_only_g1_stability.png")
