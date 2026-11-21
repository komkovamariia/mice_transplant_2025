"""Rebuild manuscript figures from completed analysis outputs."""

import json

import matplotlib.pyplot as plt
import pandas as pd

from src.article_figures import article_name, plot_fisher, plot_group_overlap
from src.differential_abundance import count_result_path
from src.eq7_figures import plot_eq7_inputs, plot_eq7_regression, plot_eq7_scores, plot_ridge_selection
from src.publication_figures import save_publication_figure
from src.repertoire_depth import plot_depth_check, plot_mds
from src.spike_in import load_snapshot_tables
from src.v_stability import plot_recovery, plot_regional_stability


CHAINS = ("TRA", "TRB")
STABILITY_TABLES = (
    "rarefaction_adaptive_stability", "balanced_mouse_adaptive_stability",
    "adaptive_only_g1_by_seed", "adaptive_only_g1_summary",
)


def save(figure, output):
    try:
        save_publication_figure(figure, output)
    finally:
        plt.close(figure)


def redraw_stability(source, output, chain, seed, depth):
    prefix = chain.lower()
    files = [source / f"{prefix}_{name}.csv" for name in STABILITY_TABLES]
    if not all(path.is_file() for path in files):
        print(f"{chain}: adaptive stability CSVs missing in {source}; run the stability stage.")
        return
    rarefaction, balance = [pd.read_csv(path) for path in files[:2]]
    matrix, summary = [pd.read_csv(path, index_col=0) for path in files[2:]]
    output = output / chain / "surplus"
    for name, figure in (
        ("rarefaction_stability", plot_recovery(rarefaction, chain, f"{depth:,}-UMI rarefaction stability",
                                              selected_seed=seed)),
        ("balanced_mouse_stability", plot_recovery(balance, chain, "Balanced-mouse subtraction stability")),
        ("adaptive_only_g1_stability", plot_regional_stability(matrix, summary, seed, chain)),
    ):
        save(figure, output / f"{chain}_g1_{name}.png")


def redraw_threshold_audit(source, output, chain):
    path = source / "v_thresholds.csv"
    if not path.is_file():
        return
    results = pd.read_csv(path).query("chain == @chain and panel == 'all_combined'")
    full = results.loc[results["level"].eq("final")]
    omitted = results.loc[results["level"].eq("mouse")]
    if len(full) != 1 or omitted.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.5), layout="constrained")
    for i, (ax, column, label) in enumerate(zip(
        axes, ("initial_share_pct", "retention_pct"), ("Initial share cutoff (%)", "Retention cutoff (%)"),
        strict=True,
    )):
        ax.set_title(f"{'AB'[i]}. Cutoffs")
        ax.axvline(full[column].iat[0], color="black", linestyle="--", lw=0.8)
        ax.scatter(omitted[column], range(len(omitted)), color="#805282", s=27)
        ax.set(xlabel=label, yticks=range(len(omitted)), yticklabels=omitted["omitted_mouse"])
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.22)
    save(fig, output / chain / "surplus" / f"{chain}_g1_threshold_audit.png")


def redraw_eq7(source, output, chain):
    summary_path = source / "calibration_summary.csv"
    if not summary_path.is_file():
        print(f"{chain}: no completed Equation 7 summary in {source}. Check after the weights stage.")
        return
    summary = pd.read_csv(summary_path).set_index("chain").loc[chain]
    prefix = chain.lower()
    names = ("calibration_rows", "coefficients", "outer_scores", "outer_status", "ridge_grid", "scored_top_v")
    files = [source / f"{prefix}_{name}.csv" for name in names]
    missing = [path.name for path in files if not path.is_file()]
    if missing:
        print(f"{chain}: Equation 7 figure inputs missing: {', '.join(missing)}.")
        return
    rows, weights, outer, outer_status, curve, scored = [
        pd.read_csv(path, dtype={"heldout_mouse": str, "outer_mouse": str})
        for path in files
    ]
    for name, figure in (
        ("inputs", plot_eq7_inputs(rows, weights, chain, summary["weight_status"])),
        ("regression", plot_eq7_regression(rows, outer, outer_status, chain, summary["weight_status"])),
        ("penalty", plot_ridge_selection(curve, chain)),
        ("selected_v_scores", plot_eq7_scores(scored, chain, summary["full_score_status"])),
    ):
        save(figure, output / chain / f"{chain}_g1_eq7_{name}.png")


def redraw_baseline(chain, baseline, cache, output, depth, metadata):
    prefix = chain.lower()
    path = baseline / f"{prefix}_depth_check.csv"
    if path.is_file():
        save(plot_depth_check(pd.read_csv(path), chain, depth), output / article_name(chain, "depth"))
    path = baseline / f"{prefix}_mds_coordinates.csv"
    if path.is_file() and metadata is not None:
        coordinates = pd.read_csv(path, index_col=0).rename_axis("sample_id")
        labels = metadata.set_index("sample_id")
        save(plot_mds(coordinates.join(labels, how="left", validate="1:1"), chain, depth),
             output / article_name(chain, "mds"))
    path = count_result_path(baseline, chain, "fisher_edger_feature_comparison.csv")
    if path.is_file():
        save(plot_fisher(pd.read_csv(path), chain), output / article_name(chain, "fisher"))
    snapshot = cache / "structural"
    if metadata is not None:
        tables = load_snapshot_tables(snapshot, chain, metadata)
        save(plot_group_overlap(tables, chain), output / article_name(chain, "overlap"))


def redraw(root, stratum, args):
    from scripts.run_v_updates import input_dir, update_threshold_figures

    cache = input_dir(root, stratum)
    output = args.figure_output or root / "figures" / "01_set_count_mdpi" / stratum
    baseline = args.baseline_results or root / "results" / "01_set_count" / stratum
    provenance = cache / "v_input_seeds.json"
    snapshot_ready = provenance.is_file() and (cache / "mouse_metadata.csv").is_file() and all(
        (cache / "structural" / "count_snapshots" / f"{chain}_manifest.json").is_file() for chain in CHAINS
    )
    seeds = json.loads(provenance.read_text()) if provenance.is_file() else {}
    if args.expect_run is not None:
        actual = [seeds.get(key) for key in ("structural_seed", "calibration_seed", "depth")]
        if actual != args.expect_run:
            raise ValueError(f"Expected cached settings {args.expect_run}; found {actual} in {provenance}.")
    if snapshot_ready:
        recorded_baseline = root / seeds["baseline_results"]
        if args.baseline_results is None:
            baseline = recorded_baseline
        elif recorded_baseline.is_dir() and baseline.resolve() != recorded_baseline.resolve():
            raise ValueError("The supplied baseline differs from the saved structural snapshots.")
        metadata = pd.read_csv(cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})
    else:
        metadata = None
    summary = baseline / "run_summary.json"
    run = json.loads(summary.read_text()) if summary.is_file() else {}
    seed = seeds["structural_seed"] if snapshot_ready else run.get("structural_downsample_seed")
    depth = seeds["depth"] if snapshot_ready else run.get("structural_downsample_umi", 15000)
    if snapshot_ready and run and (
        run.get("structural_downsample_seed") != seed or run.get("structural_downsample_umi") != depth
    ):
        raise ValueError("The baseline metadata disagree with the cached structural seed or depth.")
    if seeds:
        print(f"Structural seed {seed}, calibration seed {seeds['calibration_seed']}, depth {depth}")
    components = set(args.components or ("structural", "stability", "eq7"))

    if "structural" in components:
        if snapshot_ready:
            update_threshold_figures(root, stratum, figure_output=output, baseline_results=baseline)
        else:
            print("Structural snapshots are unavailable. Run the prepare stage to restore adaptive plots.")
        for chain in CHAINS:
            redraw_baseline(chain, baseline, cache, output / chain, depth, metadata)
    if "stability" in components:
        adaptive = root / "results" / "01_set_count_adaptive" / stratum
        for chain in CHAINS:
            source = args.stability_results
            if source is None:
                has_cached = all((adaptive / f"{chain.lower()}_{name}.csv").is_file() for name in STABILITY_TABLES)
                source = adaptive if has_cached else baseline
            plot_seed, plot_depth = seed, depth
            if (source / "run_summary.json").is_file():
                settings = json.loads((source / "run_summary.json").read_text())
                plot_seed = settings.get("structural_downsample_seed", seed)
                plot_depth = settings.get("structural_downsample_umi", depth)
            redraw_stability(source, output, chain, plot_seed, plot_depth)
            redraw_threshold_audit(adaptive, output, chain)
    if "eq7" in components:
        for chain in CHAINS:
            redraw_eq7(cache / "eq7_weights", output, chain)
    if seeds:
        output.mkdir(parents=True, exist_ok=True)
        fields = ("structural_seed", "calibration_seed", "depth")
        (output / "figure_run.json").write_text(
            json.dumps({key: seeds[key] for key in fields}, indent=2) + "\n",
        )
    print(f"Figure directory: {output}")
