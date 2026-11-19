#!/usr/bin/env python3
"""Refresh V-segment figures or calibrate Equation 7 without running other analyses."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_analysis import (
    _base_environment, archive_first_approach_figures, execute_notebook, repo_root,
)
from src.spike_in import load_snapshot_tables
from src.v_thresholds import DROPOUT_GROUPS, fit_mouse_dropout_thresholds, mouse_omissions


INPUT_CELLS = ("venn-parameters", "venn-inputs", "venn-historical-code-06", "venn-historical-code-15")
CHAINS = ("TRA", "TRB")
FIGURE_APPROACH = "01_set_count_adaptive"


def input_dir(root, stratum):
    return root / "results" / "01_v_update_inputs" / stratum


def prepare_inputs(root, stratum, baseline, environment):
    if baseline is None:
        baseline = root / "results" / "01_set_count" / stratum
        if not (baseline / "run_summary.json").is_file():
            baseline = root / "results" / FIGURE_APPROACH / stratum
    baseline = baseline.resolve()
    summary_path = baseline / "run_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Provide the previous run's stratum directory with run_summary.json: {baseline}")
    summary = json.loads(summary_path.read_text())
    if summary["stratum"] != stratum or summary["mode"] != "standard":
        raise ValueError("The baseline must be a standard run for the requested stratum.")
    execute_notebook(
        root / "venn_original.ipynb", f"v_inputs/{stratum}", resume=False, cell_ids=INPUT_CELLS,
        environment={
            **environment, "MICE_TCR_MODE": "v-inputs", "MICE_TCR_STRATUM": stratum,
            "MICE_TCR_V_BASELINE_RESULTS": str(baseline),
        },
    )


def prepare_if_missing(root, stratum, baseline, environment):
    cache = input_dir(root, stratum)
    seeds_path = cache / "v_input_seeds.json"
    complete = seeds_path.is_file() and (cache / "mouse_metadata.csv").is_file() and all(
        (cache / source / "count_snapshots" / f"{chain}_manifest.json").is_file()
        for source in ("raw", "structural") for chain in CHAINS
    )
    if complete:
        recorded = Path(json.loads(seeds_path.read_text())["baseline_results"])
        if baseline is not None and recorded != baseline.resolve():
            raise ValueError("Existing inputs were prepared from a different baseline. Use the prepare stage to replace them.")
        return
    prepare_inputs(root, stratum, baseline, environment)


def _metadata(cache):
    return pd.read_csv(cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})


def _write_csv(frame, path):
    pending = path.with_suffix(".csv.tmp")
    frame.to_csv(pending, index=False)
    pending.replace(path)


def update_threshold_figures(root, stratum):
    import matplotlib.pyplot as plt
    from src.article_figures import (
        ARTICLE_FONTS, article_dpi, article_name, plot_dumbbell, plot_localization, plot_retention, plot_subsets,
    )
    from src.figures import apply_article_style
    from src.strata import STRATUM_LABELS
    from src.structural_v import subset_v_tables, v_plot_profile

    cache = input_dir(root, stratum)
    metadata = _metadata(cache)
    results = root / "results" / FIGURE_APPROACH / stratum
    figures = root / "figures" / FIGURE_APPROACH / stratum
    results.mkdir(parents=True, exist_ok=True)
    records = []

    def estimate(chain, panel, tables, cohort=None):
        cutoffs, rows = fit_mouse_dropout_thresholds(tables, metadata, cohort_tables=cohort)
        records.extend({"chain": chain, "panel": panel, **row} for row in rows)
        return cutoffs

    def save(figure, chain, panel):
        name = article_name(chain, panel)
        path = figures / chain / name
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=article_dpi(name), facecolor="white")
        plt.close(figure)

    apply_article_style()
    plt.rcParams.update({"font.sans-serif": ARTICLE_FONTS, "axes.grid": False})
    title = "combined CD4+ and CD8+" if stratum == "all_combined" else STRATUM_LABELS[stratum]
    suffix = "" if stratum == "all_combined" else f" | {title}"
    for chain in CHAINS:
        tables = load_snapshot_tables(cache / "structural", chain, metadata)
        audit = pd.read_csv(cache / f"{chain.lower()}_structural_depth_audit.csv", dtype={"sample_id": str})
        cohort = {group: pd.DataFrame(columns=rows["sample_id"]) for group, rows in audit.groupby("group_no")}
        thresholds = estimate(chain, stratum, tables, cohort)
        profile = v_plot_profile(tables)
        save(plot_dumbbell(profile, chain, suffix, thresholds=thresholds), chain, "dumbbell")
        bubble, genes = plot_retention(profile, chain, title, thresholds=thresholds)
        save(bubble, chain, "retention")
        save(plot_localization(profile, genes, chain, suffix, thresholds=thresholds), chain, "localization")
        _write_csv(pd.DataFrame({"v": genes}), results / f"{chain.lower()}_bubble_labels.csv")
        if stratum == "all_combined":
            profiles, cutoffs = {}, {}
            for panel in ("CD8+", "CD4+"):
                ids = metadata.loc[metadata["_cell_subset"].eq(panel[:3].lower()), "sample_id"]
                subset = subset_v_tables(tables, ids)
                profiles[panel] = v_plot_profile(subset)
                cutoffs[panel] = estimate(chain, panel, subset, cohort)
            save(plot_subsets(profiles, chain, thresholds_by_subset=cutoffs), chain, "cd4_cd8")
    _write_csv(pd.DataFrame(records), results / "v_thresholds.csv")
    print(f"Adaptive thresholds written to {results / 'v_thresholds.csv'}")
    print(f"Figure archive written to {archive_first_approach_figures(output_name=FIGURE_APPROACH)}")


def calibration_dir(root, cache, chain):
    digest = hashlib.sha256()
    digest.update((cache / "raw" / "count_snapshots" / f"{chain}_manifest.json").read_bytes())
    seeds = json.loads((cache / "v_input_seeds.json").read_text())
    digest.update(repr((seeds["depth"], seeds["calibration_seed"])).encode())
    digest.update((root / "src" / "structural_sampling.py").read_bytes())
    return cache / "calibration" / chain / digest.hexdigest()


def evidence_dir(root, cache, chain):
    from src.eq7_calibration import fold_evidence, with_structural_fields
    from src.structural_v import structural_remainder, subset_v_tables, v_region_proportions

    files = [
        cache / "raw" / "count_snapshots" / f"{chain}_manifest.json",
        calibration_dir(root, cache, chain) / "count_snapshots" / f"{chain}_manifest.json",
        root / "src" / "differential_abundance.py",
        root / "src" / "count_statistics.py",
    ]
    digest = hashlib.sha256()
    for path in files:
        digest.update(path.read_bytes())
    digest.update(repr(DROPOUT_GROUPS).encode())
    for function in (
        mouse_omissions, fold_evidence, with_structural_fields,
        structural_remainder, subset_v_tables, v_region_proportions,
    ):
        digest.update(inspect.getsource(function).encode())
    return cache / "eq7_evidence" / chain / digest.hexdigest()


def build_weight_evidence(root, cache):
    from src.differential_abundance import prepare_counts
    from src.eq7_calibration import fold_evidence
    from src.runtime import configure_runtime
    from src.spike_in import snapshot_tables
    from src.structural_sampling import downsample_group_tables
    from src.structural_v import subset_v_tables

    metadata = _metadata(cache)
    seeds = json.loads((cache / "v_input_seeds.json").read_text())
    workers = configure_runtime()
    for chain in CHAINS:
        calibration = calibration_dir(root, cache, chain)
        raw = load_snapshot_tables(cache / "raw", chain, metadata)
        if (calibration / "count_snapshots" / f"{chain}_manifest.json").is_file():
            tables = load_snapshot_tables(calibration, chain, metadata)
        else:
            tables, _ = downsample_group_tables(raw, depth=seeds["depth"], seed=seeds["calibration_seed"])
            snapshot_tables(tables, calibration, chain, metadata)
        directory = evidence_dir(root, cache, chain)
        directory.mkdir(parents=True, exist_ok=True)
        if (directory / "fold_evidence.csv").is_file():
            print(f"Reusing fold evidence from {directory / 'fold_evidence.csv'}", flush=True)
            continue
        counts, units = prepare_counts(raw, metadata, pool_within_mouse=True)
        records, frames = [], []
        omissions = mouse_omissions(tables, metadata, cohort_tables=raw)
        for i, (group, mouse, training, status) in enumerate(omissions, 1):
            output = directory / f"fold_{i:03d}.csv"
            if training is not None:
                if output.is_file():
                    frames.append(pd.read_csv(output))
                    status = "cached"
                else:
                    ids = metadata.loc[metadata["mouse_id"].eq(mouse) & metadata["group_no"].eq(group), "sample_id"]
                    held = subset_v_tables({group: tables[group]}, ids)[group]
                    rows, status = fold_evidence(training, held, group, mouse, counts, units, n_workers=workers)
                    if rows is not None:
                        _write_csv(rows, output)
                        frames.append(rows)
            records.append({"heldout_group": group, "heldout_mouse": mouse, "status": status})
            print(f"{chain} omitted {mouse}: {status}", flush=True)
        _write_csv(pd.DataFrame(records), directory / "fold_status.csv")
        if not frames:
            raise ValueError(f"No usable Equation 7 folds for {chain}; see {directory / 'fold_status.csv'}")
        _write_csv(pd.concat(frames, ignore_index=True), directory / "fold_evidence.csv")
        print(f"Fold evidence written to {directory / 'fold_evidence.csv'}")


def fit_weights(root, cache, top, ridge):
    from src.eq7_calibration import calibrate_weights, score_candidates, with_structural_fields
    from src.spike_in import EDGER_SUCCESS

    metadata = _metadata(cache)
    baseline = Path(json.loads((cache / "v_input_seeds.json").read_text())["baseline_results"])
    edger_status = json.loads((baseline / "run_summary.json").read_text())["edger_status_by_chain"]
    output = cache / "eq7_weights"
    output.mkdir(parents=True, exist_ok=True)
    for chain in CHAINS:
        if edger_status.get(chain) != EDGER_SUCCESS:
            raise ValueError(f"The saved full-cohort {chain} edgeR model was not fitted successfully.")
        calibration = calibration_dir(root, cache, chain)
        if not (calibration / "count_snapshots" / f"{chain}_manifest.json").is_file():
            raise FileNotFoundError(f"Run the evidence stage for {chain} before fitting weights.")
        source = evidence_dir(root, cache, chain) / "fold_evidence.csv"
        if not source.is_file():
            raise FileNotFoundError(f"Run the evidence stage for {chain} before fitting weights.")
        weights, selected = calibrate_weights(pd.read_csv(source), top, ridge)
        tables = load_snapshot_tables(calibration, chain, metadata)
        prefix = f"{chain.lower()}v"
        ranking = with_structural_fields(pd.read_csv(baseline / f"{prefix}_ranking.csv"), tables)
        _write_csv(weights.assign(top=top, ridge=ridge), output / f"{chain.lower()}_coefficients.csv")
        _write_csv(selected, output / f"{chain.lower()}_calibration_rows.csv")
        _write_csv(score_candidates(ranking, weights, top), output / f"{chain.lower()}_scored_top_v.csv")
        print(f"{chain} weights and per-segment contributions written to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "thresholds", "evidence", "weights"))
    parser.add_argument("--baseline-results", type=Path, help="Previous all_combined results directory.")
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--metadata-csv", type=Path)
    parser.add_argument("--clonoset-index", type=Path)
    parser.add_argument("--working-dir", type=Path)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--ridge", type=float, default=0.1)
    args = parser.parse_args()
    if args.top < 1 or not 0 < args.ridge < float("inf"):
        parser.error("--top must be positive and --ridge must be finite and positive")
    root, stratum = repo_root(), "all_combined"
    cache = input_dir(root, stratum)
    if args.stage == "prepare":
        prepare_inputs(root, stratum, args.baseline_results, _base_environment(args))
    elif args.stage == "thresholds":
        prepare_if_missing(root, stratum, args.baseline_results, _base_environment(args))
        update_threshold_figures(root, stratum)
    elif args.stage == "evidence":
        build_weight_evidence(root, cache)
    else:
        fit_weights(root, cache, args.top, args.ridge)


if __name__ == "__main__":
    main()
