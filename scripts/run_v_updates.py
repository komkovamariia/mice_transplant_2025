#!/usr/bin/env python3
"""Refresh V-segment figures or calibrate Equation 7 without running other analyses."""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import sys
from itertools import combinations
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_analysis import (
    _base_environment, archive_first_approach_figures, execute_notebook, repo_root,
)
from src.spike_in import load_snapshot_tables
from src.v_thresholds import DROPOUT_GROUPS, THRESHOLD_METHOD, fit_mouse_dropout_thresholds


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
        seeds = json.loads(seeds_path.read_text())
        recorded = Path(seeds["baseline_results"])
        if baseline is not None and recorded != baseline.resolve():
            raise ValueError("Existing inputs were prepared from a different baseline. Use the prepare stage to replace them.")
        summary_path = recorded / "run_summary.json"
        if summary_path.is_file():
            summary = json.loads(summary_path.read_text())
            if (seeds["structural_seed"], seeds["depth"]) == (
                summary["structural_downsample_seed"], summary["structural_downsample_umi"],
            ):
                return
            baseline = recorded
    prepare_inputs(root, stratum, baseline, environment)


def _metadata(cache):
    return pd.read_csv(cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})


def _baseline_results(cache):
    return Path(json.loads((cache / "v_input_seeds.json").read_text())["baseline_results"])


def _write_csv(frame, path):
    pending = path.with_suffix(".csv.tmp")
    frame.to_csv(pending, index=False)
    pending.replace(path)


def update_threshold_figures(root, stratum):
    import matplotlib.pyplot as plt
    from src.article_figures import (
        ARTICLE_FONTS, article_dpi, article_name, eligible_segments, plot_dumbbell, plot_edger,
        plot_localization, plot_retention, plot_subsets,
    )
    from src.figures import apply_article_style
    from src.strata import STRATUM_LABELS
    from src.structural_v import subset_v_tables, v_plot_profile

    cache = input_dir(root, stratum)
    metadata = _metadata(cache)
    baseline = _baseline_results(cache)
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
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ARTICLE_FONTS, "axes.grid": False})
    title = "combined CD4+ and CD8+" if stratum == "all_combined" else STRATUM_LABELS[stratum]
    suffix = "" if stratum == "all_combined" else f" | {title}"
    for chain in CHAINS:
        tables = load_snapshot_tables(cache / "structural", chain, metadata)
        audit = pd.read_csv(cache / f"{chain.lower()}_structural_depth_audit.csv", dtype={"sample_id": str})
        cohort = {group: pd.DataFrame(columns=rows["sample_id"]) for group, rows in audit.groupby("group_no")}
        thresholds = estimate(chain, stratum, tables, cohort)
        profile = v_plot_profile(tables)
        eligible = eligible_segments(profile, thresholds).sort_values(["score_s", "v"], ascending=[False, True])
        _write_csv(eligible, results / f"{chain.lower()}_selected_v.csv")
        save(plot_dumbbell(profile, chain, suffix, thresholds=thresholds), chain, "dumbbell")
        bubble, genes = plot_retention(profile, chain, title, thresholds=thresholds)
        save(bubble, chain, "retention")
        save(plot_localization(profile, genes, chain, suffix, thresholds=thresholds), chain, "localization")
        _write_csv(pd.DataFrame({"v": genes}), results / f"{chain.lower()}_bubble_labels.csv")
        ranking_file = baseline / f"{chain.lower()}v_ranking.csv"
        if stratum == "all_combined" or ranking_file.is_file():
            save(plot_edger(pd.read_csv(ranking_file), chain, suffix), chain, "edger")
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
    from src.eq7_calibration import attach_mouse_response, training_evidence, with_structural_fields
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
        build_chain_evidence, attach_mouse_response, training_evidence, with_structural_fields,
        structural_remainder, subset_v_tables, v_region_proportions,
    ):
        digest.update(inspect.getsource(function).encode())
    return cache / "eq7_evidence" / chain / digest.hexdigest()


def build_chain_evidence(directory, tables, raw, metadata, workers):
    from src.differential_abundance import prepare_counts
    from src.eq7_calibration import attach_mouse_response, training_evidence
    from src.structural_v import subset_v_tables

    counts, units = prepare_counts(raw, metadata, pool_within_mouse=True)
    mice = [
        (group, mouse)
        for group in DROPOUT_GROUPS
        for mouse in sorted(metadata.loc[metadata["sample_id"].isin(raw[group].columns), "mouse_id"].unique())
    ]
    ids_by_mouse = metadata.groupby(["group_no", "mouse_id"])["sample_id"].agg(list).to_dict()
    held = [
        subset_v_tables({group: tables[group]}, ids_by_mouse[(group, mouse)])[group]
        for group, mouse in mice
    ]
    records, singles, pairs = [], [], []
    for size in (1, 2):
        for split in combinations(range(len(mice)), size):
            i, *other = split
            j = other[0] if other else -1
            excluded = [mice[index] for index in split]
            affected_groups = {group for group, _ in excluded}
            excluded_ids = {sample for key in excluded for sample in ids_by_mouse[key]}
            retained_ids = metadata.loc[~metadata["sample_id"].isin(excluded_ids), "sample_id"]
            remaining = subset_v_tables({group: tables[group] for group in affected_groups}, retained_ids)
            cache_path = directory / ("split_" + "_".join(map(str, split)) + ".csv")
            if any(held[index].empty for index in split):
                status = "no eligible libraries for an omitted mouse"
            elif any(table.empty for table in remaining.values()):
                status = "omission emptied a structural group"
            elif cache_path.is_file():
                rows = pd.read_csv(cache_path, dtype={"heldout_mouse": str})
                status = "fitted"
            else:
                training = {**tables, **remaining}
                evidence, status = training_evidence(training, excluded, counts, units, n_workers=workers)
                if evidence is not None:
                    rows = pd.concat([
                        attach_mouse_response(evidence, training, held[index], *mice[index]).assign(
                            heldout_fold=index, outer_fold=j if index == i else i,
                        )
                        for index in split
                    ], ignore_index=True)
                    _write_csv(rows, cache_path)
            if status == "fitted":
                (pairs if other else singles).append(rows)
            records.append({
                "fold": i, "paired_fold": j, "group": mice[i][0], "mouse": mice[i][1],
                "paired_group": mice[j][0] if other else "", "paired_mouse": mice[j][1] if other else "",
                "status": status,
            })
            print(f"Omitted {', '.join(mouse for _, mouse in excluded)}: {status}", flush=True)

    empty = pd.DataFrame(columns=["heldout_fold", "outer_fold"])
    _write_csv(pd.concat(singles, ignore_index=True) if singles else empty, directory / "fold_evidence.csv")
    _write_csv(pd.concat(pairs, ignore_index=True) if pairs else empty, directory / "inner_fold_evidence.csv")
    _write_csv(pd.DataFrame(records), directory / "fold_status.csv")


def build_weight_evidence(root, cache):
    from src.runtime import configure_runtime
    from src.spike_in import snapshot_tables
    from src.structural_sampling import downsample_group_tables

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
        status = directory / "fold_status.csv"
        if status.is_file():
            print(f"Reusing {chain} evidence in {directory}", flush=True)
        else:
            print(f"Preparing {chain} single and paired omissions", flush=True)
            build_chain_evidence(directory, tables, raw, metadata, workers)


def fit_outer_weights(directory, output, prefix, top, ridge):
    from src.eq7_calibration import CalibrationUnavailable, calibrate_weights, score_candidates

    outer_rows = pd.read_csv(directory / "fold_evidence.csv", dtype={"heldout_mouse": str})
    inner_rows = pd.read_csv(directory / "inner_fold_evidence.csv", dtype={"heldout_mouse": str})
    fold_status = pd.read_csv(directory / "fold_status.csv", dtype={"mouse": str})
    outcomes = ["heldout_n_v", "heldout_n_recovered", "response"]
    records, coefficients, scores = [], [], []
    for fold in fold_status.loc[fold_status["paired_fold"].eq(-1)].itertuples(index=False):
        key = {"outer_fold": fold.fold, "outer_group": fold.group, "outer_mouse": fold.mouse}
        status = fold.status
        n_candidates = n_observed = 0
        if status == "fitted":
            try:
                weights, _ = calibrate_weights(inner_rows.loc[inner_rows["outer_fold"].eq(fold.fold)], top, ridge)
            except CalibrationUnavailable as exc:
                status = str(exc)
            else:
                evidence = outer_rows.loc[outer_rows["heldout_fold"].eq(fold.fold)]
                candidates = score_candidates(evidence.drop(columns=outcomes), weights)
                candidates = candidates.join(evidence.set_index("v_gene")[outcomes], on="v_gene", validate="1:1")
                candidates = candidates.assign(**key)
                n_candidates = len(candidates)
                n_observed = int(candidates["response"].notna().sum())
                if n_observed:
                    status = "scored"
                elif n_candidates:
                    status = "No selected V segments observed in the outer mouse."
                else:
                    status = "No candidates passed the outer training cutoffs."
                coefficients.append(weights.assign(**key))
                scores.append(candidates)
        records.append({**key, "status": status, "n_candidates": n_candidates, "n_observed": n_observed})

    status = pd.DataFrame(records)
    _write_csv(status, output / f"{prefix}_outer_status.csv")
    empty_weights = pd.DataFrame(columns=["outer_fold", "outer_group", "outer_mouse", "term", "weight", "intercept"])
    empty_scores = pd.DataFrame(columns=[
        "outer_fold", "outer_group", "outer_mouse", "v_gene", "E", "response", "predicted_response",
    ])
    predictions = pd.concat(scores, ignore_index=True) if scores else empty_scores
    _write_csv(pd.concat(coefficients, ignore_index=True) if coefficients else empty_weights,
               output / f"{prefix}_outer_coefficients.csv")
    _write_csv(predictions, output / f"{prefix}_outer_scores.csv")
    return outer_rows, predictions, status


def fit_weights(root, cache, top, ridge):
    import matplotlib.pyplot as plt
    from src.article_figures import ARTICLE_FONTS
    from src.eq7_calibration import (
        CalibrationUnavailable, calibrate_weights, score_candidates, with_structural_fields,
    )
    from src.eq7_figures import plot_eq7_inputs, plot_eq7_regression, plot_eq7_scores
    from src.figures import apply_article_style
    from src.spike_in import EDGER_SUCCESS

    metadata = _metadata(cache)
    baseline = _baseline_results(cache)
    edger_status = json.loads((baseline / "run_summary.json").read_text())["edger_status_by_chain"]
    output = cache / "eq7_weights"
    output.mkdir(parents=True, exist_ok=True)
    summaries = []
    apply_article_style()
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ARTICLE_FONTS, "axes.grid": False})
    for chain in CHAINS:
        calibration = calibration_dir(root, cache, chain)
        if not (calibration / "count_snapshots" / f"{chain}_manifest.json").is_file():
            raise FileNotFoundError(f"Run the evidence stage for {chain} before fitting weights.")
        directory = evidence_dir(root, cache, chain)
        if not (directory / "fold_status.csv").is_file():
            raise FileNotFoundError(f"Run the evidence stage for {chain} before fitting weights.")
        prefix = chain.lower()
        rows, outer, outer_status = fit_outer_weights(directory, output, prefix, top, ridge)
        try:
            weights, selected = calibrate_weights(rows, top, ridge)
        except CalibrationUnavailable as exc:
            status = str(exc)
            weights = pd.DataFrame(columns=["term", "weight", "varied_in_training", "intercept"])
            selected = rows.iloc[:0]
        else:
            status = "fitted"
        _write_csv(weights.assign(top=top, ridge=ridge, threshold_method=THRESHOLD_METHOD),
                   output / f"{prefix}_coefficients.csv")
        _write_csv(selected, output / f"{prefix}_calibration_rows.csv")
        tables = load_snapshot_tables(cache / "structural", chain, metadata)
        ranking = with_structural_fields(pd.read_csv(baseline / f"{prefix}v_ranking.csv"), tables)
        score_status = status
        score_weights = weights
        if edger_status.get(chain) != EDGER_SUCCESS:
            score_status = "Full-cohort edgeR results unavailable."
            score_weights = weights.iloc[:0]
        scored = score_candidates(ranking, score_weights)
        _write_csv(scored, output / f"{prefix}_scored_top_v.csv")
        summaries.append({
            "chain": chain, "threshold_method": THRESHOLD_METHOD,
            "weight_status": status, "full_score_status": score_status,
            "n_selected": len(scored), "n_with_E": int(scored["E"].notna().sum()),
            "n_outer_folds_scored": int(outer_status["status"].eq("scored").sum()),
        })
        figure_dir = root / "figures" / FIGURE_APPROACH / "all_combined" / chain
        figure_dir.mkdir(parents=True, exist_ok=True)
        for name, figure in (
            ("inputs", plot_eq7_inputs(selected, weights, chain, status)),
            ("regression", plot_eq7_regression(selected, outer, outer_status, chain, status)),
            ("selected_v_scores", plot_eq7_scores(scored, chain, score_status)),
        ):
            figure.savefig(figure_dir / f"{chain}_g1_eq7_{name}.png", dpi=300, facecolor="white")
            plt.close(figure)
        print(f"{chain} full-cohort calibration: {status}. Reports written to {output}")
    _write_csv(pd.DataFrame(summaries), output / "calibration_summary.csv")
    print(f"Figure archive written to {archive_first_approach_figures(output_name=FIGURE_APPROACH)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("all", "prepare", "thresholds", "evidence", "weights"))
    parser.add_argument("--baseline-results", type=Path, help="Previous all_combined results directory.")
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--metadata-csv", type=Path)
    parser.add_argument("--clonoset-index", type=Path)
    parser.add_argument("--working-dir", type=Path)
    parser.add_argument("--top", type=int, help="Optional cap on eligible V segments per calibration fold.")
    parser.add_argument("--ridge", type=float, default=0.1)
    args = parser.parse_args()
    if (args.top is not None and args.top < 1) or not 0 < args.ridge < float("inf"):
        parser.error("--top must be positive and --ridge must be finite and positive")
    root, stratum = repo_root(), "all_combined"
    cache = input_dir(root, stratum)
    if args.stage == "prepare":
        prepare_inputs(root, stratum, args.baseline_results, _base_environment(args))
    elif args.stage in {"all", "thresholds"}:
        prepare_if_missing(root, stratum, args.baseline_results, _base_environment(args))
        update_threshold_figures(root, stratum)
        if args.stage == "all":
            build_weight_evidence(root, cache)
            fit_weights(root, cache, args.top, args.ridge)
    elif args.stage == "evidence":
        build_weight_evidence(root, cache)
    else:
        fit_weights(root, cache, args.top, args.ridge)


if __name__ == "__main__":
    main()
