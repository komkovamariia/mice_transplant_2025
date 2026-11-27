#!/usr/bin/env python3
"""Compare aaV-positive fractions with existing Eq7 evidence on cached mouse splits.

No rarefaction, mouse splitting, structural subtraction or calibration preparation
is performed. The 75 recorded triplets per chain are read from the 1249 cache.
Only new raw-count aaV-level edgeR/Fisher tests are performed and cached per fold.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_continuous import (
    JOINT_CONTRAST, JOINT_VALUE, ORDERS, check_input_recovery, load_inputs,
    prediction_summary, ranking_rows, ranking_summary, split_mice,
)
from src.eq7_aav_fraction import FRACTION, FRACTION_MODEL, fraction_evidence, prepare_aav_counts
from src.eq7_calibration import FEATURES
from src.eq7_continuous import BINARY_MODEL, CONTINUOUS_MODEL, JointRecoveryPanel
from src.eq7_weighted import GROUPS, contrasts, fit_ridge, split_shortlist, training_contrasts
from src.spike_in import load_snapshot_tables


PENALTIES = (1.0, 10.0)
ABLATIONS = {
    "aav_fraction_without_edgeR": (FRACTION[1], FRACTION[2]),
    "aav_fraction_without_Fisher": (FRACTION[0], FRACTION[2]),
    "aav_fraction_without_structure": FRACTION[:2],
}
SINGLE_PREDICTOR_MODELS = {
    "structure_only": (FRACTION[2],),
    "edgeR_only_binary": (FEATURES[0],),
    "edgeR_only_fraction": (FRACTION[0],),
}
MODELS = {FRACTION_MODEL: FRACTION, **ABLATIONS, **SINGLE_PREDICTOR_MODELS}
RANKERS = {
    "structural_score": "n_g1_exclusive_aav",
    "aav_edgeR_positive_fraction": FRACTION[0],
    "aav_Fisher_positive_fraction": FRACTION[1],
}


def cache_signature(cache, chain, root):
    digest = hashlib.sha256()
    digest.update((cache / "raw" / "count_snapshots" / f"{chain}_manifest.json").read_bytes())
    for relative in (
        "src/eq7_aav_fraction.py",
        "src/differential_abundance.py",
        "src/count_statistics.py",
    ):
        digest.update((root / relative).read_bytes())
    from rpy2.robjects import r
    digest.update(str(r('as.character(packageVersion("edgeR"))')[0]).encode())
    digest.update(b"g1-vs-g5g6;BH<0.05;logFC>0;denominator=tested;aaV")
    return digest.hexdigest()


def cached_fractions(root, cache, output, chain, raw, metadata, folds, workers):
    """Resume from per-triplet CSVs if present; never rebuild original folds."""
    counts, units = prepare_aav_counts(raw, metadata)
    directory = output / "effect_cache" / chain / cache_signature(cache, chain, root)
    directory.mkdir(parents=True, exist_ok=True)
    total = len(folds)
    for index, (split, rows) in enumerate(folds.items(), 1):
        mice = split_mice(rows)
        label = hashlib.sha256(repr(mice).encode()).hexdigest()[:16]
        path = directory / f"{split}_{label}.csv"
        if path.is_file():
            effects = pd.read_csv(path)
            source = "cached"
        else:
            effects = fraction_evidence(counts, units, mice, workers=workers)
            pending = path.with_suffix(".tmp")
            effects.to_csv(pending, index=False)
            pending.replace(path)
            source = "computed"
        if not set(FRACTION[:2]).issubset(effects.columns):
            raise ValueError(f"Invalid aaV feature cache: {path}")
        if index % 5 == 0 or index == 1:
            print(f"{chain}: aaV evidence {index}/{total} ({source})", flush=True)
        yield split, effects


def attach_fractions(shortlist, effects):
    result = shortlist.merge(effects, on="v_gene", how="left", validate="one_to_one")
    if result[list(FRACTION[:2])].isna().any().any():
        missing = result.loc[result[FRACTION[0]].isna(), "v_gene"].tolist()
        raise ValueError(f"Missing aaV evidence for shortlisted segments: {missing}")
    return result


def evaluate(chain, panel, folds, effects):
    lookup = {key: index for index, key in enumerate(panel.mice)}
    inner_measures, ranks, predictions, for_final_fit = {}, [], [], []
    for split, fold in folds.items():
        mice = split_mice(fold)
        indices = tuple(lookup[key] for key in mice)
        observed = panel.measure(panel.reference(indices), indices)
        check_input_recovery(fold, panel, indices, observed)
        eligible = attach_fractions(split_shortlist(fold), effects[split])
        if eligible.empty:
            raise ValueError(f"No eligible V segments in {chain} {split}")
        genes = panel.genes.get_indexer(eligible["v_gene"])
        if np.any(genes < 0):
            raise ValueError(f"Unknown shortlisted segment in {chain} {split}")
        truth = contrasts(observed[:, :, genes])[1]
        training = training_contrasts(
            panel, indices, inner_measures, recovery_column=JOINT_CONTRAST,
        )
        training = eligible.merge(
            training[["v_gene", JOINT_CONTRAST]], on="v_gene", validate="one_to_one",
        )
        order = {name: eligible[column].to_numpy(dtype=float) for name, column in RANKERS.items()}
        for model, features in MODELS.items():
            penalties = PENALTIES if model == FRACTION_MODEL or model in SINGLE_PREDICTOR_MODELS else (1.0,)
            for penalty in penalties:
                fit = fit_ridge(training, features, JOINT_CONTRAST, penalty)
                pred = fit.predict(eligible)
                predictions.append(pd.DataFrame({
                    "chain": chain, "split": split, "v_gene": eligible["v_gene"],
                    "g1_mouse": mice[0][1], "g5_mouse": mice[1][1], "g6_mouse": mice[2][1],
                    "model": model, "outcome": JOINT_CONTRAST, "penalty": penalty,
                    "observed": truth, "prediction": pred, "constant_baseline": fit.baseline,
                }))
                if penalty == 1:
                    order[model] = pred
        for row in ranking_rows(
            eligible, observed[:, :, genes], mice, ORDERS[chain], order,
            recovery_column=JOINT_VALUE,
        ):
            ranks.append({"chain": chain, "split": split, **row})
        for_final_fit.append(eligible.assign(**{JOINT_CONTRAST: truth}, fit_weight=1 / len(eligible)))
    final_data = pd.concat(for_final_fit, ignore_index=True)
    coefficients = []
    for penalty in PENALTIES:
        fit = fit_ridge(
            final_data, FRACTION, JOINT_CONTRAST, penalty,
            weight="fit_weight",
        )
        for feature, coef in fit.coefficients.items():
            coefficients.append({
                "chain": chain, "model": FRACTION_MODEL, "penalty": penalty,
                "term": feature, "weight": coef, "ridge_intercept": fit.intercept,
                "constant_baseline": fit.baseline,
            })
    return pd.DataFrame(ranks), pd.concat(predictions, ignore_index=True), pd.DataFrame(coefficients)


def check_reference(new_predictions, old_predictions):
    """Verify old and new models score identical held-out target values."""
    selected = new_predictions.loc[
        new_predictions["model"].eq(FRACTION_MODEL),
        ["chain", "split", "v_gene", "penalty", "observed"],
    ]
    reference = old_predictions.loc[
        old_predictions["model"].isin((BINARY_MODEL, CONTINUOUS_MODEL)),
        ["chain", "split", "v_gene", "penalty", "observed"],
    ]
    matched = reference.merge(
        selected, on=["chain", "split", "v_gene", "penalty"],
        how="left", indicator=True, validate="many_to_one", suffixes=("_old", "_new"),
    )
    if not matched["_merge"].eq("both").all():
        raise ValueError("Reference predictions use different mouse folds or candidate V segments.")
    np.testing.assert_allclose(
        matched["observed_old"], matched["observed_new"], rtol=0, atol=1e-9,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True, help="Existing eq7_benchmark_1249/_inputs only.")
    parser.add_argument("--output", type=Path, default=Path("results/eq7_aav_fraction_1249"))
    parser.add_argument("--reference-results", type=Path, help="Existing eq7_continuous_1249 results to compare unchanged baseline fits.")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers for new per-aaV Fisher tests.")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if not args.cache.joinpath("v_input_seeds.json").is_file():
        parser.error("The prepared seed-1249 _inputs cache is required; no preparation is allowed.")
    provenance = json.loads((args.cache / "v_input_seeds.json").read_text())
    if tuple(provenance[k] for k in ("structural_seed", "calibration_seed", "depth")) != (1249, 1031, 15000):
        parser.error("Expected the existing structural=1249, calibration=1031, depth=15000 cache.")
    root = Path(__file__).resolve().parents[1]
    metadata = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})
    args.output.mkdir(parents=True, exist_ok=True)
    all_ranks, all_predictions, all_coefficients = [], [], []
    for chain in ("TRA", "TRB"):
        panel, folds = load_inputs(root, args.cache, chain, metadata, panel_type=JointRecoveryPanel)
        if len(folds) != 75:
            raise ValueError(f"Expected 75 existing triplets for {chain}, got {len(folds)}.")
        raw = load_snapshot_tables(args.cache / "raw", chain, metadata)
        features = dict(cached_fractions(
            root, args.cache, args.output, chain, raw, metadata, folds, args.workers,
        ))
        print(f"{chain}: scoring pre-existing {len(folds)} triplets", flush=True)
        ranks, predictions, coefficients = evaluate(chain, panel, folds, features)
        all_ranks.append(ranks)
        all_predictions.append(predictions)
        all_coefficients.append(coefficients)
    ranks = pd.concat(all_ranks, ignore_index=True)
    predictions = pd.concat(all_predictions, ignore_index=True)
    if args.reference_results:
        src = args.reference_results
        old_predictions = pd.read_csv(src / "internal_predictions.csv")
        old_ranks = pd.read_csv(src / "internal_ranking_by_mouse_split.csv")
        check_reference(predictions, old_predictions)
        old_predictions = old_predictions.loc[old_predictions["model"].isin(
            (BINARY_MODEL, CONTINUOUS_MODEL)
        )]
        old_ranks = old_ranks.loc[old_ranks["model"].isin((BINARY_MODEL, CONTINUOUS_MODEL))]
        predictions = pd.concat([predictions, old_predictions], ignore_index=True)
        ranks = pd.concat([ranks, old_ranks], ignore_index=True)
    ranks.to_csv(args.output / "internal_ranking_by_mouse_split.csv", index=False)
    predictions.to_csv(args.output / "internal_predictions.csv", index=False)
    ranking_summary(
        ranks, ["chain", "model", "top_k"], recovery_column=JOINT_VALUE,
        contrast_column=JOINT_CONTRAST,
    ).to_csv(args.output / "internal_ranking_summary.csv", index=False)
    prediction_summary(
        predictions, ["chain", "model", "outcome", "penalty"], repeated_splits=True,
    ).to_csv(args.output / "internal_error_summary.csv", index=False)
    pd.concat(all_coefficients, ignore_index=True).to_csv(
        args.output / "frozen_coefficients.csv", index=False,
    )
    (args.output / "run.json").write_text(json.dumps({
        "structural_seed": 1249, "calibration_seed": 1031, "depth": 15000,
        "n_triplets_per_chain": 75, "penalties": list(PENALTIES),
        "workers": args.workers, "new_rarefaction": False,
        "reuse_original_evidence": True, "reference_results": str(args.reference_results),
        "fraction_denominator": "tested_aaV_per_method_within_V",
        "positive_significance": "BH_FDR<0.05 and signed_log2FC>0",
    }, indent=2) + "\n")
    print(f"Done. New and comparison tables: {args.output}", flush=True)


if __name__ == "__main__":
    main()
