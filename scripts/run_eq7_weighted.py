#!/usr/bin/env python3
"""Compare Equation 7 with count and abundance-weighted recovery contrasts."""

from __future__ import annotations

import argparse
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_external import external_raw
from scripts.run_v_updates import calibration_dir, evidence_dir
from src.eq7_calibration import FEATURES
from src.eq7_weighted import (
    COUNT_CONTRAST, GROUPS, MOLECULE_CONTRAST, RecoveryPanel, contrasts,
    fit_ridge, full_training_profiles, pool_mice, split_shortlist, top_weights,
    training_contrasts,
)
from src.spike_in import load_snapshot_tables
from src.structural_sampling import downsample_group_tables
from src.structural_v import subset_v_tables


CHAINS = ("TRA", "TRB")
ORDERS = {"TRA": (5, 10), "TRB": (1, 2)}
PENALTIES = (1.0, 10.0)
MOLECULE_MODEL = "weighted_contrast_ridge"
COUNT_MODEL = "count_contrast_ridge"
RANKERS = {
    "structural_score": "n_g1_exclusive_aav",
    "edger_indicator": FEATURES[0],
    "fisher_indicator": FEATURES[1],
}


def fold_path(directory, parent, filename):
    if directory is not None and (directory / filename).is_file():
        return directory
    matches = [path for path in parent.iterdir() if (path / filename).is_file()]
    if len(matches) != 1:
        raise ValueError(f"Expected one cache containing {filename} under {parent}.")
    return matches[0]


def load_inputs(root, cache, chain, metadata, *, panel_type=RecoveryPanel):
    expected_calibration = calibration_dir(root, cache, chain)
    calibration = fold_path(
        expected_calibration, cache / "calibration" / chain, f"count_snapshots/{chain}_manifest.json",
    )
    tables = load_snapshot_tables(calibration, chain, metadata)
    expected_evidence = evidence_dir(root, cache, chain) if calibration == expected_calibration else None
    evidence = fold_path(
        expected_evidence, cache / "eq7_evidence" / chain, "stratified_fold_evidence.csv",
    )
    rows = pd.read_csv(
        evidence / "stratified_fold_evidence.csv", dtype={"heldout_mouse": str},
    )
    panel = panel_type(pool_mice(tables, metadata))
    folds = dict(tuple(rows.groupby("evidence_split")))
    by_group = panel.by_group()
    seen = {
        tuple(fold.groupby("heldout_group")["heldout_mouse"].first().reindex(GROUPS))
        for fold in folds.values()
    }
    expected = set(product(*([panel.mice[i][1] for i in by_group[group]] for group in GROUPS)))
    if len(folds) != len(expected) or seen != expected:
        raise ValueError(f"The {chain} evidence must contain one fold for each eligible mouse triplet.")
    return panel, folds


def model_specs():
    yield "legacy_score_refit", FEATURES, "signed_recovery_pct", 10.0, True
    for penalty in PENALTIES:
        yield COUNT_MODEL, FEATURES, COUNT_CONTRAST, penalty, False
        yield MOLECULE_MODEL, FEATURES, MOLECULE_CONTRAST, penalty, False
        yield "weighted_contrast_nonnegative", FEATURES, MOLECULE_CONTRAST, penalty, True
    for omitted, name in zip(FEATURES, ("edgeR", "Fisher", "structural_count"), strict=True):
        terms = tuple(term for term in FEATURES if term != omitted)
        yield f"weighted_without_{name}", terms, MOLECULE_CONTRAST, 1.0, False


def check_input_recovery(rows, panel, query, observed):
    for idx, mouse in enumerate(query):
        group, mouse_id = panel.mice[mouse]
        frozen = rows.loc[rows["heldout_mouse"].eq(mouse_id)].set_index("v_gene")["response"]
        reconstructed = pd.Series(observed[idx, 0], index=panel.genes).reindex(frozen.index)
        np.testing.assert_allclose(reconstructed, (100 if group == "g1" else -100) * frozen, atol=1e-12)


def ranking_rows(rows, measured, mice, sizes, ranks, *, recovery_column="weighted_molecules_per_million"):
    for name, score in ranks.items():
        for size in sizes:
            keep = top_weights(score, size)
            for idx, (group, mouse) in enumerate(mice):
                clone, molecules = measured[idx]
                defined = np.isfinite(clone)
                denom = keep[defined].sum()
                yield {
                    "model": name, "top_k": size, "group": group, "mouse": mouse,
                    "count_recovery_pct": clone[defined] @ keep[defined] / denom if denom else np.nan,
                    recovery_column: molecules @ keep / keep.sum(),
                    "effective_k": min(size, len(rows)),
                }


def evaluate_internal(chain, panel, folds, device):
    lookup = {key: idx for idx, key in enumerate(panel.mice)}
    measured_inner, ranks, predictions = {}, [], []
    for name, fold in folds.items():
        mice_by_group = fold[["heldout_group", "heldout_mouse"]].drop_duplicates().set_index("heldout_group")
        mice = [(group, mice_by_group.at[group, "heldout_mouse"]) for group in GROUPS]
        query = tuple(lookup[key] for key in mice)
        observed = panel.measure(panel.reference(query), query)
        check_input_recovery(fold, panel, query, observed)
        eligible = split_shortlist(fold)
        genes = panel.genes.get_indexer(eligible["v_gene"])
        if len(eligible) == 0 or np.any(genes < 0):
            raise ValueError(f"Missing structural candidates in {chain} {name}.")
        test = observed[:, :, genes]
        truth = dict(zip((COUNT_CONTRAST, MOLECULE_CONTRAST), contrasts(test), strict=True))
        training = eligible.merge(training_contrasts(panel, query, measured_inner), on="v_gene", validate="one_to_one")
        rankings = {model: eligible[term].to_numpy(dtype=float) for model, term in RANKERS.items()}
        for model, terms, target, penalty, nonnegative in model_specs():
            fit = fit_ridge(training, terms, target, penalty, nonnegative=nonnegative, device=device)
            predicted = fit.predict(eligible)
            if penalty == 1 or model == "legacy_score_refit":
                rankings[model] = predicted
            if target in truth:
                predictions.append(pd.DataFrame({
                    "chain": chain, "split": name, "v_gene": eligible["v_gene"],
                    "g1_mouse": mice[0][1], "g5_mouse": mice[1][1], "g6_mouse": mice[2][1],
                    "model": model, "outcome": target, "penalty": penalty, "observed": truth[target],
                    "prediction": predicted, "constant_baseline": fit.baseline,
                }))
        for row in ranking_rows(eligible, test, mice, ORDERS[chain], rankings):
            ranks.append({"chain": chain, "split": name, **row})
    return pd.DataFrame(ranks), pd.concat(predictions, ignore_index=True)


def external_validation(root, cache, chain, metadata, args, fit, *,
                        panel_type=RecoveryPanel, augment=None, model_names=None,
                        recovery_column="weighted_molecules_per_million", rankers=RANKERS, raw=None):
    if raw is None:
        raw = load_snapshot_tables(cache / "raw", chain, metadata)
    predictions, recoveries = [], []
    for study, directory, subsets in (
        ("PRJNA954716", args.skin, ("cd4",)),
        ("PRJNA650133", args.marrow, ("cd4", "cd8")),
    ):
        if directory is None or (study == "PRJNA954716" and chain != "TRB"):
            continue
        validation = args.article_inputs / study / "validation"
        if study == "PRJNA650133":
            validation /= chain
        candidates = pd.read_csv(validation / "train_reference_scores.csv")
        matched = metadata.loc[metadata["_tissue"].eq("spleen") & metadata["_cell_subset"].isin(subsets)
                               & metadata["group_no"].isin(("g1", "g5"))]
        train_tables = subset_v_tables({group: raw[group] for group in ("g1", "g5")}, matched["sample_id"])
        sampling = {"depth": 15000, "seed": 1031}
        internal, _ = downsample_group_tables(train_tables, **sampling)
        ref_panel = panel_type(pool_mice(internal, matched))
        external_meta = pd.read_csv(root / "config" / f"{study.lower()}_samples.csv", dtype={"mouse_id": str})
        counts, qc = external_raw(directory, external_meta, chain)
        expected_qc = pd.read_csv(validation / "external_library_qc.csv").set_index("sample_id")
        actual_qc = qc.set_index("sample_id")["report_sha256"]
        if not expected_qc.loc[actual_qc.index, "report_sha256"].eq(actual_qc).all():
            raise ValueError(f"{study} TRUST4 reports disagree with the recorded run checksums.")
        external, _ = downsample_group_tables(counts, **sampling)
        query = panel_type(pool_mice(external, external_meta), align_to=ref_panel)
        genes = ref_panel.genes.get_indexer(candidates["v_gene"])
        if np.any(genes < 0):
            raise ValueError(f"{study} candidate V segments disagree with the matched internal panel.")
        measured = query.measure(ref_panel.reference(), range(len(query.mice)))[:, :, genes]
        locked = pd.read_csv(validation / "external_mouse_predictions.csv")
        for idx, (group, mouse) in enumerate(query.mice):
            recorded = locked.loc[locked["heldout_mouse"].eq(mouse)].set_index("v_gene")["response"]
            np.testing.assert_allclose(
                measured[idx, 0], (100 if group == "g1" else -100) * recorded.reindex(candidates["v_gene"]),
                atol=1e-12,
            )
        if augment is not None:
            candidates = augment(candidates, train_tables, matched)
        ranks = {name: candidates[term].to_numpy() for name, term in rankers.items()}
        ranks["legacy_score_frozen"] = candidates["E"].to_numpy()
        for target, models in fit.items():
            for penalty, model in models.items():
                prediction = model.predict(candidates)
                model_name = (model_names[target] if model_names is not None else
                              COUNT_MODEL if target == COUNT_CONTRAST else MOLECULE_MODEL)
                if penalty == 1:
                    ranks[model_name] = prediction
                observations = np.full(len(candidates), np.nan)
                if study == "PRJNA954716":
                    offset = 0 if target == COUNT_CONTRAST else 1
                    groups = query.by_group()
                    observations = measured[list(groups["g1"]), offset].mean(axis=0) - measured[list(groups["g5"]), offset].mean(axis=0)
                predictions.append(pd.DataFrame({
                    "study": study, "chain": chain, "v_gene": candidates["v_gene"],
                    "model": model_name, "outcome": target if model_names is None else args.outcome,
                    "penalty": penalty, "observed": observations,
                    "prediction": prediction, "constant_baseline": model.baseline,
                    "status": "observed_two_arm_contrast" if study == "PRJNA954716" else "graft_arm_only",
                }))
        for row in ranking_rows(candidates, measured, query.mice, ORDERS[chain], ranks, recovery_column=recovery_column):
            recoveries.append({"study": study, "chain": chain, **row})
    return pd.DataFrame(recoveries), predictions


def ranking_summary(ranks, keys, *, recovery_column="weighted_molecules_per_million",
                    contrast_column=MOLECULE_CONTRAST):
    endpoints = ["count_recovery_pct", recovery_column]
    mice = ranks.groupby([*keys, "group", "mouse"])[endpoints].mean().reset_index()
    groups = mice.groupby([*keys, "group"])[endpoints].mean()
    results = []
    for key, frame in groups.groupby(level=keys):
        means = frame.droplevel(keys)
        row = dict(zip(keys, key, strict=True))
        for endpoint in endpoints:
            for group in GROUPS:
                row[f"{group}_{endpoint}"] = means[endpoint].get(group, np.nan)
        if {"g1", "g5"} <= set(means.index):
            ctl = means.loc[means.index != "g1"].mean()
            row[COUNT_CONTRAST], row[contrast_column] = means.loc["g1"] - ctl
        results.append(row)
    return pd.DataFrame(results)


def prediction_summary(predictions, keys, *, repeated_splits=False):
    scored = predictions.dropna(subset=["observed"]).copy()
    scored["mae"] = (scored["prediction"] - scored["observed"]).abs()
    scored["constant_mae"] = (scored["constant_baseline"] - scored["observed"]).abs()
    errors = scored.groupby([*keys, "split"] if repeated_splits else keys)[["mae", "constant_mae"]].mean()
    return errors.groupby(level=keys).mean().reset_index()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path, help="Existing manuscript input or benchmark _inputs directory.")
    parser.add_argument("--output", type=Path, default=Path("results/eq7_weighted_1249"))
    parser.add_argument("--article-inputs", type=Path, help="Overleaf analysis_outputs directory with existing external results.")
    parser.add_argument("--skin", type=Path, help="PRJNA954716 run directory containing trust4_reports.")
    parser.add_argument("--marrow", type=Path, help="PRJNA650133 run directory containing trust4_reports.")
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args()
    if (args.skin or args.marrow) and not args.article_inputs:
        parser.error("External evaluation requires --article-inputs.")
    if args.device == "cuda":
        try:
            import torch
        except ImportError:
            parser.error("--device cuda requires a CUDA-enabled PyTorch installation.")
        if not torch.cuda.is_available():
            parser.error("--device cuda requires a visible CUDA device.")

    root = Path(__file__).resolve().parents[1]
    provenance = json.loads((args.cache / "v_input_seeds.json").read_text())
    seeds = {key: provenance[key] for key in ("structural_seed", "calibration_seed", "depth")}
    if tuple(seeds.values()) != (1249, 1031, 15000):
        parser.error("Use the manuscript's 1249/1031/15000 input cache.")
    metadata = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})
    args.output.mkdir(parents=True, exist_ok=True)
    internal_ranks, internal_predictions, external_ranks, external_predictions = [], [], [], []
    frozen_coefficients = []
    for chain in CHAINS:
        panel, folds = load_inputs(root, args.cache, chain, metadata)
        print(f"{chain}: evaluating {len(folds)} held-out triplets", flush=True)
        ranks, predictions = evaluate_internal(chain, panel, folds, args.device)
        internal_ranks.append(ranks)
        internal_predictions.append(predictions)
        fits = full_training_profiles(panel, folds, PENALTIES, args.device)
        for outcome, models in fits.items():
            for penalty, model in models.items():
                for term, weight in model.coefficients.items():
                    frozen_coefficients.append({
                        "chain": chain, "outcome": outcome, "penalty": penalty, "term": term,
                        "weight": weight, "ridge_intercept": model.intercept, "constant_baseline": model.baseline,
                    })
        if args.skin or args.marrow:
            ranks, predictions = external_validation(root, args.cache, chain, metadata, args, fits)
            external_ranks.append(ranks)
            external_predictions.extend(predictions)
    for prefix, ranks, predictions, keys in (
        ("internal", internal_ranks, internal_predictions, ["chain"]),
        ("external", external_ranks, external_predictions, ["study", "chain"]),
    ):
        if not ranks:
            continue
        ranks = pd.concat(ranks, ignore_index=True)
        predictions = pd.concat(predictions, ignore_index=True)
        summaries = (
            ("ranking_by_mouse_split", ranks), ("predictions", predictions),
            ("ranking_summary", ranking_summary(ranks, [*keys, "model", "top_k"])),
            ("error_summary", prediction_summary(
                predictions, [*keys, "model", "outcome", "penalty"], repeated_splits=prefix == "internal",
            )),
        )
        for name, frame in summaries:
            frame.to_csv(args.output / f"{prefix}_{name}.csv", index=False)
    pd.DataFrame(frozen_coefficients).to_csv(args.output / "frozen_coefficients.csv", index=False)
    (args.output / "run.json").write_text(json.dumps({
        **seeds, "device": args.device, "fixed_penalties": PENALTIES,
        "external_skin": args.skin is not None, "external_marrow": args.marrow is not None,
    }, indent=2) + "\n")
    print(f"Equation 7 comparison written to {args.output}")


if __name__ == "__main__":
    main()
