#!/usr/bin/env python3
"""Compare binary and continuous predictors against the same joint recovery target."""

import argparse
import hashlib
import inspect
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_weighted import (
    CHAINS, ORDERS, check_input_recovery, external_validation, load_inputs,
    prediction_summary, ranking_rows, ranking_summary,
)
from src.differential_abundance import fit_edger, fit_fisher
from src.eq7_calibration import FEATURES
from src.eq7_continuous import (
    BINARY_MODEL, CONTINUOUS, CONTINUOUS_MODEL, JOINT_CONTRAST, JointRecoveryPanel,
    continuous_evidence, prepare_v_evidence,
)
from src.eq7_weighted import GROUPS, contrasts, fit_ridge, split_shortlist, training_contrasts
from src.spike_in import load_snapshot_tables


JOINT_VALUE = "joint_recovery"
PENALTIES = (1.0, 10.0)
MODEL_TERMS = {BINARY_MODEL: FEATURES, CONTINUOUS_MODEL: CONTINUOUS}
RANKERS = {
    "structural_score": "n_g1_exclusive_aav",
    "continuous_edger": CONTINUOUS[0],
    "continuous_fisher": CONTINUOUS[1],
}


def attach_effects(rows, effects):
    augmented = rows.merge(effects, on="v_gene", how="left", validate="one_to_one")
    if augmented[list(CONTINUOUS[:2])].isna().any().any():
        raise ValueError("A structural candidate is missing from the focal-group count table.")
    return augmented


def split_mice(rows):
    mice = rows[["heldout_group", "heldout_mouse"]].drop_duplicates().set_index("heldout_group")
    return tuple((group, mice.at[group, "heldout_mouse"]) for group in GROUPS)


def fold_effects(cache, output, chain, raw, metadata, folds):
    from rpy2.robjects import r

    snapshot = cache / "raw" / "count_snapshots" / f"{chain}_manifest.json"
    digest = hashlib.sha256(snapshot.read_bytes())
    for function in (prepare_v_evidence, fit_edger, fit_fisher, continuous_evidence):
        digest.update(inspect.getsource(function).encode())
    digest.update(str(r('as.character(packageVersion("edgeR"))')[0]).encode())
    directory = output / "effect_cache" / chain / digest.hexdigest()
    directory.mkdir(parents=True, exist_ok=True)
    counts, units = prepare_v_evidence(raw, metadata)
    for index, (name, fold) in enumerate(folds.items(), start=1):
        mice = split_mice(fold)
        suffix = hashlib.sha256(repr(mice).encode()).hexdigest()[:16]
        path = directory / f"{name}_{suffix}.csv"
        if not path.is_file():
            effects = continuous_evidence(counts, units, mice)
            staged = path.with_suffix(".tmp")
            effects.to_csv(staged, index=False)
            staged.replace(path)
        if index % 15 == 0:
            print(f"{chain}: quantitative evidence {index}/{len(folds)}", flush=True)
        yield name, pd.read_csv(path)


def predictions(rows, estimate, model, penalty, observed):
    return pd.DataFrame({
        "v_gene": rows["v_gene"], "model": model, "outcome": JOINT_CONTRAST,
        "penalty": penalty, "observed": observed, "prediction": estimate.predict(rows),
        "constant_baseline": estimate.baseline,
    })


def evaluate(chain, panel, folds, effects, device):
    lookup = {key: idx for idx, key in enumerate(panel.mice)}
    inner, fits, ranks, refits = {}, [], [], []
    ablations = {f"continuous_without_{name}": tuple(term for term in CONTINUOUS if term != omitted)
                 for name, omitted in zip(("edgeR", "Fisher", "structure"), CONTINUOUS, strict=True)}
    for name, fold in folds.items():
        mice = split_mice(fold)
        left_out = tuple(lookup[mouse] for mouse in mice)
        values = panel.measure(panel.reference(left_out), left_out)
        check_input_recovery(fold, panel, left_out, values)
        eligible = attach_effects(split_shortlist(fold), effects[name])
        indices = panel.genes.get_indexer(eligible["v_gene"])
        if np.any(indices < 0) or not len(indices):
            raise ValueError(f"{chain} {name} has an invalid shortlist.")
        observed = contrasts(values[:, :, indices])[1]
        outcomes = training_contrasts(panel, left_out, inner, recovery_column=JOINT_CONTRAST)
        train = eligible.merge(outcomes[["v_gene", JOINT_CONTRAST]], on="v_gene", validate="one_to_one")
        refits.append(eligible.assign(**{JOINT_CONTRAST: observed}, fit_weight=1 / len(eligible)))
        order = {model: eligible[term].to_numpy() for model, term in RANKERS.items()}
        for model, features in {**MODEL_TERMS, **ablations}.items():
            for penalty in PENALTIES if model in MODEL_TERMS else (1.0,):
                estimate = fit_ridge(train, features, JOINT_CONTRAST, penalty, device=device)
                pred = predictions(eligible, estimate, model, penalty, observed)
                fits.append(pred.assign(chain=chain, split=name, **{f"{g}_mouse": m for g, m in mice}))
                if penalty == 1:
                    order[model] = pred["prediction"].to_numpy()
        for row in ranking_rows(eligible, values[:, :, indices], mice, ORDERS[chain], order,
                                recovery_column=JOINT_VALUE):
            ranks.append({"chain": chain, "split": name, **row})

    full_train = pd.concat(refits, ignore_index=True)
    frozen = {model: {p: fit_ridge(full_train, terms, JOINT_CONTRAST, p, weight="fit_weight", device=device)
                      for p in PENALTIES} for model, terms in MODEL_TERMS.items()}
    return pd.DataFrame(ranks), pd.concat(fits, ignore_index=True), frozen


def save(output, prefix, ranks, scored, keys):
    ranks = pd.concat(ranks, ignore_index=True)
    scored = pd.concat(scored, ignore_index=True)
    rank = ranking_summary(ranks, [*keys, "model", "top_k"], recovery_column=JOINT_VALUE,
                           contrast_column=JOINT_CONTRAST)
    errors = prediction_summary(scored, [*keys, "model", "outcome", "penalty"],
                                repeated_splits=prefix == "internal")
    for suffix, table in (("ranking_by_mouse_split", ranks), ("predictions", scored),
                          ("ranking_summary", rank), ("error_summary", errors)):
        table.to_csv(output / f"{prefix}_{suffix}.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=Path("results/eq7_continuous_1249"))
    parser.add_argument("--article-inputs", type=Path)
    parser.add_argument("--skin", type=Path, help="PRJNA954716 directory with trust4_reports.")
    parser.add_argument("--marrow", type=Path, help="PRJNA650133 directory with trust4_reports.")
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
    parser.add_argument("--prepare", choices=CHAINS, help="Only prepare evidence for the specified chain.")
    args = parser.parse_args()
    if (args.skin or args.marrow) and not args.article_inputs:
        parser.error("External evaluation requires --article-inputs.")
    args.outcome = JOINT_CONTRAST
    provenance = json.loads((args.cache / "v_input_seeds.json").read_text())
    if tuple(provenance[key] for key in ("structural_seed", "calibration_seed", "depth")) != (1249, 1031, 15000):
        parser.error("Use the manuscript's 1249/1031/15000 input cache.")
    root = Path(__file__).resolve().parents[1]
    metadata = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})
    args.output.mkdir(parents=True, exist_ok=True)
    ranks, scored, external_ranks, external_scored, coefficients = [], [], [], [], []
    for chain in (args.prepare,) if args.prepare else CHAINS:
        base, folds = load_inputs(root, args.cache, chain, metadata, panel_type=JointRecoveryPanel)
        raw = load_snapshot_tables(args.cache / "raw", chain, metadata)
        effects = dict(fold_effects(args.cache, args.output, chain, raw, metadata, folds))
        if args.prepare:
            continue
        print(f"{chain}: fitting and scoring {len(folds)} held-out triplets", flush=True)
        rank, pred, frozen = evaluate(chain, base, folds, effects, args.device)
        ranks.append(rank)
        scored.append(pred)
        for model, penalties in frozen.items():
            for penalty, fit in penalties.items():
                for term, weight in fit.coefficients.items():
                    coefficients.append((chain, model, penalty, term, weight, fit.intercept, fit.baseline))
        del base, folds, effects
        if args.skin or args.marrow:
            print(f"{chain}: evaluating available external cohorts", flush=True)
            external_rank, external_predictions = external_validation(
                root, args.cache, chain, metadata, args, frozen, panel_type=JointRecoveryPanel,
                augment=lambda candidates, tables, meta: attach_effects(
                    candidates, continuous_evidence(*prepare_v_evidence(tables, meta))),
                model_names={model: model for model in MODEL_TERMS},
                recovery_column=JOINT_VALUE, rankers=RANKERS, raw=raw,
            )
            if external_predictions:
                external_ranks.append(external_rank)
                external_scored.extend(external_predictions)
        del raw
    if not args.prepare:
        save(args.output, "internal", ranks, scored, ["chain"])
        if external_scored:
            save(args.output, "external", external_ranks, external_scored, ["study", "chain"])
        pd.DataFrame(coefficients, columns=[
            "chain", "model", "penalty", "term", "weight", "ridge_intercept", "constant_baseline",
        ]).to_csv(args.output / "frozen_coefficients.csv", index=False)
        (args.output / "run.json").write_text(json.dumps({
            **provenance, "device": args.device, "fixed_penalties": PENALTIES,
            "external_skin": args.skin is not None, "external_marrow": args.marrow is not None,
        }, indent=2) + "\n")


if __name__ == "__main__":
    main()
