#!/usr/bin/env python3
"""Fisher odds-ratio direction sensitivity and external aaV-fraction validation.

Requires existing 1249/1031/15000 inputs and completed aaV-fraction effect cache.
The 75 held-out triplets and original edgeR values are re-used unchanged.
Only Fisher statistics with the odds-ratio direction are newly computed.
External outcomes and structural candidates are read from the published
Overleaf analysis archive; no FASTQ, TRUST4 or external rarefaction is run.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_aav_fraction import PENALTIES, attach_fractions
from scripts.run_eq7_continuous import (
    JOINT_CONTRAST, JOINT_VALUE, ORDERS, check_input_recovery,
    load_inputs, prediction_summary, ranking_rows, ranking_summary, split_mice,
)
from src.eq7_aav_fraction import FRACTION, FRACTION_MODEL, fraction_evidence, prepare_aav_counts
from src.eq7_continuous import JointRecoveryPanel
from src.eq7_fisher_or import ODDS_FISHER, RAW_FISHER, fisher_fraction_sensitivity
from src.eq7_weighted import contrasts, fit_ridge, split_shortlist, top_weights, training_contrasts
from src.spike_in import load_snapshot_tables
from src.structural_v import subset_v_tables


ODDS_MODEL = "aav_fisher_odds_ridge"
ODDS_FEATURES = (FRACTION[0], ODDS_FISHER, FRACTION[2])
ABLATIONS = {
    "odds_without_edgeR": ODDS_FEATURES[1:],
    "odds_without_Fisher": (ODDS_FEATURES[0], ODDS_FEATURES[2]),
    "odds_without_structure": ODDS_FEATURES[:2],
}
MODELS = {FRACTION_MODEL: FRACTION, ODDS_MODEL: ODDS_FEATURES, **ABLATIONS}
RANKS = {
    "structural_score": "n_g1_exclusive_aav",
    "aaV_edgeR_fraction": FRACTION[0],
    "aaV_Fisher_raw_fraction": RAW_FISHER,
    "aaV_Fisher_odds_fraction": ODDS_FISHER,
}


def old_fraction_file(previous: Path, chain: str, split: str, mice):
    label = hashlib.sha256(repr(mice).encode()).hexdigest()[:16]
    matches = list((previous / "effect_cache" / chain).glob(f"*/{split}_{label}.csv"))
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one completed original feature file for {chain} {split}; found {len(matches)}.")
    return matches[0]


def fold_fractions(cache, previous, output, chain, folds, counts, units, workers):
    root = Path(__file__).resolve().parents[1]
    signature = hashlib.sha256()
    signature.update((cache / "raw" / "count_snapshots" / f"{chain}_manifest.json").read_bytes())
    signature.update((root / "src" / "eq7_fisher_or.py").read_bytes())
    signature.update(b"fisher-haldane-anscombe-odds-direction-2026-11")
    directory = output / "fisher_cache" / chain / signature.hexdigest()
    directory.mkdir(parents=True, exist_ok=True)
    for i, (split, fold) in enumerate(folds.items(), 1):
        mice = split_mice(fold)
        label = hashlib.sha256(repr(mice).encode()).hexdigest()[:16]
        old = pd.read_csv(old_fraction_file(previous, chain, split, mice))
        cache_file = directory / f"{split}_{label}.csv"
        if cache_file.is_file():
            current = pd.read_csv(cache_file)
            status = "cached"
        else:
            current = fisher_fraction_sensitivity(counts, units, mice, workers=workers)
            staged = cache_file.with_suffix(".tmp")
            current.to_csv(staged, index=False)
            staged.replace(cache_file)
            status = "computed"
        merged = old.merge(
            current[[ "v_gene", RAW_FISHER, ODDS_FISHER, "n_Fisher_odds_positive"]],
            on="v_gene", validate="one_to_one", suffixes=("_original", "_recomputed"),
        )
        if len(merged) != len(old):
            raise ValueError(f"Missing Fisher fraction genes in {chain} {split}")
        np.testing.assert_allclose(
            merged[f"{RAW_FISHER}_original"], merged[f"{RAW_FISHER}_recomputed"],
            rtol=0, atol=1e-12,
            err_msg=f"Recomputed raw Fisher fraction disagrees with old cache: {chain} {split}",
        )
        merged[RAW_FISHER] = merged.pop(f"{RAW_FISHER}_original")
        merged.drop(columns=[f"{RAW_FISHER}_recomputed"], inplace=True)
        if i % 5 == 0 or i == 1:
            print(f"{chain}: Fisher odds direction {i}/{len(folds)} ({status})", flush=True)
        yield split, merged


def evaluate_chain(chain, panel, folds, features):
    lookup = {mouse: idx for idx, mouse in enumerate(panel.mice)}
    inner_measures, rankings, predictions, final_rows = {}, [], [], []
    for split, fold in folds.items():
        mice = split_mice(fold)
        left_out = tuple(lookup[m] for m in mice)
        observations = panel.measure(panel.reference(left_out), left_out)
        check_input_recovery(fold, panel, left_out, observations)
        eligible = attach_fractions(split_shortlist(fold), features[split])
        gene_ids = panel.genes.get_indexer(eligible["v_gene"])
        if not len(gene_ids) or (gene_ids < 0).any():
            raise ValueError(f"Missing structural V gene in {chain} {split}")
        truth = contrasts(observations[:, :, gene_ids])[1]
        training = training_contrasts(
            panel, left_out, inner_measures, recovery_column=JOINT_CONTRAST,
        )
        training = eligible.merge(training[["v_gene", JOINT_CONTRAST]], on="v_gene", validate="one_to_one")
        final_rows.append(eligible.assign(**{
            JOINT_CONTRAST: truth,
            "fit_weight": 1.0 / len(eligible),
        }))
        rank_scores = {name: eligible[column].to_numpy(dtype=float) for name, column in RANKS.items()}
        for model, terms in MODELS.items():
            for penalty in PENALTIES if model in (FRACTION_MODEL, ODDS_MODEL) else (1.0,):
                fit = fit_ridge(training, terms, JOINT_CONTRAST, penalty)
                predicted = fit.predict(eligible)
                predictions.append(pd.DataFrame({
                    "chain": chain, "split": split, "v_gene": eligible["v_gene"],
                    "g1_mouse": mice[0][1], "g5_mouse": mice[1][1], "g6_mouse": mice[2][1],
                    "model": model, "outcome": JOINT_CONTRAST, "penalty": penalty,
                    "observed": truth, "prediction": predicted,
                    "constant_baseline": fit.baseline,
                }))
                if penalty == 1.0:
                    rank_scores[model] = predicted
        for row in ranking_rows(
            eligible, observations[:, :, gene_ids], mice, ORDERS[chain],
            rank_scores, recovery_column=JOINT_VALUE,
        ):
            rankings.append({"chain": chain, "split": split, **row})
    all_training = pd.concat(final_rows, ignore_index=True)
    frozen = {
        model: {penalty: fit_ridge(
            all_training, terms, JOINT_CONTRAST, penalty, weight="fit_weight",
        ) for penalty in PENALTIES}
        for model, terms in ((FRACTION_MODEL, FRACTION), (ODDS_MODEL, ODDS_FEATURES))
    }
    weights = []
    for name, fits in frozen.items():
        for penalty, fit in fits.items():
            for term, value in fit.coefficients.items():
                weights.append({
                    "chain": chain, "model": name, "penalty": penalty, "term": term,
                    "weight": value, "ridge_intercept": fit.intercept,
                    "constant_baseline": fit.baseline,
                })
    return pd.DataFrame(rankings), pd.concat(predictions, ignore_index=True), frozen, pd.DataFrame(weights)


def verify_old_model(new, previous, old_frozen):
    old = pd.read_csv(previous / "internal_predictions.csv")
    select = lambda x: x.loc[x["model"].eq(FRACTION_MODEL)]
    keys = ["chain", "split", "v_gene", "model", "penalty"]
    merged = select(new).merge(
        select(old), on=keys, validate="one_to_one", indicator=True,
        suffixes=("_new", "_prior"),
    )
    if not merged["_merge"].eq("both").all() or len(merged) != len(select(old)):
        raise ValueError("The old and new models differ in held-out split/candidate keys.")
    for name in ("observed", "prediction", "constant_baseline"):
        np.testing.assert_allclose(
            merged[f"{name}_new"], merged[f"{name}_prior"], rtol=0, atol=1e-9,
            err_msg=f"Changed {name} for old model; aborting comparison.",
        )
    saved = pd.read_csv(previous / "frozen_coefficients.csv")
    check = old_frozen.merge(
        saved, on=["chain", "model", "penalty", "term"],
        suffixes=("_new", "_prior"), validate="one_to_one", indicator=True,
    )
    if not check["_merge"].eq("both").all() or len(check) != len(saved):
        raise ValueError("Frozen coefficients do not match the original run.")
    for col in ("weight", "ridge_intercept", "constant_baseline"):
        np.testing.assert_allclose(
            check[f"{col}_new"], check[f"{col}_prior"], rtol=0, atol=1e-8,
        )
    return len(merged)


def archive_csv(archive, name):
    name = f"analysis_outputs/{name}"
    try:
        data = archive.read(name)
    except KeyError as exc:
        raise FileNotFoundError(f"Archive missing {name}") from exc
    return pd.read_csv(io.BytesIO(data))


def external_check(archive_path, raw_by_chain, metadata, frozen_by_chain, workers):
    predictions, ranking = [], []
    with zipfile.ZipFile(archive_path) as archive:
        reference = archive_csv(archive, "eq7_continuous_1249/external_predictions.csv")
        for study, chain, subsets in (
            ("PRJNA954716", "TRB", ("cd4",)),
            ("PRJNA650133", "TRA", ("cd4", "cd8")),
            ("PRJNA650133", "TRB", ("cd4", "cd8")),
        ):
            suffix = f"/{chain}" if study == "PRJNA650133" else ""
            candidates = archive_csv(
                archive, f"{study}/validation{suffix}/train_reference_scores.csv",
            )
            matched = metadata.loc[
                metadata["_tissue"].eq("spleen")
                & metadata["_cell_subset"].isin(subsets)
                & metadata["group_no"].isin(("g1", "g5"))
            ]
            raw = raw_by_chain[chain]
            training_tables = subset_v_tables(
                {group: raw[group] for group in ("g1", "g5")}, matched["sample_id"],
            )
            pooled, units = prepare_aav_counts(training_tables, matched)
            effects = fraction_evidence(pooled, units, (), workers=workers)
            changed = fisher_fraction_sensitivity(pooled, units, (), workers=workers)
            both = effects.merge(
                changed[["v_gene", RAW_FISHER, ODDS_FISHER]],
                on="v_gene", suffixes=("_original", "_recomputed"), validate="one_to_one",
            )
            np.testing.assert_allclose(
                both[f"{RAW_FISHER}_original"], both[f"{RAW_FISHER}_recomputed"],
                atol=1e-12, rtol=0,
            )
            both[RAW_FISHER] = both.pop(f"{RAW_FISHER}_original")
            both.drop(columns=f"{RAW_FISHER}_recomputed", inplace=True)
            candidates = attach_fractions(candidates, both)
            old = reference.loc[
                reference["study"].eq(study) & reference["chain"].eq(chain)
                & reference["model"].eq("continuous_evidence_ridge")
                & reference["penalty"].eq(1.0)
            ].copy()
            observed = old.set_index("v_gene")["observed"]
            status = old.set_index("v_gene")["status"]
            if set(observed.index) != set(candidates["v_gene"]) or observed.index.duplicated().any():
                raise ValueError(f"External candidate mismatch: {study}/{chain}")
            for model, models in frozen_by_chain[chain].items():
                for penalty, fit in models.items():
                    current = pd.DataFrame({
                        "study": study, "chain": chain, "v_gene": candidates["v_gene"],
                        "model": model, "outcome": JOINT_CONTRAST, "penalty": penalty,
                        "observed": observed.reindex(candidates["v_gene"]).to_numpy(),
                        "prediction": fit.predict(candidates),
                        "constant_baseline": fit.baseline,
                        "status": status.reindex(candidates["v_gene"]).to_numpy(),
                    })
                    predictions.append(current)
                    if study == "PRJNA954716" and penalty == 1:
                        score = current["prediction"].to_numpy()
                        truth = current["observed"].to_numpy()
                        for k in ORDERS[chain]:
                            weight = top_weights(score, k)
                            ranking.append({
                                "study": study, "chain": chain, "model": model, "top_k": k,
                                "observed_joint_contrast": float(truth @ weight / weight.sum()),
                            })
            if study == "PRJNA954716":
                truth = observed.reindex(candidates["v_gene"]).to_numpy()
                for label, values in (
                    ("structural_score", candidates["n_g1_exclusive_aav"].to_numpy(dtype=float)),
                    ("aaV_Fisher_raw_fraction", candidates[RAW_FISHER].to_numpy()),
                    ("aaV_Fisher_odds_fraction", candidates[ODDS_FISHER].to_numpy()),
                ):
                    for k in ORDERS[chain]:
                        weight = top_weights(values, k)
                        ranking.append({
                            "study": study, "chain": chain, "model": label, "top_k": k,
                            "observed_joint_contrast": float(truth @ weight / weight.sum()),
                        })
    return pd.concat(predictions, ignore_index=True), pd.DataFrame(ranking)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--previous-results", type=Path, required=True,
                        help="Completed original results/eq7_aav_fraction_1249 folder.")
    parser.add_argument("--output", type=Path, default=Path("results/eq7_fisher_or_external_1249"))
    parser.add_argument("--article-archive", type=Path,
                        help="Overleaf zip containing matched reference scores and continuous external outcomes.")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Workers must be >= 1.")
    seeds_file = args.cache / "v_input_seeds.json"
    if not seeds_file.is_file():
        parser.error("Missing EXISTING prepared _inputs; this script never prepares them.")
    seeds = json.loads(seeds_file.read_text())
    if tuple(seeds[k] for k in ("structural_seed", "calibration_seed", "depth")) != (1249, 1031, 15000):
        parser.error("Only the 1249/1031/15000 prepared cache is supported.")
    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    metadata = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})
    ranks, predictions, coefficient_tables, all_frozen = [], [], [], {}
    raw_by_chain = {}
    for chain in ("TRA", "TRB"):
        panel, folds = load_inputs(root, args.cache, chain, metadata, panel_type=JointRecoveryPanel)
        if len(folds) != 75:
            raise ValueError(f"Expected 75 recorded held-out triplets for {chain}")
        raw = load_snapshot_tables(args.cache / "raw", chain, metadata)
        raw_by_chain[chain] = raw
        counts, units = prepare_aav_counts(raw, metadata)
        features = dict(fold_fractions(
            args.cache, args.previous_results, args.output, chain, folds,
            counts, units, args.workers,
        ))
        print(f"{chain}: fitting both variants on the SAME 75 old splits", flush=True)
        rank, pred, frozen, coeff = evaluate_chain(chain, panel, folds, features)
        ranks.append(rank)
        predictions.append(pred)
        coefficient_tables.append(coeff)
        all_frozen[chain] = frozen
    rank = pd.concat(ranks, ignore_index=True)
    pred = pd.concat(predictions, ignore_index=True)
    coefficients = pd.concat(coefficient_tables, ignore_index=True)
    checked = verify_old_model(
        pred, args.previous_results,
        coefficients.loc[coefficients["model"].eq(FRACTION_MODEL)],
    )
    print(f"Confirmed unchanged original predictions: {checked} rows", flush=True)
    pred.to_csv(args.output / "internal_predictions.csv", index=False)
    rank.to_csv(args.output / "internal_ranking_by_mouse_split.csv", index=False)
    coefficients.to_csv(args.output / "frozen_coefficients.csv", index=False)
    prediction_summary(
        pred, ["chain", "model", "outcome", "penalty"], repeated_splits=True,
    ).to_csv(args.output / "internal_error_summary.csv", index=False)
    ranking_summary(
        rank, ["chain", "model", "top_k"], recovery_column=JOINT_VALUE,
        contrast_column=JOINT_CONTRAST,
    ).to_csv(args.output / "internal_ranking_summary.csv", index=False)
    if args.article_archive:
        print("External validation from matched original reference scores and stored outcomes", flush=True)
        external, external_ranks = external_check(
            args.article_archive, raw_by_chain, metadata, all_frozen, args.workers,
        )
        external.to_csv(args.output / "external_predictions.csv", index=False)
        external_ranks.to_csv(args.output / "external_topk_summary.csv", index=False)
        prediction_summary(
            external, ["study", "chain", "model", "outcome", "penalty"],
        ).to_csv(args.output / "external_error_summary.csv", index=False)
    (args.output / "run.json").write_text(json.dumps({
        "structural_seed": 1249, "calibration_seed": 1031, "depth": 15000,
        "n_old_triplets_per_chain": 75,
        "new_rarefaction": False, "new_splits": False, "new_edgeR_per_fold": False,
        "new_fisher_per_fold": True, "workers": args.workers,
        "previous_results": str(args.previous_results),
        "article_archive": str(args.article_archive) if args.article_archive else None,
        "external_uses_archived_outcomes": bool(args.article_archive),
    }, indent=2) + "\n")
    print(f"Finished sensitivity and external comparison: {args.output}", flush=True)


if __name__ == "__main__":
    main()
