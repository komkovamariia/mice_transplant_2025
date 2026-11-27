#!/usr/bin/env python3
"""Evaluate frozen aaV-fraction ridge coefficients on the two public panels."""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_continuous import JOINT_VALUE, save
from scripts.run_eq7_weighted import external_validation
from src.eq7_aav_fraction import FRACTION, FRACTION_MODEL, fraction_evidence, prepare_aav_counts
from src.eq7_continuous import JOINT_CONTRAST, JointRecoveryPanel
from src.eq7_weighted import RidgeFit


def read_fits(directory, chain):
    table = pd.read_csv(directory / "frozen_coefficients.csv")
    selected = table.loc[table.chain.eq(chain) & table.model.eq(FRACTION_MODEL)]
    fits = {}
    for penalty, rows in selected.groupby("penalty"):
        weights = rows.set_index("term")["weight"]
        if set(weights.index) != set(FRACTION) or len(weights) != len(FRACTION):
            raise ValueError(f"Invalid frozen terms for {chain}, penalty {penalty}.")
        if rows.ridge_intercept.nunique() != 1 or rows.constant_baseline.nunique() != 1:
            raise ValueError(f"Inconsistent frozen intercept or baseline for {chain}, penalty {penalty}.")
        fits[penalty] = RidgeFit(weights, rows.ridge_intercept.iat[0], rows.constant_baseline.iat[0])
    if set(fits) != {1, 10}:
        raise ValueError(f"Expected frozen fits at penalties 1 and 10 for {chain}.")
    return {FRACTION_MODEL: fits}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--fraction-results", required=True, type=Path)
    parser.add_argument("--article-inputs", required=True, type=Path)
    parser.add_argument("--skin", type=Path)
    parser.add_argument("--marrow", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path, default=Path("results/eq7_aav_fraction_external_1249"))
    args = parser.parse_args()
    if not (args.skin or args.marrow):
        parser.error("Supply at least one study directory with TRUST4 reports.")
    if args.workers < 1:
        parser.error("--workers must be positive.")
    keys = ("structural_seed", "calibration_seed", "depth")
    for path in (args.cache / "v_input_seeds.json", args.fraction_results / "run.json"):
        data = json.loads(path.read_text())
        if tuple(data[k] for k in keys) != (1249, 1031, 15000):
            parser.error(f"Unexpected input seeds or depth in {path}.")
    if data.get("fisher_effect") != "log2_odds_ratio_0.5":
        parser.error("Refit the aaV-fraction model using Fisher log2 odds ratio first.")
    args.outcome = JOINT_CONTRAST
    root = Path(__file__).resolve().parents[1]
    metadata = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})

    def augment(candidates, tables, mice):
        counts, units = prepare_aav_counts(tables, mice)
        effects = fraction_evidence(counts, units, workers=args.workers)
        rows = candidates.merge(effects, on="v_gene", how="left", validate="one_to_one")
        if rows[list(FRACTION)].isna().any().any():
            raise ValueError("Matched internal predictors are missing for a frozen candidate.")
        return rows

    rankings, predictions = [], []
    for chain in (("TRA", "TRB") if args.marrow else ("TRB",)):
        ranked, scored = external_validation(
            root, args.cache, chain, metadata, args, read_fits(args.fraction_results, chain),
            panel_type=JointRecoveryPanel, augment=augment, model_names={FRACTION_MODEL: FRACTION_MODEL},
            recovery_column=JOINT_VALUE, rankers={},
        )
        if scored:
            rankings.append(ranked.loc[ranked.model.eq(FRACTION_MODEL)])
            predictions.extend(scored)
    args.output.mkdir(parents=True, exist_ok=True)
    save(args.output, "external", rankings, predictions, ["study", "chain"])
    print(f"External evaluation written to {args.output}")


if __name__ == "__main__":
    main()
