#!/usr/bin/env python3
"""Compare seven equal-size V-segment selection rules in single-mouse holdouts."""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_weighted import fold_path
from scripts.run_v_updates import calibration_dir, evidence_dir
from src.differential_abundance import fit_edger
from src.eq7_calibration import shortlist
from src.eq7_continuous import prepare_v_evidence
from src.eq7_weighted import top_weights
from src.spike_in import EDGER_SUCCESS, aa_v_key, load_snapshot_tables


def selection_rows(chain, folds, raw, metadata):
    counts, units = prepare_v_evidence(raw, metadata)
    for split, rows in folds.groupby("evidence_split"):
        group, mouse = rows[["heldout_group", "heldout_mouse"]].iloc[0]
        k = len(shortlist(rows))
        omit = units.index[units["group_no"].eq(group) & units["mouse_id"].eq(mouse)]
        retained = units.drop(index=omit)
        test, status = fit_edger(counts.loc[:, retained.index], retained)
        if status != EDGER_SUCCESS:
            raise ValueError(f"{chain} {split}: {status}")
        test["v_gene"] = test["feature_id"].map(lambda key: aa_v_key(key)[1])
        v = rows.merge(test[["v_gene", "logFC", "FDR"]], on="v_gene", how="left", validate="1:1")
        fold_change = v["logFC"].fillna(0).clip(lower=0)
        significant = v["FDR"].lt(.05)
        rules = {
            "structural_screen": v["v_gene"].isin(shortlist(rows)["v_gene"]).astype(int),
            "S_ungated": v["n_g1_exclusive_aav"],
            "initial_share_ungated": v["n_full"],
            "edgeR_binary_ungated": v["edgeR_g1_enriched"],
            "Fisher_binary_ungated": v["fisher_g1_enriched"],
            "V_usage_edgeR_logFC": fold_change,
            "V_usage_edgeR_significant": (fold_change * -np.log10(v["FDR"].fillna(1).clip(lower=1e-300))).where(significant, 0),
        }
        observed = v["response"].notna().to_numpy()
        for method, scores in rules.items():
            weight = top_weights(scores.to_numpy(), k)
            recovery = v["response"].abs().fillna(0).to_numpy()
            denom = weight[observed].sum()
            yield dict(chain=chain, fold=rows.heldout_fold.iat[0], group=group, mouse=mouse,
                       method=method, k_train_structural=k, expected_observed_v=denom,
                       mean_raw_recovery=weight @ recovery / denom if denom else np.nan)


def contrast_summary(details):
    means = details.groupby(["chain", "method", "group"])["mean_raw_recovery"].mean().unstack()
    means["g1_minus_mean_g5g6"] = means["g1"] - (means["g5"] + means["g6"]) / 2
    return means.reset_index()


def paired_primary_interval(details, draws=2000, seed=20261122):
    """Paired within-group bootstrap, conditional on recorded held-out fits."""
    records = []
    for chain, frame in details.groupby("chain"):
        rng = np.random.default_rng(seed)
        mice = frame.pivot(index=["group", "mouse"], columns="method", values="mean_raw_recovery")
        delta = mice["structural_screen"] - mice["edgeR_binary_ungated"]
        bootstrap = {}
        observed = {}
        for group in ("g1", "g5", "g6"):
            values = delta.loc[group].to_numpy()
            if not np.isfinite(values).all():
                raise ValueError(f"{chain}: missing paired responses in {group}")
            sampled = rng.integers(0, len(values), size=(draws, len(values)))
            bootstrap[group] = values[sampled].mean(axis=1)
            observed[group] = values.mean()
        differences = bootstrap["g1"] - (bootstrap["g5"] + bootstrap["g6"]) / 2
        lo, hi = np.quantile(differences, [.025, .975])
        records.append(dict(chain=chain, contrast=observed["g1"] - (observed["g5"] + observed["g6"]) / 2,
                            bootstrap_lo=lo, bootstrap_hi=hi, draws=draws, seed=seed))
    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    metadata = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})
    records = []
    for chain in ("TRA", "TRB"):
        calibration = fold_path(calibration_dir(root, args.cache, chain),
                                args.cache / "calibration" / chain,
                                f"count_snapshots/{chain}_manifest.json")
        directory = evidence_dir(root, args.cache, chain, calibration=calibration)
        folds = pd.read_csv(directory / "fold_evidence.csv", dtype={"heldout_mouse": str})
        raw = load_snapshot_tables(args.cache / "raw", chain, metadata)
        records.extend(selection_rows(chain, folds, raw, metadata))
    details = pd.DataFrame(records)
    args.output.mkdir(parents=True, exist_ok=True)
    details.to_csv(args.output / "internal_1249_selection_by_mouse.csv", index=False)
    contrast_summary(details).to_csv(args.output / "internal_1249_selection_summary.csv", index=False)
    paired_primary_interval(details).to_csv(args.output / "primary_paired_bootstrap.csv", index=False)


if __name__ == "__main__":
    main()
