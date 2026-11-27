#!/usr/bin/env python3
"""Update cached Fisher indicators without changing mouse splits or edgeR evidence."""

import argparse
from pathlib import Path
import sys

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_v_updates import calibration_dir, evidence_dir
from src.differential_abundance import fit_fisher, prepare_counts
from src.spike_in import aa_v_key, load_snapshot_tables


def fisher_indicators(counts, units, omitted, workers):
    keys = pd.MultiIndex.from_frame(units[["group_no", "mouse_id"]])
    retained = units.loc[~keys.isin(omitted)]
    contrast = retained.index[retained["contrast_group"].isin(("g1", "allogeneic"))]
    train = counts.loc[:, retained.index]
    train = train.loc[train[contrast].sum(axis=1).gt(0)]
    result = fit_fisher(train, retained, n_workers=workers)
    significant = result.loc[result["Fisher_FDR_g1"].lt(.05)].copy()
    significant["v_gene"] = significant["feature_id"].map(lambda key: aa_v_key(key)[1])
    return significant.groupby("v_gene")["Fisher_log2OR_g1"].mean().gt(0)


def refresh_chain(root, cache, chain, metadata, workers, source):
    target = evidence_dir(root, cache, chain, calibration=source[0])
    old = source[1]
    target.mkdir(parents=True, exist_ok=True)
    if target == old:
        return target
    tables = load_snapshot_tables(cache / "raw", chain, metadata)
    counts, units = prepare_counts(tables, metadata, pool_within_mouse=True)
    for filename in ("fold_evidence.csv", "inner_fold_evidence.csv", "stratified_fold_evidence.csv"):
        path = old / filename
        if not path.is_file():
            raise FileNotFoundError(path)
        rows = pd.read_csv(path, dtype={"heldout_mouse": str})
        for split, split_rows in rows.groupby("evidence_split", sort=False):
            pending = target / f"{split}.csv"
            if not pending.is_file():
                omitted = tuple(split_rows[["heldout_group", "heldout_mouse"]]
                                .drop_duplicates().itertuples(index=False, name=None))
                indicators = fisher_indicators(counts, units, omitted, workers)
                updated = split_rows.assign(
                    fisher_g1_enriched=split_rows["v_gene"].map(indicators).eq(True).astype(int)
                )
                temp = pending.with_suffix(".tmp")
                updated.to_csv(temp, index=False)
                temp.replace(pending)
                print(f"{chain}: {split}", flush=True)
        result = pd.concat([pd.read_csv(target / f"{split}.csv", dtype={"heldout_mouse": str})
                            for split in rows["evidence_split"].unique()], ignore_index=True)
        result.to_csv(target / filename, index=False)
    for name in ("fold_status.csv", "stratified_fold_status.csv"):
        (target / name).write_bytes((old / name).read_bytes())
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--chain", choices=("TRA", "TRB"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    metadata = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"sample_id": str, "mouse_id": str})
    for chain in (args.chain,) if args.chain else ("TRA", "TRB"):
        calibration = calibration_dir(root, args.cache, chain)
        if not (calibration / "count_snapshots").is_dir():
            matches = list((args.cache / "calibration" / chain).glob("*/count_snapshots"))
            if len(matches) != 1:
                raise ValueError(f"Cannot identify unique calibration snapshots for {chain}")
            calibration = matches[0].parent
        target = evidence_dir(root, args.cache, chain, calibration=calibration)
        matches = [folder.parent for folder in (args.cache / "eq7_evidence" / chain).glob("*/stratified_fold_evidence.csv")
                   if folder.parent != target]
        if (target / "stratified_fold_evidence.csv").is_file():
            print(f"{chain}: odds-ratio cache already available at {target}")
            continue
        if len(matches) != 1:
            raise ValueError(f"Cannot identify unique source evidence for {chain}: {matches}")
        print(refresh_chain(root, args.cache, chain, metadata, args.workers, (calibration, matches[0])))


if __name__ == "__main__":
    main()
