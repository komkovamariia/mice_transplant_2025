#!/usr/bin/env python3
"""Train and evaluate the three joint-recovery regressions by splenic T-cell subset."""

import argparse
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_continuous import JOINT_VALUE, PENALTIES, RANKERS, predictions, save
from scripts.run_eq7_external import external_raw
from scripts.run_eq7_weighted import CHAINS, ORDERS, ranking_rows
from src.differential_abundance import fit_edger, fit_fisher, prepare_counts, rank_v_segments
from src.eq7_aav_fraction import summarize_positive_fraction
from src.eq7_calibration import FEATURES, profile_coordinates, shortlist, with_structural_fields
from src.eq7_continuous import (
    BINARY_MODEL, CONTINUOUS, CONTINUOUS_MODEL, JOINT_CONTRAST, JointRecoveryPanel,
    continuous_evidence, prepare_v_evidence,
)
from src.eq7_weighted import GROUPS, contrasts, fit_ridge, pool_mice, training_contrasts
from src.spike_in import EDGER_SUCCESS, load_snapshot_tables
from src.structural_sampling import downsample_group_tables
from src.structural_v import structural_remainder, subset_v_tables
from src.v_thresholds import PlotThresholds

DEPTH = 15_000
SEED = 1031
SUBSETS = {"cd4": ("cd4",), "cd8": ("cd8",), "pooled": ("cd4", "cd8")}
FRACTION_MODEL = "aaV_fraction_ridge"
FRACTIONS = ("edgeR_aav_positive_fraction", "fisher_aav_positive_fraction", FEATURES[2])
MODELS = {BINARY_MODEL: FEATURES, CONTINUOUS_MODEL: CONTINUOUS, FRACTION_MODEL: FRACTIONS}


def matched_internal(raw, metadata, subset):
    """Restrict all three panels of a chain to the same mice and both passing phenotypes."""
    libraries = pd.concat([table.sum(axis=0) for table in raw.values()]).rename("functional_umi")
    audit = metadata.merge(libraries, left_on="sample_id", right_index=True, validate="one_to_one")
    audit = audit.loc[audit["_tissue"].eq("spleen") & audit["_cell_subset"].isin(SUBSETS["pooled"])].copy()
    passes = audit["functional_umi"].ge(DEPTH)
    mouse_key = ["group_no", "mouse_id"]
    paired_mice = audit.loc[passes].groupby(mouse_key)["_cell_subset"].nunique().eq(2)
    audit["paired_mouse"] = pd.MultiIndex.from_frame(audit[mouse_key]).isin(paired_mice[paired_mice].index)
    audit["included"] = audit["paired_mouse"] & passes & audit["_cell_subset"].isin(SUBSETS[subset])
    return subset_v_tables(raw, audit.loc[audit["included"], "sample_id"]), audit


def training_predictors(tables, counts, units, v_counts, v_units, mice):
    keys = pd.MultiIndex.from_frame(units[["group_no", "mouse_id"]])
    units = units.loc[~keys.isin(mice)]
    train = counts.loc[:, units.index]
    contrast = units.index[units["contrast_group"].isin(("g1", "allogeneic"))]
    train = train.loc[train[contrast].sum(axis=1).gt(0)]
    edge, status = fit_edger(train, units)
    if status != EDGER_SUCCESS:
        raise ValueError(status)
    fisher = fit_fisher(train, units, n_workers=1)
    remainder = structural_remainder(tables)
    indicators = rank_v_segments(remainder, edge, fisher, "spleen", "Splenic training mice")
    structure = with_structural_fields(indicators, tables, remainder=remainder)
    effects = continuous_evidence(v_counts, v_units, mice)
    rows = structure.merge(effects, on="v_gene", how="left", validate="one_to_one")
    for test, fdr, effect, name in (
        (edge, "FDR", "logFC", FRACTIONS[0]),
        (fisher, "Fisher_FDR_g1", "Fisher_log2FC_g1", FRACTIONS[1]),
    ):
        rows[name], _, _ = summarize_positive_fraction(test, fdr, effect, rows["v_gene"])
    if rows[list(CONTINUOUS[:2])].isna().any().any():
        raise ValueError("A structural candidate is missing from the raw V-segment counts.")
    return rows


def model_specs():
    for model, terms in MODELS.items():
        for penalty in PENALTIES:
            yield model, terms, penalty
        if model != BINARY_MODEL:
            for name, omitted in zip(("edgeR", "Fisher", "structure"), terms, strict=True):
                yield f"{model}_without_{name}", tuple(term for term in terms if term != omitted), 1.0


def train(chain, metadata, raw, tables, output):
    counts, units = prepare_counts(raw, metadata, pool_within_mouse=True)
    v_counts, v_units = prepare_v_evidence(raw, metadata)
    panel = JointRecoveryPanel(pool_mice(tables, metadata))
    sample_ids = metadata.groupby(["group_no", "mouse_id"])["sample_id"].agg(list)
    all_samples = set(sample_ids.explode())
    partitions = list(product(*(panel.by_group()[group] for group in GROUPS)))
    cached_measurements, fits, ranks, refits = {}, [], [], []
    for index, triplet in enumerate(partitions, start=1):
        mice = tuple(panel.mice[i] for i in triplet)
        split = "split_" + "_".join(map(str, triplet))
        omitted = {sample for key in mice for sample in sample_ids.loc[key]}
        training_tables = subset_v_tables(tables, all_samples - omitted)
        eligible = shortlist(training_predictors(training_tables, counts, units, v_counts, v_units, mice))
        if eligible.empty:
            raise ValueError(f"Empty structural shortlist for {chain}/{split}.")
        genes = panel.genes.get_indexer(eligible["v_gene"])
        measured = panel.measure(panel.reference(triplet), triplet)[:, :, genes]
        observed = contrasts(measured)[1]
        outcomes = training_contrasts(panel, triplet, cached_measurements, recovery_column=JOINT_CONTRAST)
        training = eligible.merge(outcomes[["v_gene", JOINT_CONTRAST]], on="v_gene", validate="one_to_one")
        refits.append(eligible.assign(**{JOINT_CONTRAST: observed}, fit_weight=1 / len(eligible)))
        order = {model: eligible[term].to_numpy() for model, term in RANKERS.items()}
        for model, terms, penalty in model_specs():
            fit = fit_ridge(training, terms, JOINT_CONTRAST, penalty)
            pred = predictions(eligible, fit, model, penalty, observed)
            fits.append(pred.assign(chain=chain, split=split, **{f"{g}_mouse": m for g, m in mice}))
            if penalty == 1 and model in MODELS:
                order[model] = pred["prediction"].to_numpy()
        for row in ranking_rows(eligible, measured, mice, ORDERS[chain], order, recovery_column=JOINT_VALUE):
            ranks.append({"chain": chain, "split": split, **row})
        if index == 1 or index % 5 == 0:
            print(f"{output.name}: {index}/{len(partitions)} triplets", flush=True)

    save(output, "internal", [pd.DataFrame(ranks)], fits, ["chain"])
    training = pd.concat(refits, ignore_index=True)
    frozen = {(model, penalty): fit_ridge(training, terms, JOINT_CONTRAST, penalty, weight="fit_weight")
              for model, terms, penalty in model_specs() if model in MODELS}
    coefficients = [(model, penalty, term, weight, fit.intercept, fit.baseline)
                    for (model, penalty), fit in frozen.items() for term, weight in fit.coefficients.items()]
    pd.DataFrame(coefficients, columns=["model", "penalty", "term", "weight", "intercept",
                                         "constant_baseline"]).to_csv(output / "frozen_coefficients.csv", index=False)
    rows = profile_coordinates(training_predictors(tables, counts, units, v_counts, v_units, ()))
    candidates = shortlist(rows)
    rows.assign(selected=rows["v_gene"].isin(candidates["v_gene"])).to_csv(
        output / "full_internal_structure.csv", index=False)
    candidates.to_csv(output / "frozen_candidates.csv", index=False)
    return panel, frozen, candidates


def external_test(root, chain, study, reports, expected_qc, panel, frozen, candidates, output):
    meta = pd.read_csv(root / "config" / f"{study.lower()}_samples.csv", dtype={"mouse_id": str})
    counts, qc = external_raw(reports, meta, chain)
    if expected_qc is not None:
        reference = pd.read_csv(expected_qc).set_index("sample_id")["report_sha256"]
        actual = qc.set_index("sample_id")["report_sha256"]
        if not reference.reindex(actual.index).eq(actual).all():
            raise ValueError(f"{study} TRUST4 checksums differ from the archived inputs.")
    tables, audit = downsample_group_tables(counts, depth=DEPTH, seed=SEED)
    qc.merge(audit, on="sample_id").to_csv(output / "external_library_qc.csv", index=False)
    query = JointRecoveryPanel(pool_mice(tables, meta), align_to=panel)
    genes = panel.genes.get_indexer(candidates["v_gene"])
    measures = query.measure(panel.reference(), range(len(query.mice)))[:, :, genes]
    per_mouse = [dict(group=group, mouse=mouse, v_gene=gene, count_recovery_pct=count, joint_recovery=value)
                 for (group, mouse), measure in zip(query.mice, measures, strict=True)
                 for gene, count, value in zip(candidates["v_gene"], *measure, strict=True)]
    pd.DataFrame(per_mouse).to_csv(output / "external_mouse_recovery.csv", index=False)
    if study == "PRJNA650133":
        return
    groups = query.by_group()
    observed = measures[list(groups["g1"]), 1].mean(axis=0) - measures[list(groups["g5"]), 1].mean(axis=0)
    scored = [predictions(candidates, fit, model, penalty, observed).assign(study=study, chain=chain)
              for (model, penalty), fit in frozen.items()]
    ranks = {model: fit.predict(candidates) for (model, penalty), fit in frozen.items() if penalty == 1}
    rec = list(ranking_rows(candidates, measures, query.mice, ORDERS[chain], ranks, recovery_column=JOINT_VALUE))
    save(output, "external", [pd.DataFrame(rec).assign(study=study, chain=chain)], scored, ["study", "chain"])

    # Only intact and graft recipients exist in this cohort. Apply the frozen numerical cutoffs.
    intact, graft = tables["g1"].index, tables["g5"].index
    full = pd.Series([key[1] for key in intact]).value_counts()
    residual = pd.Series([key[1] for key in intact.difference(graft)], dtype=str).value_counts()
    structure = pd.DataFrame({"n_full": full, "n_g1_exclusive_aav": residual}).fillna(0)
    structure = profile_coordinates(structure.rename_axis("v_gene").reset_index())
    cutoffs = PlotThresholds(*candidates.iloc[0][["training_share_cutoff_pct", "training_retention_cutoff_pct"]])
    structure["selected"] = cutoffs.includes(structure["initial_share_pct"], structure["retention_pct"])
    structure.to_csv(output / "external_structure.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chain", choices=CHAINS, required=True)
    parser.add_argument("--subset", choices=SUBSETS, required=True)
    parser.add_argument("--skin", type=Path)
    parser.add_argument("--marrow", type=Path)
    parser.add_argument("--article-inputs", type=Path, help="Optional directory with archived TRUST4 QC checksums.")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    meta = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"mouse_id": str})
    candidates = meta.loc[meta["chain"].eq(args.chain) & meta["_tissue"].eq("spleen"), "sample_id"]
    raw = load_snapshot_tables(args.cache / "raw", args.chain, meta, sample_ids=candidates)
    raw, audit = matched_internal(raw, meta, args.subset)
    args.output.mkdir(parents=True, exist_ok=True)
    audit.to_csv(args.output / "internal_libraries.csv", index=False)
    tables, _ = downsample_group_tables(raw, depth=DEPTH, seed=SEED)
    panel, frozen, candidates = train(args.chain, meta, raw, tables, args.output)
    for study, reports, compatible, qc_subdir in (
        ("PRJNA954716", args.skin, (args.chain, args.subset) == ("TRB", "cd4"), Path()),
        ("PRJNA650133", args.marrow, args.subset == "pooled", Path(args.chain)),
    ):
        if reports and compatible:
            expected = (args.article_inputs / study / "validation" / qc_subdir / "external_library_qc.csv"
                        if args.article_inputs else None)
            external_test(root, args.chain, study, reports, expected, panel, frozen, candidates, args.output)
    run = {"chain": args.chain, "subset": args.subset, "tissue": "spleen", "depth": DEPTH,
           "seed": SEED, "penalties": PENALTIES, "mice": {g: len(i) for g, i in panel.by_group().items()}}
    (args.output / "run.json").write_text(json.dumps(run, indent=2) + "\n")


if __name__ == "__main__":
    main()
