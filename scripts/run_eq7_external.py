#!/usr/bin/env python3
"""Prepare external repertoires or test frozen Equation 7 weights."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from http.client import HTTPException
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen

import numpy as np
import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_analysis import repo_root
from scripts.run_v_updates import _metadata, _write_csv, calibration_input_sha256, input_dir
from src.eq7_calibration import FEATURES, CalibrationUnavailable, score_candidates
from src.eq7_external import (
    MATCHED_GROUPS, compare_shortlists, matched_panel, mouse_prediction_losses, mouse_predictions,
    panel_cutoffs, summarize_transfer,
)
from src.runtime import configure_runtime
from src.spike_in import load_snapshot_tables
from src.structural_sampling import downsample_group_tables


STUDY = "PRJNA954716"
ZHENG = "PRJNA650133"
ASSEMBLER = "TRUST4 v1.1.7"
CHAIN = "TRB"
ENA = "https://ftp.sra.ebi.ac.uk/vol1/fastq"
PRODUCTIVE_CDR3 = r"[ACDEFGHIKLMNPQRSTVWY]+"
V_GENES = {
    "TRA": re.compile(r"(?:TRA|TRD)V\d+[DN]?(?:-\d+)?(?:[/\-]DV\d+(?:-\d+)?)?"),
    "TRB": re.compile(r"TRBV\d+(?:-\d+)?"),
}


def checksum(path, algorithm="sha256"):
    with path.open("rb") as source:
        return hashlib.file_digest(source, algorithm).hexdigest()


def fetch_read(reads, run, mate, expected):
    path = reads / f"{run}_{mate}.fastq.gz"
    if path.is_file() and checksum(path, "md5") == expected:
        return path
    url = f"{ENA}/{run[:6]}/0{run[-2:]}/{run}/{path.name}"
    pending = path.with_suffix(".tmp")
    for attempt in range(3):
        try:
            with urlopen(url, timeout=90) as response, pending.open("wb") as target:
                shutil.copyfileobj(response, target)
            if checksum(pending, "md5") != expected:
                raise ValueError(f"ENA checksum mismatch for {path.name}.")
        except (OSError, HTTPException, ValueError):
            pending.unlink(missing_ok=True)
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)
        else:
            pending.replace(path)
            return path


def prepare_public(work, metadata, trust4, reference, threads, chains=(CHAIN,)):
    reads, reports = work / "reads", work / "trust4_reports"
    reads.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    for sample in metadata.itertuples(index=False):
        run = sample.sample_id
        report_path = reports / f"{run}_report.tsv"
        fastqs = [
            fetch_read(reads, run, mate, expected)
            for mate, expected in enumerate((sample.fastq_md5_r1, sample.fastq_md5_r2), start=1)
        ]
        provenance = {
            "assembler": ASSEMBLER, "reference_sha256": checksum(reference),
            "run_trust4_sha256": checksum(trust4),
            "trust4_binary_sha256": checksum(trust4.with_name("trust4")),
            "fastq_md5": [sample.fastq_md5_r1, sample.fastq_md5_r2],
        }
        provenance_path = reports / f"{run}_assembly.json"
        existing = json.loads(provenance_path.read_text()) if provenance_path.is_file() else {}
        if report_path.is_file() and existing == {**provenance, "report_sha256": checksum(report_path)}:
            for chain in chains:
                read_trust4_report(report_path, run, chain)
            continue
        print(f"Assembling {run} with {ASSEMBLER}", flush=True)
        with TemporaryDirectory(prefix=f"{run}_", dir=reports) as temporary:
            command = [
                str(trust4), "-f", str(reference), "--ref", str(reference),
                "-1", str(fastqs[0]), "-2", str(fastqs[1]), "--repseq", "-o", run,
                "--od", temporary, "--clean", "1", "-t", str(threads),
            ]
            with (reports / f"{run}.log").open("w") as log:
                subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
            if ASSEMBLER not in (reports / f"{run}.log").read_text():
                raise ValueError(f"Expected {ASSEMBLER} in the assembly log for {run}.")
            assembly = Path(temporary) / report_path.name
            for chain in chains:
                read_trust4_report(assembly, run, chain)
            assembly.replace(report_path)
        _write_json({**provenance, "report_sha256": checksum(report_path)}, provenance_path)


def unambiguous_v(call, chain=CHAIN):
    genes = {hit.split("*", 1)[0].strip().replace("/", "-") for hit in re.split(r"[|,;]", call)}
    if len(genes) != 1:
        return ""
    gene = genes.pop()
    return gene if V_GENES[chain].fullmatch(gene) else ""


def read_trust4_report(path, run, chain=CHAIN):
    report = pd.read_csv(path, sep="\t", usecols=["#count", "CDR3aa", "V", "J"], keep_default_na=False)
    numeric = pd.to_numeric(report["#count"], errors="raise")
    if not (numeric.gt(0) & numeric.mod(1).eq(0)).all():
        raise ValueError(f"Expected positive integer read counts for {run}.")
    report["v_gene"] = report["V"].map(lambda call: unambiguous_v(call, chain))
    functional = (
        report["CDR3aa"].str.fullmatch(PRODUCTIVE_CDR3) & report["v_gene"].ne("")
        & report["J"].str.startswith(f"{chain}J")
    )
    counts = report.loc[functional].groupby(["CDR3aa", "v_gene"])["#count"].sum().astype("int64")
    if counts.empty:
        raise ValueError(f"No productive {chain} aaV with a unique V call in {run}.")
    counts.index = pd.Index(counts.index.to_list(), tupleize_cols=False)
    qc = {
        "sample_id": run, "n_trust4_rows": len(report), f"n_productive_{chain.lower()}_aav": len(counts),
        f"n_productive_{chain.lower()}_read_fragments": int(counts.sum()), "report_sha256": checksum(path),
    }
    return counts.rename(run), qc


def external_raw(work, metadata, chain=CHAIN):
    raw, audit = {}, []
    for group, rows in metadata.groupby("group_no"):
        libraries = []
        for run in rows["sample_id"]:
            path = work / "trust4_reports" / f"{run}_report.tsv"
            counts, qc = read_trust4_report(path, run, chain)
            libraries.append(counts)
            audit.append(qc)
        raw[group] = pd.concat(libraries, axis=1).fillna(0).astype("int64").sort_index()
    return raw, pd.DataFrame(audit)


def _write_json(record, path):
    pending = path.with_suffix(".json.tmp")
    pending.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    pending.replace(path)


def frozen_weights(cache, chain):
    coefficients = cache / "eq7_weights" / f"{chain.lower()}_coefficients.csv"
    if not coefficients.is_file():
        raise FileNotFoundError("Run scripts/run_v_updates.py all to fit the internal weights before testing.")
    weights = pd.read_csv(coefficients)
    if set(weights["term"]) != set(FEATURES) or len(weights) != len(FEATURES):
        raise CalibrationUnavailable(f"A complete frozen {chain} fit is required for the external test.")
    if (not np.isfinite(weights[["weight", "intercept", "baseline_prediction"]]).all().all()
            or weights["weight"].lt(0).any()):
        raise CalibrationUnavailable(f"The frozen {chain} fit has invalid coefficients.")
    if ("calibration_input_sha256" not in weights or
            not weights["calibration_input_sha256"].eq(calibration_input_sha256(cache, chain)).all()):
        raise CalibrationUnavailable(f"The frozen {chain} weights do not match the current internal input snapshots.")
    return weights, checksum(coefficients)


def test_public(root, work, metadata):
    from src.article_figures import ARTICLE_FONTS
    from src.eq7_figures import plot_external_comparison, plot_external_predictions
    from src.figures import apply_article_style
    from src.publication_figures import save_publication_figure
    import matplotlib.pyplot as plt

    cache = input_dir(root, "all_combined")
    weights, weights_hash = frozen_weights(cache, CHAIN)
    seeds = json.loads((cache / "v_input_seeds.json").read_text())
    training = _metadata(cache)
    internal_raw = load_snapshot_tables(cache / "raw", CHAIN, training)
    matched = training.loc[
        training["_tissue"].eq("spleen") & training["_cell_subset"].eq("cd4")
        & training["group_no"].isin(MATCHED_GROUPS)
    ]
    settings = {"depth": seeds["depth"], "seed": seeds["calibration_seed"], "n_workers": configure_runtime()}
    internal, reference, train_depth = matched_panel(internal_raw, matched, **settings)
    cutoffs = panel_cutoffs(internal)
    if cutoffs is None:
        raise CalibrationUnavailable("The matched internal reference did not yield valid cutoffs.")
    train_scores = score_candidates(internal, weights, frozen_cutoffs=cutoffs)
    external_counts, qc = external_raw(work, metadata)
    external, test_tables, test_depth = matched_panel(external_counts, metadata, **settings)
    test_scores = score_candidates(external, weights, frozen_cutoffs=cutoffs)
    predictions = mouse_predictions(reference, train_scores, test_tables, metadata, weights["baseline_prediction"].iat[0])
    losses = mouse_prediction_losses(predictions, metadata)
    compared = compare_shortlists(train_scores, test_scores, internal, external)

    output = work / "validation"
    output.mkdir(parents=True, exist_ok=True)
    for name, frame in (
        ("train_reference_scores", train_scores), ("external_scores", test_scores),
        ("external_mouse_predictions", predictions), ("external_mouse_losses", losses),
        ("external_segment_comparison", compared), ("external_library_qc", metadata.merge(qc, on="sample_id")),
        ("train_reference_depth", train_depth),
        ("external_depth", test_depth.rename(columns={
            "functional_umi_before": "read_fragments_before", "target_umi": "target_read_fragments",
            "functional_umi_after": "read_fragments_after",
        })),
    ):
        _write_csv(frame, output / f"{name}.csv")
    _write_csv(weights, output / "frozen_trb_coefficients.csv")
    test_summary = {
        "study": STUDY, "chain": CHAIN, "tra_status": "not_sequenced_in_external_cohort",
        "external_conditions": {"g1": "intact BALB/c", "g5": "C57BL/6 skin graft in BALB/c"},
        "panel": "splenic CD4+ T cells", "external_counts": "read_fragments", "internal_counts": "UMI",
        "depth": seeds["depth"], "rarefaction_seed": seeds["calibration_seed"],
        "frozen_gate": cutoffs._asdict(), "training_coefficients_sha256": weights_hash,
        **summarize_transfer(losses, compared),
    }
    _write_json(test_summary, output / "external_summary.json")
    apply_article_style()
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ARTICLE_FONTS, "axes.grid": False})
    for name, figure in (
        ("external_mouse_predictions", plot_external_predictions(predictions, losses, test_summary)),
        ("external_segment_comparison", plot_external_comparison(compared, test_summary)),
    ):
        save_publication_figure(figure, output / f"{name}.png")
        plt.close(figure)
    print(f"External validation outputs written to {output}")


def test_zheng(root, work, metadata):
    """Evaluate negative-arm recovery with the locked internal g1 versus g5 reference."""
    cache = input_dir(root, "all_combined")
    seeds = json.loads((cache / "v_input_seeds.json").read_text())
    training = _metadata(cache)
    matched = training.loc[
        training["_tissue"].eq("spleen") & training["_cell_subset"].isin(("cd4", "cd8"))
        & training["group_no"].isin(MATCHED_GROUPS)
    ]
    sampling = {"depth": seeds["depth"], "seed": seeds["calibration_seed"]}
    workers = configure_runtime()
    summaries = {}
    for chain in V_GENES:
        weights, weights_hash = frozen_weights(cache, chain)
        internal_raw = load_snapshot_tables(cache / "raw", chain, training)
        evidence, reference, internal_depth = matched_panel(internal_raw, matched, n_workers=workers, **sampling)
        cutoffs = panel_cutoffs(evidence)
        if cutoffs is None:
            raise CalibrationUnavailable(f"The matched internal {chain} reference did not yield valid cutoffs.")
        candidates = score_candidates(evidence, weights, frozen_cutoffs=cutoffs)
        counts, qc = external_raw(work, metadata, chain)
        if any(table.sum(axis=0).ge(sampling["depth"]).any() for table in counts.values()):
            sampled, depth = downsample_group_tables(counts, **sampling)
        else:
            sampled = {"g5": counts["g5"].iloc[:0, :0]}
            read_totals = counts["g5"].sum(axis=0)
            depth = pd.DataFrame({
                "sample_id": read_totals.index, "group_no": "g5", "functional_umi_before": read_totals.values,
                "target_umi": sampling["depth"], "kept_for_structural_analysis": False, "functional_umi_after": 0,
            })
        depth = depth.rename(columns={
            "functional_umi_before": "read_fragments_before", "target_umi": "target_read_fragments",
            "functional_umi_after": "read_fragments_after",
        })
        predictions = mouse_predictions(reference, candidates, sampled, metadata, weights["baseline_prediction"].iat[0])
        losses = mouse_prediction_losses(predictions, metadata)
        losses["n_frozen_candidates"] = len(candidates)
        output = work / "validation" / chain
        output.mkdir(parents=True, exist_ok=True)
        for name, frame in (
            ("frozen_coefficients", weights), ("train_reference_scores", candidates), ("train_reference_depth", internal_depth),
            ("external_library_qc", metadata.merge(qc, on="sample_id")), ("external_depth", depth),
            ("external_mouse_predictions", predictions), ("external_mouse_losses", losses),
        ):
            _write_csv(frame, output / f"{name}.csv")
        summaries[chain] = {
            "n_frozen_candidates": len(candidates),
            "n_mice_depth_eligible": int(depth["kept_for_structural_analysis"].sum()),
            "n_mice_scored": int(losses["mae"].notna().sum()),
            "g5_mouse_mae": float(losses["mae"].mean()) if losses["mae"].notna().any() else None,
            "g5_baseline_mae": float(losses["baseline_mae"].mean()) if losses["baseline_mae"].notna().any() else None,
            "frozen_gate": cutoffs._asdict(), "training_coefficients_sha256": weights_hash,
        }
    _write_json({
        "study": ZHENG, "endpoint": "single_arm_g5_signed_recovery",
        "external_condition": "C57BL/6 marrow allograft in BALB/c, splenocytes at day 9",
        "external_counts": "read_fragments", "internal_counts": "UMI", "depth": sampling["depth"],
        "rarefaction_seed": sampling["seed"], "chains": summaries,
    }, work / "validation" / "external_summary.json")
    print(f"External validation outputs written to {work / 'validation'}")


def main():
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "test"))
    parser.add_argument("--study", choices=(ZHENG, STUDY), default=ZHENG)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--trust4", type=Path, help="Compiled TRUST4 v1.1.7 run-trust4 path.")
    parser.add_argument("--reference", type=Path, help="Mouse IMGT reference bundled with TRUST4.")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    work = (args.work_dir or root / "results" / "eq7_external" / args.study).resolve()
    metadata = pd.read_csv(root / "config" / f"{args.study.lower()}_samples.csv")
    if args.stage == "prepare":
        if args.trust4 is None or args.threads < 1:
            parser.error("prepare requires --trust4 and a positive --threads value")
        trust4 = args.trust4.resolve()
        reference = args.reference or trust4.parent / "mouse" / "mouse_IMGT+C.fa"
        chains = tuple(V_GENES) if args.study == ZHENG else (CHAIN,)
        prepare_public(work, metadata, trust4, reference.resolve(), args.threads, chains)
    elif args.study == ZHENG:
        test_zheng(root, work, metadata)
    else:
        test_public(root, work, metadata)


if __name__ == "__main__":
    main()
