"""Fit the notebook's edgeR model on synthetic mouse-level libraries."""

import importlib.util
import shutil

import numpy as np
import pandas as pd
import pytest


@pytest.mark.skipif(
    shutil.which("Rscript") is None or importlib.util.find_spec("rpy2") is None,
    reason="Requires R, edgeR and rpy2; exercised by the edgeR integration CI job.",
)
def test_notebook_edger_recovers_injected_aav_with_real_r():
    from src.spike_in import apply_spike_to_tables, combine_count_tables

    rng = np.random.default_rng(1031)
    # independent mouse columns; shared background supports dispersion estimation
    keys = pd.Index([(f"CASS{i:04d}", "TRAV1") for i in range(200)], tupleize_cols=False)
    # identifiers here label synthetic features; they are not asserted to be biological sequences
    tables = {}
    for group in ("g1", "g2", "g4", "g5", "g6"):
        columns = [f"{group}_m{i}" for i in range(6)]
        tables[group] = pd.DataFrame(rng.negative_binomial(20, 0.2, (200, 6)),
                                     index=keys, columns=columns)
    selection = {"v_gene": "TRAV9", "clonotypes": [["CASSQ", "TRAV9"], ["CASSR", "TRAV9"]]}
    changed, _ = apply_spike_to_tables(tables, selection, .05)
    counts = combine_count_tables({g: changed[g] for g in ("g1", "g5", "g6")})
    metadata = pd.DataFrame({"contrast_group": [
        "g1" if sample.startswith("g1_") else "allogeneic" for sample in counts.columns
    ]}, index=counts.columns)
    from src.differential_abundance import fit_edger

    result, status = fit_edger(counts, metadata)
    assert status == "edgeR quasi-likelihood model fitted successfully."
    injected = result[result["feature_id"].isin([("CASSQ", "TRAV9"), ("CASSR", "TRAV9")])]
    assert len(injected) == 2
    assert injected["logFC"].gt(0).all()
    assert injected["FDR"].lt(.05).all()
