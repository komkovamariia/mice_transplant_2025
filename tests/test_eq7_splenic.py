import pandas as pd

from scripts.run_eq7_splenic import DEPTH, SUBSETS, matched_internal, snapshot_totals
from src.spike_in import load_snapshot_tables, snapshot_tables


def test_paired_mouse_filter_is_identical_across_splenic_subsets():
    samples = pd.DataFrame([
        ("a4", "g1", "a", "spleen", "cd4", DEPTH),
        ("a8", "g1", "a", "spleen", "cd8", DEPTH),
        ("b4", "g1", "b", "spleen", "cd4", DEPTH),
        ("b8", "g1", "b", "spleen", "cd8", DEPTH - 1),
        ("bt8", "g1", "b", "thymus", "cd8", DEPTH),
    ], columns=["sample_id", "group_no", "mouse_id", "_tissue", "_cell_subset", "count"])
    keys = pd.Index([("CASSF", "TRBV2")], tupleize_cols=False)
    raw = {"g1": pd.DataFrame([samples["count"].tolist()], index=keys, columns=samples["sample_id"])}
    for subset, expected in zip(SUBSETS, (["a4"], ["a8"], ["a4", "a8"]), strict=True):
        tables, audit = matched_internal(raw, samples, subset)
        assert tables["g1"].columns.tolist() == expected
        assert audit.loc[audit["included"], "mouse_id"].unique().tolist() == ["a"]


def test_selective_snapshot_read_matches_the_full_loaded_slice(tmp_path):
    samples = pd.DataFrame([
        ("a4", "g1", "a", "spleen", "cd4"),
        ("a8", "g1", "a", "spleen", "cd8"),
    ], columns=["sample_id", "group_no", "mouse_id", "source", "subtype"])
    keys = pd.Index([("CASSF", "TRBV2"), ("CASRF", "TRBV19")], tupleize_cols=False)
    full = pd.DataFrame({"a4": [1, 0], "a8": [0, 2]}, index=keys)
    snapshot_tables({"g1": full}, tmp_path, "TRB", samples)
    totals = snapshot_totals(tmp_path, "TRB", samples["sample_id"])
    pd.testing.assert_series_equal(totals, full.sum(axis=0))
    selected = load_snapshot_tables(tmp_path, "TRB", samples, sample_ids=["a8"])
    pd.testing.assert_frame_equal(selected["g1"], full.loc[full["a8"].gt(0), ["a8"]])


def test_combined_panel_requires_both_phenotypes_in_both_tissues():
    samples = pd.DataFrame([
        (f"{mouse}{tissue}{subset}", "g1", mouse, tissue, subset,
         DEPTH - 1 if (mouse, tissue, subset) == ("b", "thymus", "cd8") else DEPTH)
        for mouse in ("a", "b") for tissue in ("spleen", "thymus") for subset in ("cd4", "cd8")
    ], columns=["sample_id", "group_no", "mouse_id", "_tissue", "_cell_subset", "count"])
    keys = pd.Index([("CASSF", "TRBV2")], tupleize_cols=False)
    raw = {"g1": pd.DataFrame([samples["count"].tolist()], index=keys, columns=samples["sample_id"])}
    tables, audited = matched_internal(raw, samples, "cd4", "combined")
    assert list(tables["g1"].columns) == ["aspleencd4", "athymuscd4"]
    assert audited.loc[audited["included"], "mouse_id"].unique().tolist() == ["a"]


def test_vectorized_fisher_matches_scipy_on_discrete_pooled_counts():
    import numpy as np
    from scipy.stats import fisher_exact

    from src.count_statistics import pooled_fisher_arrays

    a = np.array([0, 1, 5, 12, 450, 2970, 9, 9])
    b = np.array([0, 2, 9, 12, 1767, 4800, 0, 0])
    total_a, total_b = a.sum(), b.sum()
    reference = np.array([fisher_exact([[x, y], [total_a-x, total_b-y]]) for x, y in zip(a, b)])
    actual = np.column_stack(pooled_fisher_arrays(a, b))
    np.testing.assert_allclose(actual, reference, atol=1e-12, rtol=1e-10)
