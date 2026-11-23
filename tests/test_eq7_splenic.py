import pandas as pd

from scripts.run_eq7_splenic import DEPTH, SUBSETS, matched_internal
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
    selected = load_snapshot_tables(tmp_path, "TRB", samples, sample_ids=["a8"])
    pd.testing.assert_frame_equal(selected["g1"], full.loc[full["a8"].gt(0), ["a8"]])
