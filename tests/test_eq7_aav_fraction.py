import numpy as np
import pandas as pd
import pytest

from src import eq7_aav_fraction
from src.eq7_aav_fraction import FRACTION, fraction_evidence, summarize_positive_fraction
from src.spike_in import EDGER_SUCCESS


KEYS = pd.Index(
    [("C1", "V1"), ("C2", "V1"), ("C3", "V2")], tupleize_cols=False,
)


def test_fraction_is_positive_significant_over_all_tested_aa_v():
    tests = pd.DataFrame({
        "feature_id": list(KEYS),
        "FDR": [0.01, 0.01, 0.2],
        "logFC": [2.0, -3.0, 5.0],
    })
    frac, n_pos, n_tested = summarize_positive_fraction(
        tests, "FDR", "logFC", pd.Index(["V1", "V2", "V3"]),
    )
    np.testing.assert_array_equal(frac, [0.5, 0, 0])
    np.testing.assert_array_equal(n_pos, [1, 0, 0])
    np.testing.assert_array_equal(n_tested, [2, 1, 0])


def test_omitted_mouse_counts_cannot_affect_fold_fractions(monkeypatch):
    counts = pd.DataFrame({
        "g1|A": [100, 100, 500],
        "g1|B": [10, 20, 3],
        "g1|C": [15, 18, 2],
        "g5|D": [4, 9, 3],
        "g6|E": [4, 8, 2],
    }, index=KEYS)
    units = pd.DataFrame({
        "group_no": ["g1", "g1", "g1", "g5", "g6"],
        "mouse_id": ["A", "B", "C", "D", "E"],
        "contrast_group": ["g1", "g1", "g1", "allogeneic", "allogeneic"],
    }, index=counts.columns)
    seen = []

    def fake_edger(train, metadata):
        assert "g1|A" not in train.columns
        assert "g1|A" not in metadata.index
        seen.append(train.to_numpy().copy())
        return pd.DataFrame({
            "feature_id": [KEYS[0], KEYS[2]], "FDR": [0.01, 0.02],
            "logFC": [2.0, -1.0],
        }), EDGER_SUCCESS

    def fake_fisher(train, metadata, n_workers):
        assert n_workers == 2
        return pd.DataFrame({
            "feature_id": list(KEYS),
            "Fisher_FDR_g1": [0.01, 0.01, 0.5],
            "Fisher_log2FC_g1": [1.0, -2.0, 1.0],
        })

    monkeypatch.setattr(eq7_aav_fraction, "fit_edger", fake_edger)
    monkeypatch.setattr(eq7_aav_fraction, "fit_fisher", fake_fisher)
    first = fraction_evidence(counts, units, (("g1", "A"),), workers=2)
    mutated = counts.copy()
    mutated["g1|A"] = [1_000_000, 0, 0]
    second = fraction_evidence(mutated, units, (("g1", "A"),), workers=2)
    pd.testing.assert_frame_equal(first, second)
    np.testing.assert_array_equal(seen[0], seen[1])
    row = first.set_index("v_gene")
    assert row.loc["V1", FRACTION[0]] == pytest.approx(1)
    assert row.loc["V1", FRACTION[1]] == pytest.approx(0.5)
    assert row.loc["V2", FRACTION[0]] == pytest.approx(0)
    assert row.loc["V2", FRACTION[1]] == pytest.approx(0)


def test_non_significant_positive_effect_is_zero_fraction():
    tests = pd.DataFrame({
        "feature_id": [KEYS[0], KEYS[1]],
        "FDR": [0.08, 0.9],
        "logFC": [4.0, 5.0],
    })
    fractions, positive, tested = summarize_positive_fraction(
        tests, "FDR", "logFC", pd.Index(["V1"]),
    )
    assert fractions.tolist() == [0.0]
    assert positive.tolist() == [0]
    assert tested.tolist() == [2]
