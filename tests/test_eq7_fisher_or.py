import numpy as np
import pandas as pd
import pytest

from src import eq7_fisher_or
from src.eq7_fisher_or import (
    ODDS_FISHER, RAW_FISHER, add_log2_odds_ratio,
    fisher_fraction_sensitivity,
)


def test_odds_direction_adjusts_for_group_depth():
    # The raw count ratio is negative (10 < 20), while the within-group
    # proportions are 10/100 and 20/1000 and the odds ratio is positive.
    fisher = pd.DataFrame({
        "feature_id": [("C1", "V1"), ("C2", "V1")],
        "fisher_g1_count": [10, 90],
        "fisher_allogeneic_count": [20, 980],
        "Fisher_FDR_g1": [0.01, 0.8],
        "Fisher_log2FC_g1": [np.log2(11/21), np.log2(91/981)],
    })
    result = add_log2_odds_ratio(fisher)
    assert result["Fisher_log2FC_g1"].iat[0] < 0
    assert result["Fisher_log2OR_g1"].iat[0] > 0
    assert result["Fisher_log2OR_g1"].iat[1] < 0


def test_fisher_fraction_changes_sign_without_refitting_on_test_mouse(monkeypatch):
    keys = pd.Index([("C1", "V1"), ("C2", "V1")], tupleize_cols=False)
    counts = pd.DataFrame({
        "g1|A": [100, 100],
        "g1|B": [5, 50],
        "g1|C": [5, 40],
        "g5|D": [20, 980],
    }, index=keys)
    units = pd.DataFrame({
        "group_no": ["g1", "g1", "g1", "g5"],
        "mouse_id": ["A", "B", "C", "D"],
        "contrast_group": ["g1", "g1", "g1", "allogeneic"],
    }, index=counts.columns)
    seen = []

    def fake_fisher(training, retained, n_workers):
        assert "g1|A" not in training.columns
        assert "g1|A" not in retained.index
        seen.append(training.copy())
        # Test statistics on the retained mice, same for both runs.
        return pd.DataFrame({
            "feature_id": list(keys),
            "fisher_g1_count": [10, 90],
            "fisher_allogeneic_count": [20, 980],
            "Fisher_log2FC_g1": [np.log2(11/21), np.log2(91/981)],
            "Fisher_FDR_g1": [0.01, 0.8],
        })

    monkeypatch.setattr(eq7_fisher_or, "fit_fisher", fake_fisher)
    first = fisher_fraction_sensitivity(counts, units, (("g1", "A"),), workers=4)
    changed = counts.copy()
    changed["g1|A"] = [2_000_000, 0]
    second = fisher_fraction_sensitivity(changed, units, (("g1", "A"),), workers=4)
    pd.testing.assert_frame_equal(first, second)
    pd.testing.assert_frame_equal(seen[0], seen[1])
    assert first[RAW_FISHER].iat[0] == pytest.approx(0)
    assert first[ODDS_FISHER].iat[0] == pytest.approx(0.5)
    assert first["n_Fisher_tested"].iat[0] == 2


def test_zero_pool_is_rejected():
    fisher = pd.DataFrame({
        "feature_id": [("C1", "V1")],
        "fisher_g1_count": [0],
        "fisher_allogeneic_count": [4],
    })
    with pytest.raises(ValueError, match="Both groups"):
        add_log2_odds_ratio(fisher)
