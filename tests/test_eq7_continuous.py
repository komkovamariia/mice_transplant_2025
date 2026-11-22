import numpy as np
import pandas as pd
import pytest

from src import eq7_continuous
from src.eq7_continuous import CONTINUOUS, JointRecoveryPanel, continuous_evidence, prepare_v_evidence
from src.spike_in import EDGER_SUCCESS


def test_joint_recovery_changes_with_breadth_at_fixed_umi():
    one = pd.Series([20, 20], index=[("A", "V1"), ("B", "V1")])
    two = pd.Series([20, 10, 10], index=[("A", "V1"), ("B", "V1"), ("C", "V1")])
    panel = JointRecoveryPanel({("g1", "one"): one, ("g1", "two"): two})
    measured = panel.measure(np.array([1, 0, 0]), [0, 1])
    np.testing.assert_allclose(measured[:, 1, 0], np.sqrt([500_000 / 2, 500_000 / 3]))
    increased = one.reindex(panel.keys, fill_value=0).copy()
    increased.iloc[:] = [30, 10, 0]
    larger = JointRecoveryPanel({("g1", "larger"): increased}, align_to=panel)
    assert larger.measure(np.array([1, 0, 0]), [0])[0, 1, 0] > measured[0, 1, 0]


def test_quantitative_features_keep_sign_and_omit_all_mouse_libraries(monkeypatch):
    tables = {
        "g1": pd.DataFrame({"m1_a": [1000, 2000], "m1_b": [1000, 2000], "m2": [5, 8]},
                           index=[("A", "V1"), ("B", "V1")]),
        "g5": pd.DataFrame({"m3": [20, 70]}, index=[("A", "V1"), ("B", "V2")]),
    }
    metadata = pd.DataFrame({
        "sample_id": ["m1_a", "m1_b", "m2", "m3"], "group_no": ["g1", "g1", "g1", "g5"],
        "mouse_id": ["m1", "m1", "m2", "m3"], "mouse_no": [1, 1, 2, 3],
    })
    counts, units = prepare_v_evidence(tables, metadata)

    def fake_edge(train, retained):
        assert train.columns.tolist() == retained.index.tolist() == ["g1|m2", "g5|m3"]
        assert train.loc[[("V", "V1")]].iloc[0].tolist() == [13, 20]
        return pd.DataFrame({"feature_id": [("V", "V1")], "logFC": [-2], "FDR": [.2]}), EDGER_SUCCESS

    monkeypatch.setattr(eq7_continuous, "fit_edger", fake_edge)
    got = continuous_evidence(counts, units, (("g1", "m1"),)).set_index("v_gene")
    assert got.at["V1", CONTINUOUS[0]] == pytest.approx(-1.6)
    assert got.at["V1", CONTINUOUS[1]] > 0
    assert got.at["V2", CONTINUOUS[0]] == 0
    assert got.at["V2", CONTINUOUS[1]] < 0
