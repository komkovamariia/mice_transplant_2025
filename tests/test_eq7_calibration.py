import pandas as pd
import pytest

from src.eq7_calibration import FEATURES, calibrate_weights, score_candidates, with_structural_fields


def test_calibration_weights_groups_equally_and_reports_each_term():
    rows = []
    for group, n_mice, sign in (("g1", 2, 1), ("g5", 3, -1), ("g6", 4, -1)):
        for i in range(n_mice):
            for j in range(3):
                rows.append({
                    "heldout_group": group, "heldout_mouse": f"mouse_{i}", "v_gene": f"TRAV{j}",
                    "n_full": 20, "only_g1_share": 0.9 - j / 10,
                    FEATURES[0]: int(group == "g1"), FEATURES[1]: 0, FEATURES[2]: j,
                    "response": sign * (0.6 - j / 10),
                })
    data = pd.DataFrame(rows)
    fit, selected = calibrate_weights(data, top=2, ridge=0.1)
    group_weights = selected.groupby("heldout_group")["regression_weight"].sum()
    assert group_weights.to_list() == pytest.approx([1 / 3] * 3)
    assert len(selected) == 18
    assert fit.loc[0, "weight"] > 0
    assert not fit.loc[1, "varied_in_training"]
    assert fit.loc[1, "weight"] == 0
    assert fit["weight"].ge(0).all()

    scored = score_candidates(data.iloc[:3], fit, top=2)
    assert len(scored) == 2
    expected = scored[[f"contribution_{term}" for term in FEATURES]].sum(axis=1)
    pd.testing.assert_series_equal(scored["calibrated_score"], expected, check_names=False)


def test_only_g1_shortlist_and_equation_7_remainder_use_distinct_backgrounds():
    keys = [("A", "TRBV1"), ("B", "TRBV1"), ("C", "TRBV2")]
    index = pd.Index(keys, tupleize_cols=False)
    g1 = pd.DataFrame({"g1_1": [8, 2, 4]}, index=index)
    tables = {"g1": g1, "g2": g1.iloc[[0]], "g4": g1.iloc[:0], "g5": g1.iloc[[1]], "g6": g1.iloc[:0]}
    verdicts = pd.DataFrame({"v_gene": ["TRBV1"], FEATURES[0]: [True], FEATURES[1]: [False]})
    result = with_structural_fields(verdicts, tables).set_index("v_gene")
    assert result.loc["TRBV1", "only_g1_share"] == pytest.approx(0.8)
    assert result.loc["TRBV1", "n_g1_exclusive_aav"] == 0
    assert result.loc["TRBV2", "n_g1_exclusive_aav"] == 1
    assert result.loc["TRBV2", FEATURES[0]] == 0
