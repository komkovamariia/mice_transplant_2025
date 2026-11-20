import pandas as pd
import pytest

from src.eq7_calibration import (
    FEATURES, STANDARDIZED_FEATURES, calibrate_weights, score_candidates, select_calibration_rows, with_structural_fields,
)


def test_calibration_weights_groups_equally_and_reports_each_term():
    rows = []
    for group, n_mice, sign in (("g1", 2, 1), ("g5", 3, -1), ("g6", 4, -1)):
        for i in range(n_mice):
            for j, (full, remaining) in enumerate([
                (1, 1), (1, 1), (10, 1), (10, 2), (10, 8), (10, 9), (10, 8), (10, 9),
            ]):
                rows.append({
                    "heldout_group": group, "heldout_mouse": f"mouse_{i}", "v_gene": f"TRAV{j}",
                    "n_full": full, "n_g1_exclusive_aav": remaining, "only_g1_share": 0.9 - j / 20,
                    FEATURES[0]: int(group == "g1"), FEATURES[1]: 0, FEATURES[2]: j,
                    "response": sign * (0.6 - j / 10),
                })
    data = pd.DataFrame(rows)
    fit, selected = calibrate_weights(data, top=2, ridge=0.1)
    group_weights = selected.groupby("heldout_group")["regression_weight"].sum()
    assert group_weights.to_list() == pytest.approx([1 / 3] * 3)
    assert len(selected) == 18
    assert set(selected["v_gene"]) == {"TRAV5", "TRAV7"}
    assert fit.loc[0, "weight"] > 0
    assert not fit.loc[1, "varied_in_training"]
    assert fit.loc[1, "weight"] == 0
    assert fit["weight"].ge(0).all()
    assert selected[list(STANDARDIZED_FEATURES)].to_numpy() @ fit["standardized_weight"].to_numpy() == (
        pytest.approx(selected["fitted_response"].to_numpy()
                      - (selected["response"] * selected["regression_weight"]).sum())
    )

    scored = score_candidates(data.iloc[:8], fit)
    assert len(scored) == 2
    expected = scored[[f"contribution_{term}" for term in FEATURES]].sum(axis=1)
    pd.testing.assert_series_equal(scored["calibrated_score"], expected, check_names=False)
    assert scored["E"].to_list() == pytest.approx(expected.to_list())
    assert scored["predicted_response"].to_list() == pytest.approx((scored["E"] + fit["intercept"].iat[0]).to_list())


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


def test_held_mouse_changes_only_evaluation_outcomes_in_its_outer_fold(tmp_path, monkeypatch):
    from scripts.run_v_updates import build_chain_evidence, fit_outer_weights
    from src.spike_in import EDGER_SUCCESS

    sizes = [(1, 1), (1, 1), (10, 1), (10, 2), (10, 8), (10, 9), (10, 8), (10, 9)]
    genes = [f"TRAV{j}" for j in range(len(sizes))]
    full = pd.Index([(str(k), gene) for gene, (n, _) in zip(genes, sizes, strict=True) for k in range(n)],
                    tupleize_cols=False)
    overlap = pd.Index([(str(k), gene) for gene, (n, retained) in zip(genes, sizes, strict=True)
                        for k in range(n - retained)], tupleize_cols=False)
    tables, metadata = {}, []
    for group, n_mice in (("g1", 4), ("g2", 1), ("g4", 1), ("g5", 2), ("g6", 2)):
        ids = []
        for i in range(n_mice):
            for library in ("cd4", "cd8"):
                sample = f"{group}_{i}_{library}"
                ids.append(sample)
                metadata.append({"sample_id": sample, "group_no": group, "mouse_id": str(i), "mouse_no": i})
        tables[group] = pd.DataFrame(1, index=full if group == "g1" else overlap, columns=ids)
    metadata = pd.DataFrame(metadata)

    def edge_r(counts, units):
        g1 = units.index[units["group_no"].eq("g1")]
        return counts.loc[:, g1].to_numpy().sum(), EDGER_SUCCESS

    def rank(_remaining, edger, _fisher, _stratum, _label):
        return pd.DataFrame({"v_gene": genes, FEATURES[0]: edger > 300, FEATURES[1]: False})

    monkeypatch.setattr("src.eq7_calibration.fit_edger", edge_r)
    monkeypatch.setattr("src.eq7_calibration.fit_fisher", lambda *_: None)
    monkeypatch.setattr("src.eq7_calibration.rank_v_segments", rank)

    changed = {**tables, "g1": tables["g1"].copy()}
    candidate = changed["g1"].index.map(lambda key: key[1] == "TRAV5" and key[0] != "0")
    changed["g1"].loc[candidate, ["g1_0_cd4", "g1_0_cd8"]] = 0

    reports = []
    for name, counts in (("original", tables), ("changed", changed)):
        directory = tmp_path / name
        directory.mkdir()
        build_chain_evidence(directory, counts, counts, metadata, workers=1)
        triples = pd.read_csv(directory / "stratified_fold_status.csv")
        assert len(triples) == 4 * 2 * 2
        assert triples[["omitted_g1", "omitted_g5", "omitted_g6"]].drop_duplicates().shape[0] == len(triples)
        fit_outer_weights(directory, directory, "tra", top=1, ridge=0.1)
        reports.append(tuple(
            pd.read_csv(directory / f"tra_outer_{report}.csv").query("outer_fold == 0").reset_index(drop=True)
            for report in ("coefficients", "scores")
        ))

    (before_weights, before_scores), (after_weights, after_scores) = reports
    assert len(before_weights) == 3
    assert len(before_scores) == 2
    pd.testing.assert_frame_equal(before_weights, after_weights)
    outcomes = ["heldout_n_v", "heldout_n_recovered", "response"]
    pd.testing.assert_frame_equal(before_scores.drop(columns=outcomes), after_scores.drop(columns=outcomes))
    v = before_scores["v_gene"].eq("TRAV5")
    assert before_scores.loc[v, "response"].item() > after_scores.loc[v, "response"].item()


def test_distinct_fold_cutoffs_and_equal_weight_for_repeated_mouse():
    rows = []
    for group in ("g1", "g5", "g6"):
        for split, remaining in (("first", [1, 5, 17, 27]), ("second", [3, 18, 29, 3])):
            for gene, full, retained in zip("ABCD", (10, 20, 30, 30), remaining, strict=True):
                rows.append({
                    "heldout_group": group, "heldout_mouse": "1", "evidence_split": split,
                    "v_gene": gene, "n_full": full, "n_g1_exclusive_aav": retained,
                    "only_g1_share": 0.9, "response": 0.2 if group == "g1" else -0.2,
                })
    selected = select_calibration_rows(pd.DataFrame(rows), None)
    assert selected.groupby("evidence_split")["v_gene"].unique().map(list).to_dict() == {
        "first": ["D"], "second": ["C"],
    }
    assert selected.groupby(["heldout_group", "evidence_split"])["regression_weight"].sum().to_list() == (
        pytest.approx([1 / 6] * 6)
    )


def test_one_se_ridge_choice_excludes_target_losses():
    from scripts.run_v_updates import select_ridge
    from src.eq7_calibration import RIDGE_GRID

    rows = []
    for ridge, score in zip(RIDGE_GRID, (0.50, 0.30, 0.32, 0.34, 0.39), strict=True):
        for group in ("g5", "g6"):
            for err in (score - 0.06, score + 0.06):
                rows.append({"ridge": ridge, "outer_group": group, "mse": err, "baseline_mse": 0.55})
        rows.append({"ridge": ridge, "outer_group": "g1", "mse": 1.0 if ridge < 10 else 0, "baseline_mse": 0.4})
    curve, chosen, status = select_ridge(pd.DataFrame(rows), None)
    assert status == "one_se_rule"
    assert chosen == 1.0
    assert curve.loc[curve["selected"], "mse"].item() == pytest.approx(0.34)
