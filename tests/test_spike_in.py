import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.spike_in import (
    apply_spike_to_tables, combine_count_tables, load_snapshot_tables, measure_case,
    select_targets, sensitivity_summary, snapshot_tables, validate_counts,
    summarize_experiment,
)

ROOT = Path(__file__).resolve().parents[1]


def matrix(keys, columns, values):
    return pd.DataFrame(values, index=pd.Index(keys, tupleize_cols=False), columns=columns)


def example():
    tables = {
        "g1": matrix([("CAAA", "TRAV1")], ["a", "b"], [[10000, 20000]]),
        "g3": matrix([("CASS", "TRAV9"), ("CATS", "TRAV9")], ["donor"], [[10], [20]]),
    }
    for g in ("g2", "g4", "g5", "g6"):
        tables[g] = matrix([("CAVV", "TRAV2")], [g], [[100]])
    retention = pd.DataFrame({"v": ["TRAV1"], "frequency_in_full_g1": [1.0],
                              "retained_cdr3_share": [1.0]})
    ranking = pd.DataFrame({"v_gene": ["TRAV1"], "edgeR_g1_enriched": [True],
                            "fisher_g1_enriched": [False]})
    selection = select_targets(tables, retention, ranking, n_clones=2)
    return tables, retention, ranking, selection


def test_selection_is_reproducible_and_excludes_candidates_and_control_shared_clones():
    tables, retention, ranking, selection = example()
    assert selection == select_targets(tables, retention, ranking, n_clones=2)
    assert selection["v_gene"] == "TRAV9"
    assert selection["baseline_absent_in_g1"] == 2
    tables["g2"] = tables["g3"].rename(columns={"donor": "control"})
    with pytest.raises(ValueError, match="No noncandidate"):
        select_targets(tables, retention, ranking, n_clones=2)
    with pytest.raises(ValueError, match="TRAV1 is a baseline candidate"):
        select_targets(tables, retention, ranking, v_gene="TRAV1", n_clones=1)


def test_addition_is_integer_family_dose_per_sample_and_preserves_controls():
    tables, _, _, selection = example()
    originals = {g: t.copy() for g, t in tables.items()}
    changed, dose = apply_spike_to_tables(tables, selection, 0.001)
    assert changed["g1"].sum().tolist() == [10010, 20020]
    assert dose.groupby("sample_id")["added_umi"].sum().to_dict() == {"a": 10, "b": 20}
    assert np.issubdtype(changed["g1"].values.dtype, np.integer)
    for g in tables:
        pd.testing.assert_frame_equal(tables[g], originals[g])
        if g != "g1":
            pd.testing.assert_frame_equal(changed[g], originals[g])


def test_zero_and_sub_umi_doses_do_not_create_phantom_presence_or_compound_doses():
    tables, _, _, selection = example()
    zero, _ = apply_spike_to_tables(tables, selection, 0)
    pd.testing.assert_frame_equal(zero["g1"], tables["g1"])
    tiny, dose = apply_spike_to_tables(tables, selection, 0.00001)
    assert dose["added_umi"].sum() == 0
    assert len(tiny["g1"]) == 1
    first, _ = apply_spike_to_tables(tables, selection, 0.001)
    second, _ = apply_spike_to_tables(tables, selection, 0.01)
    assert first["g1"].sum().tolist() == [10010, 20020]
    assert second["g1"].sum().tolist() == [10100, 20200]


@pytest.mark.parametrize("fraction", [-0.1, float("nan"), float("inf"), 1])
def test_invalid_doses_are_rejected(fraction):
    tables, _, _, selection = example()
    with pytest.raises(ValueError, match="fraction"):
        apply_spike_to_tables(tables, selection, fraction)


def test_rounding_does_not_force_one_umi_per_clone():
    tables, _, _, selection = example()
    changed, dose = apply_spike_to_tables(tables, selection, 0.00005)
    assert dose.groupby("sample_id")["added_umi"].sum().to_dict() == {"a": 1, "b": 1}
    assert len(changed["g1"]) == 2  # one family member remains truly absent


def test_snapshot_roundtrip_detects_changed_counts_and_metadata(tmp_path):
    tables, _, _, _ = example()
    md = pd.DataFrame({"sample_id": ["a"], "group_no": ["g1"], "mouse_id": ["m1"],
                       "source": ["thymus"], "subtype": ["cd4"]})
    snapshot_tables(tables, tmp_path, "TRA", md)
    restored = load_snapshot_tables(tmp_path, "TRA", md)
    for g in tables:
        pd.testing.assert_frame_equal(restored[g], tables[g])
    with pytest.raises(ValueError, match="metadata changed"):
        load_snapshot_tables(tmp_path, "TRA", md.assign(mouse_id="m2"))
    path = tmp_path / "count_snapshots/TRA_g1.csv.gz"
    path.write_bytes(b"changed")
    with pytest.raises(ValueError, match="snapshot changed"):
        load_snapshot_tables(tmp_path, "TRA")


def test_count_union_retains_injected_rows_for_inference_and_rejects_duplicate_samples():
    tables, _, _, selection = example()
    changed, _ = apply_spike_to_tables(tables, selection, 0.001)
    combined = combine_count_tables({g: changed[g] for g in ("g1", "g2", "g5", "g6")})
    assert combined["a"].sum() == 10010
    assert ("CASS", "TRAV9") in combined.index
    assert combined.at[("CASS", "TRAV9"), "g5"] == 0
    with pytest.raises(ValueError, match="more than one"):
        combine_count_tables({"one": tables["g1"], "two": tables["g1"]})


def _cell(identifier):
    notebook = json.loads((ROOT / "venn_original.ipynb").read_text())
    return "".join(next(c for c in notebook["cells"] if c["id"] == identifier)["source"])


def test_injected_counts_reach_real_notebook_exact_sets_and_pooled_inference(tmp_path):
    tables, _, _, selection = example()
    changed, _ = apply_spike_to_tables(tables, selection, 0.001)
    scope = {"pd": pd, "np": np}
    tree = ast.parse(_cell("venn-historical-code-09"))
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    exec(compile(ast.Module(body=functions, type_ignores=[]), "notebook_functions", "exec"), scope)
    comparison, _, remaining = scope["compare_full_vs_remaining_after_subtractions"](
        changed["g1"], [changed[g] for g in ("g2", "g4", "g5", "g6")], "g1")
    assert set(map(tuple, selection["clonotypes"])) <= remaining
    assert comparison.loc[comparison["v"].eq("TRAV9"), "retained_cdr3_share"].item() == 1
    md = pd.DataFrame([
        {"sample_id": sample, "group_no": g, "mouse_id": f"{g}_m1", "mouse_no": 1}
        for g, table in changed.items() for sample in table.columns
    ])
    scope.update({f"ct_{g}_aaV_tra": table for g, table in changed.items()})
    scope.update(metadata=md, POOL_WITHIN_MOUSE=True, RESULT_DIR=tmp_path,
                 combine_count_tables=combine_count_tables, display=lambda *a: None)
    # the cell's scipy/statsmodels imports belong to the next Fisher step
    code = _cell("venn-differential-input")
    code = "\n".join(line for line in code.splitlines() if not line.startswith("from scipy")
                     and not line.startswith("from statsmodels"))
    exec(code, scope)
    pooled = scope["ct_tra_aaV"]
    assert pooled["g1|g1_m1"].sum() == 30030
    assert pooled.at[("CASS", "TRAV9"), "g1|g1_m1"] == 15
    assert scope["sample_metadata"].index.is_unique
    assert pooled.shape[1] == 4


def test_recovery_reports_filtering_separately_from_failed_fdr():
    tables, _, _, selection = example()
    changed, _ = apply_spike_to_tables(tables, selection, 0.001)
    retention = pd.DataFrame({"v": ["TRAV9"], "cdr3_in_full_g1": [2],
        "cdr3_remaining_in_g1": [2], "retained_cdr3_share": [1.0],
        "frequency_in_full_g1": [2 / 3]})
    edger = pd.DataFrame({"feature_id": [("CASS", "TRAV9")], "logFC": [2.0], "FDR": [0.01]})
    ranking = pd.DataFrame({"v_gene": ["TRAV9"], "stratum_rank": [1]})
    metrics, recovery = measure_case(selection, 0.001, changed["g1"],
        set(changed["g1"].index), retention, edger, ranking)
    assert metrics["n_recovered_remainder"] == 2
    assert metrics["n_edger_tested_selected"] == 1
    assert metrics["n_edger_positive_fdr_selected"] == 1
    assert recovery["edger_tested"].tolist().count(False) == 1
    assert metrics["v_umi_share_g1"] == pytest.approx(30 / 30030)


def test_sensitivity_limits_distinguish_pass_fail_transitions_from_maximum_tested_dose():
    frame = pd.DataFrame({"fraction": [0, .0001, .001, .01],
        "all_selected_in_remainder": [False, True, True, True],
        "any_edger_positive_fdr": [False, False, True, True],
        "all_edger_positive_fdr": [False, False, True, False]})
    summary = sensitivity_summary(frame)
    results = {r["criterion"]: r for r in summary["criteria"]}
    any_positive = results["any_edger_positive_fdr"]
    assert any_positive["lowest_tested_passing_fraction"] == .001
    assert any_positive["largest_tested_failure_below_first_pass"] == .0001
    assert not any_positive["upper_failure_observed"]
    assert results["all_edger_positive_fdr"]["upper_failure_observed"]


def test_later_notebook_stratum_cells_call_analysis_instead_of_returning_lambdas():
    for path in (ROOT / "notebooks").glob("*.ipynb"):
        calls = []
        scope = {"STRATA": ["all_combined", "cd4_thymus"], "context": object(),
                 "repertoire": object(), "display": lambda x: None,
                 "run_stratum": lambda *args: calls.append(args[-1])}
        for cell in json.loads(path.read_text())["cells"]:
            source = "".join(cell["source"])
            if cell["cell_type"] == "code" and source.startswith("if "):
                exec(source, scope)
        assert calls == ["all_combined", "cd4_thymus"]


def test_summary_records_entry_from_unranked_baseline_without_inventing_a_rank(tmp_path):
    baseline = {"fraction": 0, "v_gene": "TRAV9", "n_selected": 2,
        "v_remaining_g1": 0, "v_retention": 0, "v_unique_share_g1": 0,
        "v_umi_share_g1": 0, "v_article_score": 0, "v_final_rank": None,
        "n_edger_positive_fdr_selected": 0, "recovery_fraction": 0,
        "edger_positive_recovery_fraction": 0, "all_selected_in_remainder": False}
    dose = {**baseline, "fraction": .001, "v_final_rank": 3,
            "recovery_fraction": 1, "all_selected_in_remainder": True}
    for case, metrics in [("baseline", baseline), ("dose_01", dose)]:
        folder = tmp_path / case
        folder.mkdir()
        (folder / "spike_metrics.json").write_text(json.dumps(metrics))
    summary = summarize_experiment(tmp_path, tmp_path / "figures", ["baseline", "dose_01"])
    assert pd.isna(summary.loc[0, "v_final_rank"])
    assert pd.isna(summary.loc[1, "rank_improvement"])
    assert summary.loc[1, "entered_ranking"]
    assert summary.loc[1, "rank_improved"]
    assert (tmp_path / "figures/spike_in_sensitivity.pdf").is_file()
