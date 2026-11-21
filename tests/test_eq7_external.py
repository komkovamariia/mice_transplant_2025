import pandas as pd
import pytest

from scripts.run_eq7_external import read_trust4_report
from src.eq7_calibration import FEATURES, score_candidates
from src.eq7_external import compare_shortlists, mouse_prediction_losses, mouse_predictions, panel_cutoffs, summarize_transfer


def test_external_shortlist_uses_only_the_frozen_internal_gate():
    internal = pd.DataFrame({
        "v_gene": list("ABCD"), "n_full": [10] * 4, "n_g1_exclusive_aav": [1, 5, 7, 9],
        "only_g1_share": [0.1, 0.5, 0.7, 0.9], FEATURES[0]: [0] * 4,
        FEATURES[1]: [0] * 4, FEATURES[2]: [1.0] * 4,
    })
    external = internal.assign(n_g1_exclusive_aav=[1, 1, 3, 5])
    weights = pd.DataFrame({"term": FEATURES, "weight": [0, 0, 1], "intercept": 0})
    assert panel_cutoffs(external) != panel_cutoffs(internal)
    assert list(score_candidates(internal, weights)["v_gene"]) == ["D"]
    assert list(score_candidates(external, weights)["v_gene"]) == ["D"]
    assert score_candidates(external, weights, frozen_cutoffs=panel_cutoffs(internal)).empty


def test_union_loss_charges_membership_once_and_handles_zero_scores():
    internal = pd.DataFrame({"v_gene": list("ABCD"), "E": [2.0, 0.0, 4.0, 1.0]})
    external = pd.DataFrame({"v_gene": list("ABCE"), "E": [1.0, 0.0, 4.0, 0.2]})
    profile = pd.DataFrame({
        "v_gene": list("ABCDE"), "n_full": [10] * 5, "n_g1_exclusive_aav": [8] * 5,
    })
    compared = compare_shortlists(internal, external, profile, profile).set_index("v_gene")
    assert compared["union_loss"].sort_index().to_list() == pytest.approx([0.5, 0, 0, 1, 1])
    assert pd.isna(compared.loc["D", "absolute_delta_E"])
    losses = pd.DataFrame({
        "heldout_group": ["g1", "g1", "g5"], "mae": [0.2, 0.3, 0.1], "baseline_mae": [0.4, 0.3, 0.2],
    })
    summary = summarize_transfer(losses, compared)
    assert summary["jaccard"] == 3 / 5
    assert summary["mean_union_loss"] == 0.5
    assert summary["balanced_mouse_mae"] == pytest.approx(0.175)


def test_trust4_counts_keep_productive_unambiguous_trb(tmp_path):
    report = tmp_path / "report.tsv"
    report.write_text(
        "#count\tCDR3aa\tV\tJ\n"
        "40\tCASSF\tTRBV2*01|TRBV2*02\tTRBJ1-1*01\n"
        "10\tCASSF\tTRBV2\tTRBJ1-1*01\n"
        "30\tCASSF\tTRBV2*01|TRBV4*02\tTRBJ1-1*01\n"
        "70\tCASS*\tTRBV2*01\tTRBJ1-1*01\n"
        "11\tCASSF\tTRAV2*01\tTRAJ1*01\n"
        "20\tCAVF\tTRAV14D-3/DV8*01\tTRAJ1*01\n"
        "50\tCAVF\tTRAV14D-3/DV8*01\tTRDJ1*01\n"
    )
    counts, qc = read_trust4_report(report, "test")
    assert counts.to_dict() == {("CASSF", "TRBV2"): 50}
    assert qc["n_productive_trb_read_fragments"] == 50
    counts, qc = read_trust4_report(report, "test", "TRA")
    assert counts.to_dict() == {("CASSF", "TRAV2"): 11, ("CAVF", "TRAV14D-3-DV8"): 20}
    assert qc["n_productive_tra_read_fragments"] == 31


def test_single_arm_leaves_unobserved_segments_and_low_depth_mice_unscored():
    index = pd.Index([("CAVF", "TRAV2"), ("CAWF", "TRAV2")], tupleize_cols=False)
    g1 = pd.DataFrame({"internal": [1]}, index=index[:1])
    empty = g1.iloc[:0, :0]
    reference = {"g1": g1, "g5": empty, "g6": empty}
    external = {"g5": pd.DataFrame({"passed": [1, 1]}, index=index)}
    metadata = pd.DataFrame({
        "group_no": ["g5", "g5"], "mouse_id": ["one", "two"], "sample_id": ["passed", "low_depth"],
    })
    candidates = pd.DataFrame({"v_gene": ["TRAV2", "TRAV3"], "predicted_response": [-0.4, 0.2]})
    predictions = mouse_predictions(reference, candidates, external, metadata, baseline_prediction=0)
    assert predictions["response"].iat[0] == -0.5
    assert pd.isna(predictions["response"].iat[1])
    losses = mouse_prediction_losses(predictions, metadata)
    assert losses["n_observed"].tolist() == [1, 0]
    assert losses["mae"].iat[0] == pytest.approx(0.1)
    assert pd.isna(losses["mae"].iat[1])
