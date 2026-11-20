import pandas as pd
import pytest

from scripts.run_eq7_external import read_trust4_report
from src.eq7_calibration import FEATURES, score_candidates
from src.eq7_external import compare_shortlists, panel_cutoffs, summarize_transfer


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
        "#count\tCDR3aa\tV\n"
        "40\tCASSF\tTRBV2*01|TRBV2*02\n"
        "10\tCASSF\tTRBV2\n"
        "30\tCASSF\tTRBV2*01|TRBV4*02\n"
        "70\tCASS*\tTRBV2*01\n"
        "11\tCASSF\tTRAV2*01\n"
    )
    counts, qc = read_trust4_report(report, "test")
    assert counts.to_dict() == {("CASSF", "TRBV2"): 50}
    assert qc["n_productive_trb_read_fragments"] == 50
