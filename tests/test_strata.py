import pandas as pd
import pytest

from src.strata import (
    _infer_cell_subset,
    analysis_metadata,
    select_stratum,
    validate_strata,
)


def _annotated_rows():
    return pd.DataFrame(
        [
            {
                "sample_id": "g1_m1_thymus_cd4",
                "mouse_id": "m1",
                "group": "g1",
                "_cell_subset": "cd4",
                "_tissue": "thymus",
                "cdr3": "CASSL",
                "v_gene": "TRAV1",
                "umi": 10,
                "ckey": "CASSL|TRAV1",
            },
            {
                "sample_id": "g1_m1_spleen_cd4",
                "mouse_id": "m1",
                "group": "g1",
                "_cell_subset": "cd4",
                "_tissue": "spleen",
                "cdr3": "CASSP",
                "v_gene": "TRAV2",
                "umi": 8,
                "ckey": "CASSP|TRAV2",
            },
            {
                "sample_id": "g5_m1_thymus_cd4",
                "mouse_id": "m1",
                "group": "g5",
                "_cell_subset": "cd4",
                "_tissue": "thymus",
                "cdr3": "CASSQ",
                "v_gene": "TRAV3",
                "umi": 7,
                "ckey": "CASSQ|TRAV3",
            },
            {
                "sample_id": "g5_m1_spleen_cd4",
                "mouse_id": "m1",
                "group": "g5",
                "_cell_subset": "cd4",
                "_tissue": "spleen",
                "cdr3": "CASSR",
                "v_gene": "TRAV4",
                "umi": 6,
                "ckey": "CASSR|TRAV4",
            },
        ]
    )


def test_combined_stratum_pools_tissues_within_group_and_mouse():
    combined = select_stratum(_annotated_rows(), "cd4_combined")

    assert combined["sample_id"].nunique() == 4
    assert combined["analysis_unit"].nunique() == 2
    assert set(combined["analysis_unit"]) == {
        "g1|m1|cd4|combined",
        "g5|m1|cd4|combined",
    }


def test_tissue_stratum_keeps_sample_level_units():
    thymus = select_stratum(_annotated_rows(), "cd4_thymus")

    assert thymus["analysis_unit"].tolist() == thymus["sample_id"].tolist()


def test_analysis_metadata_rejects_conflicting_units():
    frame = select_stratum(_annotated_rows(), "cd4_thymus")
    frame.loc[frame.index[-1], "analysis_unit"] = frame.iloc[0]["analysis_unit"]

    with pytest.raises(ValueError, match="inconsistent"):
        analysis_metadata(frame)


def test_validate_strata_accepts_requested_subset():
    report = validate_strata(
        _annotated_rows(),
        ("cd4_thymus", "cd4_spleen", "cd4_combined"),
    )

    assert report["stratum"].tolist() == [
        "cd4_thymus",
        "cd4_spleen",
        "cd4_combined",
    ]
    assert (
        report.loc[
            report["stratum"].eq("cd4_combined"),
            "n_analysis_units",
        ].item()
        == 2
    )


def test_ambiguous_explicit_subset_is_not_overridden_by_sample_name():
    frame = pd.DataFrame(
        {
            "subtype": ["CD4 CD8"],
            "sample_id": ["g1_m1_thymus_cd4"],
        }
    )

    assert _infer_cell_subset(frame).isna().all()
