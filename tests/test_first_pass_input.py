from pathlib import Path

import pandas as pd
import pytest

from src.count_statistics import fisher_log2_odds_ratio
from src.first_pass_input import count_table_to_repertoire, prepare_sample_index


def _write_inputs(tmp_path: Path):
    metadata = pd.DataFrame(
        [
            {
                "chain": "alpha",
                "sample_no": 100,
                "sample_id": "g1_m1_thymus_cd4_alpha",
                "group_no": "g1",
                "mouse_no": 1,
                "source": "thymus",
                "subtype": "cd4",
            },
            {
                "chain": "alpha",
                "sample_no": 101,
                "sample_id": "g5_m2_thymus_cd4_alpha",
                "group_no": "g5",
                "mouse_no": 2,
                "source": "thymus",
                "subtype": "cd4",
            },
            {
                "chain": "beta",
                "sample_no": 100,
                "sample_id": "g1_m1_thymus_cd4_beta",
                "group_no": "g1",
                "mouse_no": 1,
                "source": "thymus",
                "subtype": "cd4",
            },
        ]
    )
    metadata_path = tmp_path / "metadata.csv"
    metadata.to_csv(metadata_path, index=False)

    export_1 = tmp_path / "sample1.tsv"
    export_2 = tmp_path / "sample2.tsv"
    export_3 = tmp_path / "sample3.tsv"
    export_1.touch()
    export_2.touch()
    export_3.touch()
    index = pd.DataFrame(
        [
            {"sample_id": "alpha-100", "filename": export_1, "chain": "TRA"},
            {"sample_id": "alpha-101", "filename": export_2, "chain": "TRA"},
            {"sample_id": "beta-100", "filename": export_3, "chain": "TRB"},
        ]
    )
    index_path = tmp_path / "clonosets.csv"
    index.to_csv(index_path, index=False)
    return metadata_path, index_path


def test_prepare_sample_index_reproduces_historical_id_merge(tmp_path):
    metadata_path, index_path = _write_inputs(tmp_path)

    sample_index = prepare_sample_index(metadata_path, index_path)

    assert sample_index["sample_id"].tolist() == [
        "g1_m1_thymus_cd4_alpha",
        "g5_m2_thymus_cd4_alpha",
    ]
    assert sample_index["group_no"].tolist() == ["g1", "g5"]
    assert sample_index["chain"].tolist() == ["TRA", "TRA"]


def test_count_table_conversion_keeps_positive_functional_trav_rows(tmp_path):
    metadata_path, index_path = _write_inputs(tmp_path)
    sample_index = prepare_sample_index(metadata_path, index_path)
    counts = pd.DataFrame(
        {
            "g1_m1_thymus_cd4_alpha": [4, 0, 7],
            "g5_m2_thymus_cd4_alpha": [1, 3, 2],
        },
        index=[
            ("CAVRDSNYQLIW", "TRAV4-2"),
            ("CAVRGSALGRLHF", "TRAV7D-3"),
            ("CAV*", "TRAV1"),
        ],
    )

    repertoire = count_table_to_repertoire(counts, sample_index)

    assert set(repertoire["v_gene"]) == {"TRAV4-2", "TRAV7D-3"}
    assert repertoire["umi"].sum() == 8
    assert set(repertoire["group"]) == {"g1", "g5"}
    assert repertoire["ckey"].str.startswith("TRA|").all()


def test_fisher_effect_uses_relative_abundance():
    assert fisher_log2_odds_ratio(100, 200, 10_000, 40_000) > 0
    assert fisher_log2_odds_ratio(1, 9, 100, 100) < 0
    assert fisher_log2_odds_ratio(0, 0, 100, 200) == 0
    assert fisher_log2_odds_ratio(7, 3, 10, 10) == pytest.approx(
        2 * __import__("math").log2(7.5 / 3.5)
    )



def test_fisher_p_value_uses_uncorrected_counts():
    from scipy.stats import fisher_exact
    from src.count_statistics import pooled_fisher_rows

    g1 = pd.Series([100, 9900])
    graft = pd.Series([200, 39800])
    rows = pooled_fisher_rows(["a", "b"], g1, graft, n_jobs=1)
    assert rows[0]["Fisher_log2OR_g1"] > 0
    assert rows[0]["Fisher_p_g1"] == pytest.approx(
        fisher_exact([[100, 200], [9900, 39800]])[1]
    )
