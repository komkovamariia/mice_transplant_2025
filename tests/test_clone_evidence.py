import numpy as np
import pandas as pd

import src.approach3_clone_alloreactivity as approach3


def _row(group, mouse, sample, cdr3, v_gene, umi):
    return {
        "group": group,
        "mouse_id": mouse,
        "sample_id": sample,
        "_cell_subset": "cd4",
        "_tissue": "thymus",
        "cdr3": cdr3,
        "v_gene": v_gene,
        "umi": umi,
        "ckey": f"{cdr3}|{v_gene}",
    }


def test_abundance_score_weights_mice_equally(tmp_path, monkeypatch):
    repertoire = pd.DataFrame(
        [
            _row("g1", "m1", "g1_m1", "CASSL", "TRAV1", 10),
            _row("g1", "m1", "g1_m1", "CASSP", "TRAV2", 90),
            _row("g5", "m2", "g5_m2", "CASSL", "TRAV1", 10),
            _row("g5", "m2", "g5_m2", "CASSP", "TRAV2", 90),
            _row("g6", "m3", "g6_m3", "CASSL", "TRAV1", 100),
            _row("g6", "m3", "g6_m3", "CASSP", "TRAV2", 900),
        ]
    )
    monkeypatch.setattr(approach3, "_output_dir", lambda: tmp_path)
    monkeypatch.setattr(
        approach3,
        "_plot_candidate_scores",
        lambda candidates, stratum: None,
    )

    approach3.run_stratum(repertoire, "cd4_thymus")
    evidence = pd.read_csv(
        tmp_path / "cd4_thymus_all_clonotype_evidence.csv"
    ).set_index("ckey")

    effect = evidence.loc[
        "CASSL|TRAV1",
        "log2_relative_abundance_fc_allogeneic_vs_g1",
    ]
    assert np.isclose(effect, 0.0)
