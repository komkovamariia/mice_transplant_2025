import numpy as np
import pandas as pd

from src.approach4_cross_approach import (
    _consensus,
    _direction_agreement,
)


def _ranking(name, effects):
    return pd.DataFrame(
        {
            "v_gene": ["TRAV1", "TRAV2"],
            "rank": [1, 2],
            f"{name}_effect_g1_vs_allogeneic": effects,
        }
    )


def test_direction_agreement_ignores_zero_and_missing_effects():
    n, agreement = _direction_agreement(
        pd.Series([1.0, -2.0, 0.0, np.nan]),
        pd.Series([3.0, 4.0, 2.0, -1.0]),
    )

    assert n == 2
    assert agreement == 0.5


def test_consensus_preserves_directional_disagreement():
    rankings = {
        "set_count": _ranking("set_count", [1.0, -1.0]),
        "sequence_embedding": _ranking(
            "sequence_embedding",
            [2.0, -2.0],
        ),
        "clone_alloreactivity": _ranking(
            "clone_alloreactivity",
            [-0.5, -3.0],
        ),
    }

    consensus = _consensus(rankings, "cd4_thymus").set_index("v_gene")

    assert consensus.loc["TRAV1", "consensus_direction"] == "g1_enriched"
    assert consensus.loc["TRAV1", "direction_agreement_fraction"] == 2 / 3
    assert consensus.loc["TRAV2", "consensus_direction"] == "allogeneic_enriched"
    assert consensus.loc["TRAV2", "direction_agreement_fraction"] == 1.0
