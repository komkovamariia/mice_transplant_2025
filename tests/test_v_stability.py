import pandas as pd
import pytest

from src.v_stability import (
    adaptive_profile, compare_shortlists, rarefaction_profiles, representative_seed,
)
from src.v_thresholds import PlotThresholds


@pytest.mark.parametrize(("chain", "size"), [("TRA", 14), ("TRB", 4)])
def test_rarefaction_retains_the_complete_adaptive_shortlist(chain, size):
    genes = [f"{chain}V{i + 1}" for i in range(size + 1)]
    features = pd.Index([(f"seq{i}", gene) for i, gene in enumerate(genes)], tupleize_cols=False)
    tables = {
        "g1": pd.DataFrame({"target": 1}, index=features),
        **{f"g{i}": pd.DataFrame({f"background{i}": len(features)}, index=features[-1:]) for i in (2, 4, 5, 6)},
    }
    profile, cutoffs = adaptive_profile(tables)
    reference = profile.loc[profile["eligible"], "v"].tolist(), cutoffs
    _, shortlists = rarefaction_profiles(tables, depth=len(features), seeds=(1201, 1202))
    comparisons = compare_shortlists(reference, shortlists)
    assert len(reference[0]) == size
    assert comparisons["n_recovered"].tolist() == [size, size]
    assert comparisons["jaccard"].tolist() == [1.0, 1.0]
    assert representative_seed(shortlists)[0] == 1201


def test_recovery_separates_unavailable_cutoffs_from_an_empty_fit():
    cutoffs = PlotThresholds(1.0, 70.0)
    result = compare_shortlists((["V1"], cutoffs), {1: ([], None), 2: ([], cutoffs)})
    assert result["status"].tolist() == ["no_valid_cutoffs", "compared"]
    assert pd.isna(result["n_recovered"].iat[0])
    assert result["n_recovered"].iat[1] == 0
