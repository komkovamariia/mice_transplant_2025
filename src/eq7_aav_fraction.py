"""Fraction of statistically supported g1-enriched aaV per V segment.

This uses existing raw UMI snapshots and held-out mouse identifiers.
It never creates rarefactions, splits or structural shortlists.
"""

import numpy as np
import pandas as pd

from .differential_abundance import fit_edger, fit_fisher, prepare_counts
from .eq7_calibration import FEATURES
from .spike_in import EDGER_SUCCESS, aa_v_key


FRACTION = ("edgeR_positive_fraction", "Fisher_positive_fraction", FEATURES[2])
FRACTION_MODEL = "aav_positive_fraction_ridge"


def prepare_aav_counts(raw_tables, metadata):
    """Pool *raw* functional UMI counts within each mouse, keeping aaV keys."""
    empty = raw_tables["g1"].iloc[:, :0]
    groups = {group: raw_tables.get(group, empty) for group in ("g1", "g2", "g5", "g6")}
    return prepare_counts(groups, metadata, pool_within_mouse=True)


def summarize_positive_fraction(results, fdr_col, effect_col, all_genes):
    """Positive FDR-significant aaV / all aaV tested by the same method.

    EdgeR's denominator excludes aaV removed by its expression filter.
    Fisher's denominator includes all aaV with positive counts in the two
    contrasted groups. Zero tested aaV gives a zero feature.
    """
    frame = results[["feature_id", fdr_col, effect_col]].copy()
    frame["v_gene"] = frame["feature_id"].map(lambda x: aa_v_key(x)[1])
    frame["positive"] = frame[fdr_col].lt(0.05) & frame[effect_col].gt(0)
    tested = frame.groupby("v_gene", observed=True).size()
    positive = frame.groupby("v_gene", observed=True)["positive"].sum()
    tested = tested.reindex(all_genes, fill_value=0).astype(int)
    positive = positive.reindex(all_genes, fill_value=0).astype(int)
    fraction = np.divide(
        positive.to_numpy(dtype=float), tested.to_numpy(dtype=float),
        out=np.zeros(len(all_genes), dtype=float), where=tested.to_numpy() > 0,
    )
    return fraction, positive.to_numpy(), tested.to_numpy()


def fraction_evidence(counts, units, excluded=(), *, workers=4):
    """Calculate fold-pure aaV-level significance fractions on retained mice."""
    keys = pd.MultiIndex.from_frame(units[["group_no", "mouse_id"]])
    remaining = units.loc[~keys.isin(excluded)]
    contrast_cols = remaining.index[remaining["contrast_group"].isin(("g1", "allogeneic"))]
    train = counts.loc[:, remaining.index]
    # Exclude aaV absent from both contrasted arms (including g2-only aaV).
    train = train.loc[train[contrast_cols].sum(axis=1).gt(0)]
    genes = pd.Index(
        sorted({aa_v_key(key)[1] for key in train.index}), name="v_gene",
    )
    edge, status = fit_edger(train, remaining)
    if status != EDGER_SUCCESS:
        raise ValueError(f"aaV-level edgeR unavailable: {status}")
    fisher = fit_fisher(train, remaining, n_workers=workers)
    edge_frac, edge_pos, edge_tested = summarize_positive_fraction(
        edge, "FDR", "logFC", genes,
    )
    fisher_frac, fisher_pos, fisher_tested = summarize_positive_fraction(
        fisher, "Fisher_FDR_g1", "Fisher_log2OR_g1", genes,
    )
    return pd.DataFrame({
        "v_gene": genes,
        FRACTION[0]: edge_frac,
        FRACTION[1]: fisher_frac,
        "n_edgeR_positive": edge_pos,
        "n_edgeR_tested": edge_tested,
        "n_Fisher_positive": fisher_pos,
        "n_Fisher_tested": fisher_tested,
    })
