"""Continuous aaV evidence and a joint abundance and diversity endpoint."""

import numpy as np
import pandas as pd

from .differential_abundance import fit_edger, fit_fisher, prepare_counts
from .eq7_calibration import FEATURES
from .eq7_weighted import RecoveryPanel
from .spike_in import EDGER_SUCCESS, aa_v_key


JOINT_CONTRAST = "joint_recovery_contrast"
CONTINUOUS = ("edgeR_v_effect", "fisher_v_effect", FEATURES[2])
BINARY_MODEL = "binary_evidence_ridge"
CONTINUOUS_MODEL = "continuous_evidence_ridge"


class JointRecoveryPanel(RecoveryPanel):
    def measure(self, reference, query):
        """Measure the breadth-weighted geometric mean at mouse level."""
        measured = super().measure(reference, query)
        # An absent segment has no recovered molecules or aaV.
        breadth = np.nan_to_num(measured[:, 0] / 100, nan=0)
        measured[:, 1] = np.sqrt(measured[:, 1] * breadth)
        return measured


def prepare_v_evidence(tables, metadata):
    """Sum raw aaV counts by V segment before pooling libraries within mouse."""
    grouped = {}
    for group in ("g1", "g5", "g6"):
        if group not in tables:
            continue
        table = tables[group]
        genes = table.index.map(lambda key: ("V", aa_v_key(key)[1]))
        grouped[group] = table.groupby(genes, sort=False).sum()
    empty = grouped["g1"].iloc[:, :0]
    analysis = {group: grouped.get(group, empty) for group in ("g1", "g2", "g5", "g6")}
    return prepare_counts(analysis, metadata, pool_within_mouse=True)


def continuous_evidence(counts, units, excluded=()):
    """Compute confidence-weighted effects of V usage on retained mice."""
    keys = pd.MultiIndex.from_frame(units[["group_no", "mouse_id"]])
    train_units = units.loc[~keys.isin(excluded)]
    train = counts.loc[:, train_units.index]
    train = train.loc[train.sum(axis=1).gt(0)]
    edge, status = fit_edger(train, train_units)
    if status != EDGER_SUCCESS:
        raise ValueError(status)
    fisher = fit_fisher(train, train_units, n_workers=1)
    a, b = fisher["fisher_g1_count"], fisher["fisher_allogeneic_count"]
    odds = np.log2(((a + 0.5) / (a.sum() - a + 0.5)) / ((b + 0.5) / (b.sum() - b + 0.5)))
    weights = pd.DataFrame({
        CONTINUOUS[0]: (edge["logFC"] * (1 - edge["FDR"])).set_axis(edge["feature_id"]),
        CONTINUOUS[1]: (odds * (1 - fisher["Fisher_FDR_g1"])).set_axis(fisher["feature_id"]),
    }).reindex(train.index).fillna(0)
    weights.index = weights.index.map(lambda key: aa_v_key(key)[1])
    return weights.rename_axis("v_gene").reset_index()
