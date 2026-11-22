"""Mouse-held-out V-segment recovery and its abundance-weighted contrast."""

from dataclasses import dataclass
from itertools import product

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from .eq7_calibration import FEATURES, shortlist
from .spike_in import aa_v_key


GROUPS = ("g1", "g5", "g6")
COUNT_CONTRAST = "count_contrast_pp"
MOLECULE_CONTRAST = "weighted_molecules_per_million_contrast"


def pool_mice(tables, metadata):
    samples = metadata.set_index("sample_id")["mouse_id"]
    pooled = {}
    for group in GROUPS:
        table = tables.get(group)
        if table is None:
            continue
        sample_mice = samples.loc[table.columns]
        for mouse, libraries in sample_mice.groupby(sample_mice):
            counts = table.loc[:, libraries.index].sum(axis=1)
            pooled[(group, mouse)] = counts.loc[counts.gt(0)]
    return pooled


class RecoveryPanel:
    """A fixed cohort of depth-standardized libraries pooled within mouse."""

    def __init__(self, counts, *, align_to=None):
        self.mice = tuple(counts)
        if align_to is None:
            focal = [counts[key].index for key in self.mice if key[0] == "g1"]
            if not focal:
                raise ValueError("The reference panel needs at least one g1 mouse.")
            self.keys = focal[0]
            for index in focal[1:]:
                self.keys = self.keys.union(index, sort=False)
            self.codes, self.genes = pd.factorize(pd.Index([aa_v_key(key)[1] for key in self.keys]))
        else:
            self.keys, self.codes, self.genes = align_to.keys, align_to.codes, align_to.genes
        self.counts = np.stack([counts[key].reindex(self.keys, fill_value=0).to_numpy() for key in self.mice])
        self.depth = np.array([counts[key].sum() for key in self.mice])
        self.unique_v = np.stack([
            pd.Index([aa_v_key(k)[1] for k in counts[key].index]).value_counts()
            .reindex(self.genes, fill_value=0).to_numpy()
            for key in self.mice
        ])
        self.present = self.counts > 0

    def by_group(self, omit=()):
        return {
            group: tuple(i for i, key in enumerate(self.mice) if key[0] == group and i not in omit)
            for group in GROUPS
        }

    def reference(self, omit=()):
        mice = self.by_group(omit)
        focal = mice["g1"]
        background = (*mice["g5"], *mice["g6"])
        if not focal or not background:
            raise ValueError("The recovery reference requires retained g1 and comparison mice.")
        prevalence = self.present[list(focal)].mean(axis=0)
        prevalence[self.present[list(background)].any(axis=0)] = 0
        return prevalence

    def measure(self, reference, query):
        """Return unique-aaV recovery in percent and weighted counts per million."""
        present = reference > 0
        codes = self.codes[present]
        n_genes = len(self.genes)
        rows = []
        for i in query:
            unique = np.bincount(codes, weights=self.present[i, present], minlength=n_genes)
            denominator = self.unique_v[i]
            fraction = np.divide(100 * unique, denominator, out=np.full(n_genes, np.nan), where=denominator > 0)
            molecules = np.bincount(
                codes, weights=self.counts[i, present] * reference[present], minlength=n_genes,
            )
            rows.append(np.stack([fraction, 1e6 * molecules / self.depth[i]]))
        return np.stack(rows)


def contrasts(measurements):
    """Use one row per group in g1, g5, g6 order."""
    return measurements[0] - measurements[1:].mean(axis=0)


def training_contrasts(panel, omitted, measurements, *, recovery_column=MOLECULE_CONTRAST):
    """Build outcomes from a symmetric triplet of retained animals at a time."""
    by_group = panel.by_group(omitted)
    triplets = list(product(*(by_group[group] for group in GROUPS)))
    if not triplets:
        raise ValueError("Every group must have an evaluable training mouse.")
    evaluations = []
    for triplet in triplets:
        exclusions = tuple(sorted((*omitted, *triplet)))
        if exclusions not in measurements:
            measurements[exclusions] = dict(zip(
                exclusions, panel.measure(panel.reference(exclusions), exclusions), strict=True,
            ))
        evaluations.append([measurements[exclusions][mouse] for mouse in triplet])
    measures = np.array(evaluations)
    valid = np.isfinite(measures)
    group_means = np.divide(
        np.nansum(measures, axis=0), valid.sum(axis=0),
        out=np.full(measures.shape[1:], np.nan), where=valid.any(axis=0),
    )
    signed = (group_means[0, 0] - group_means[1, 0] - group_means[2, 0]) / 3
    count, molecule = contrasts(group_means)
    return pd.DataFrame({
        "v_gene": panel.genes, "signed_recovery_pct": signed,
        COUNT_CONTRAST: count, recovery_column: molecule,
    })


@dataclass(frozen=True)
class RidgeFit:
    coefficients: pd.Series
    intercept: float
    baseline: float

    def predict(self, rows):
        return self.intercept + rows[self.coefficients.index].to_numpy() @ self.coefficients.to_numpy()


def fit_ridge(rows, terms, target, penalty, *, nonnegative=False, device="cpu", weight=None):
    """Weighted ridge on training-standardized predictors."""
    train = rows.loc[rows[target].notna()]
    if train.empty:
        raise ValueError(f"No training outcomes for {target}.")
    x = train[list(terms)].to_numpy(dtype=float)
    y = train[target].to_numpy(dtype=float)
    w = train[weight].to_numpy(dtype=float) if weight else np.ones(len(train))
    w = w / w.sum()
    means = w @ x
    baseline = float(w @ y)
    scales = np.sqrt(w @ (x - means) ** 2)
    active = scales > 1e-12
    beta = np.zeros(len(terms))
    if active.any():
        z = (x[:, active] - means[active]) / scales[active]
        design = z * np.sqrt(w[:, None])
        response = (y - baseline) * np.sqrt(w)
        regularizer = np.eye(active.sum())
        if nonnegative:
            slopes = nnls(
                np.vstack([design, np.sqrt(penalty) * regularizer]),
                np.r_[response, np.zeros(active.sum())],
            )[0]
        else:
            gram = design.T @ design + penalty * regularizer
            rhs = design.T @ response
            if device == "cuda":
                import torch

                slopes = torch.linalg.solve(
                    torch.as_tensor(gram, device="cuda"), torch.as_tensor(rhs, device="cuda"),
                ).cpu().numpy()
            else:
                slopes = np.linalg.solve(gram, rhs)
        beta[active] = slopes / scales[active]
    return RidgeFit(pd.Series(beta, index=list(terms)), float(baseline - means @ beta), baseline)


def top_weights(ranking, k):
    """Expected inclusion with uniform tie-breaking at the last position."""
    ordered = np.sort(ranking)[::-1]
    cutoff = ordered[min(k, len(ordered)) - 1]
    better = ranking > cutoff
    tied = ranking == cutoff
    return better.astype(float) + tied * ((min(k, len(ordered)) - better.sum()) / tied.sum())


def split_shortlist(rows):
    fixed = rows.drop(columns=["response", "heldout_n_v", "heldout_n_recovered"])
    return shortlist(fixed.drop_duplicates("v_gene")).reset_index(drop=True)


def full_training_profiles(panel, folds, penalties, device):
    """Fit the external transport model from the 75 training-only triplets."""
    model_inputs = []
    lookup = {key: idx for idx, key in enumerate(panel.mice)}
    for rows in folds.values():
        mice = rows[["heldout_group", "heldout_mouse"]].drop_duplicates().set_index("heldout_group")
        query = tuple(lookup[(g, mice.at[g, "heldout_mouse"])] for g in GROUPS)
        values = panel.measure(panel.reference(query), query)
        count, abundance = contrasts(values)
        eligible = split_shortlist(rows)
        eligible[COUNT_CONTRAST] = eligible["v_gene"].map(pd.Series(count, index=panel.genes))
        eligible[MOLECULE_CONTRAST] = eligible["v_gene"].map(pd.Series(abundance, index=panel.genes))
        model_inputs.append(eligible.assign(fit_weight=1 / len(eligible)))
    train = pd.concat(model_inputs, ignore_index=True)
    return {
        target: {penalty: fit_ridge(train, FEATURES, target, penalty, device=device, weight="fit_weight")
                 for penalty in penalties}
        for target in (COUNT_CONTRAST, MOLECULE_CONTRAST)
    }
