"""Recompute the structural, abundance and selection results for the splenic manuscript."""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_eq7_splenic import DEPTH, SEED, matched_internal, prepare_mouse_evidence, training_predictors
from src.article_figures import plot_group_overlap
from src.differential_abundance import fit_edger, fit_fisher, rank_v_segments
from src.eq7_calibration import FEATURES, shortlist
from src.eq7_continuous import JointRecoveryPanel
from src.eq7_weighted import GROUPS, pool_mice, top_weights
from src.repertoire_depth import depth_check, mds_points
from src.spike_in import EDGER_SUCCESS, SUBTRACTION_GROUPS, load_snapshot_tables
from src.structural_sampling import _stable_sample_seed, downsample_group_tables
from src.structural_v import structural_remainder, subset_v_tables, v_plot_profile
from src.v_stability import RAREFACTION_SEEDS, compare_shortlists, representative_seed, regional_tables
from src.v_thresholds import fit_mouse_dropout_thresholds, fit_plot_thresholds

DISPLAY_SEED = 1249
PANELS = ("cd4", "cd8", "pooled")
METHODS = {
    "structural_screen": "Adaptive two-cutoff shortlist",
    "S_ungated": r"Structural rank $S(v)$",
    "initial_share_ungated": "Initial V-segment share",
    "edgeR_binary_ungated": "edgeR aaV indicator",
    "Fisher_binary_ungated": "Fisher aaV indicator",
    "old_E_heuristic_ungated": "Original fixed-weight score",
    "V_usage_edgeR_logFC": "edgeR V-segment effect",
    "V_usage_edgeR_significant": "edgeR V-segment effect with FDR",
}


def prepared_profiles(raw, metadata, chain, output):
    profiles, tables_by_subset, audit_rows = {}, {}, []
    for subset in PANELS:
        selected, audit = matched_internal(raw, metadata, subset)
        tables, _ = downsample_group_tables(selected, depth=DEPTH, seed=DISPLAY_SEED)
        profile = v_plot_profile(tables)
        cutoff, rows = fit_mouse_dropout_thresholds(tables, metadata)
        audit_rows.extend({"chain": chain, "subset": subset, **row} for row in rows)
        profile["selected"] = cutoff.includes(profile["frequency_in_full_g1_pct"],
                                               profile["retained_cdr3_share_pct"])
        profile.to_csv(output / f"profile_{subset}.csv", index=False)
        profiles[subset] = profile
        tables_by_subset[subset] = tables
        if subset == "pooled":
            selected_raw = selected
            audit.to_csv(output / "libraries.csv", index=False)
    pd.DataFrame(audit_rows).to_csv(output / "thresholds.csv", index=False)
    return profiles, tables_by_subset, selected_raw


def abundance(raw, tables, metadata, output, chain):
    (counts, units), _ = prepare_mouse_evidence(raw, metadata)
    edge, status = fit_edger(counts, units)
    if status != EDGER_SUCCESS:
        raise ValueError(status)
    fisher = fit_fisher(counts, units, n_workers=1)
    remaining = structural_remainder(tables)
    ranking = rank_v_segments(remaining, edge, fisher, "spleen", "Splenic pooled analysis")
    compared = edge.merge(fisher, on="feature_id", validate="one_to_one")
    ranking.to_csv(output / "abundance.csv", index=False)
    compared.to_csv(output / "edgeR_fisher.csv", index=False)
    g1 = set(tables["g1"].index)
    stats = []
    flags = {}
    for test, q, effect, method in ((edge, "FDR", "logFC", "edgeR"),
                                    (fisher, "Fisher_FDR_g1", "Fisher_log2FC_g1", "Fisher")):
        positive = test.loc[test[q].lt(.05) & test[effect].gt(0), "feature_id"]
        keys = set(positive)
        gene_names = test["feature_id"].map(lambda feature: feature[1])
        flags[method] = [
            set(gene_names[test[q].lt(.05)]), set(gene_names[test[effect].gt(0)]),
            {key[1] for key in keys},
            set(ranking.loc[ranking[f"{method.lower()}_g1_enriched" if method == "Fisher" else "edgeR_g1_enriched"], "v_gene"]),
        ]
        stats.append({"chain": chain, "method": method, "significant": int(test[q].lt(.05).sum()),
                      "positive": len(keys), "positive_in_g1": len(keys & g1),
                      "retained_positive": len(keys & remaining)})
    pd.DataFrame(stats).to_csv(output / "signal.csv", index=False)
    structural = v_plot_profile(tables)
    limits = fit_plot_thresholds(structural["frequency_in_full_g1_pct"], structural["retained_cdr3_share_pct"])
    chosen = set(structural.loc[limits.includes(structural["frequency_in_full_g1_pct"],
                                                structural["retained_cdr3_share_pct"]), "v"])
    flags["Both"] = [a & b for a, b in zip(flags["edgeR"], flags["Fisher"], strict=True)]
    fields = ("any_significant", "any_positive", "significant_positive", "positive_mean")
    cells = [{"test": test, "criterion": field, "all": len(genes), "selected": len(genes & chosen)}
             for test, sets in flags.items() for field, genes in zip(fields, sets, strict=True)]
    pd.DataFrame(cells).to_csv(output / "abundance_criteria.csv", index=False)
    association = spearmanr(compared["logFC"], compared["Fisher_log2FC_g1"]).statistic
    return len(compared), association


def selection_benchmark(raw, metadata, output, chain):
    tables, _ = downsample_group_tables(raw, depth=DEPTH, seed=SEED)
    panel = JointRecoveryPanel(pool_mice(tables, metadata))
    (counts, units), (v_counts, v_units) = prepare_mouse_evidence(raw, metadata)
    samples = metadata.groupby(["group_no", "mouse_id"])["sample_id"].agg(set)
    all_samples = set.union(*samples)
    records = []
    for idx, (group, mouse) in enumerate(panel.mice):
        removed = (group, mouse)
        training = subset_v_tables(tables, all_samples - samples.at[removed])
        features = training_predictors(training, counts, units, v_counts, v_units, (removed,))
        candidates = shortlist(features)
        scores = panel.measure(panel.reference((idx,)), (idx,))[0, 0]
        eligible = panel.genes.get_indexer(features["v_gene"])
        control = v_units.loc[~pd.MultiIndex.from_frame(v_units[["group_no", "mouse_id"]]).isin([removed])]
        v_edge, status = fit_edger(v_counts.loc[:, control.index], control)
        if status != EDGER_SUCCESS:
            raise ValueError(status)
        v_edge["v_gene"] = v_edge["feature_id"].map(lambda feature: feature[1])
        v_edge = v_edge.set_index("v_gene").reindex(features["v_gene"])
        effect = v_edge["logFC"].fillna(0).clip(lower=0)
        significant_effect = (effect * -np.log10(v_edge["FDR"].clip(lower=np.finfo(float).tiny))
                              ).where(v_edge["FDR"].lt(.05), 0)
        rankings = (
            None, features["n_g1_exclusive_aav"], features["n_full"],
            features[FEATURES[0]], features[FEATURES[1]],
            4*features[FEATURES[0]] + 2*features[FEATURES[1]] + features[FEATURES[2]],
            effect, significant_effect,
        )
        observed = scores[eligible]
        available = np.isfinite(observed)
        for method, priority in zip(METHODS, rankings, strict=True):
            weights = (features["v_gene"].isin(candidates["v_gene"]).to_numpy(dtype=float) if priority is None else
                       top_weights(priority.to_numpy(), len(candidates)))
            denominator = weights[available].sum()
            recovery = observed[available] @ weights[available] / denominator if denominator else np.nan
            records.append({"chain": chain, "method": method, "mouse": mouse, "group": group,
                            "n_selected": len(candidates), "recovery_pct": recovery})
    by_mouse = pd.DataFrame(records)
    by_mouse.to_csv(output / "selection_by_mouse.csv", index=False)
    means = by_mouse.groupby(["method", "group"])["recovery_pct"].mean().unstack()
    effect = (means["g1"] - (means["g5"] + means["g6"]) / 2).rename("contrast_pp")
    wide = by_mouse.pivot(index=["group", "mouse"], columns="method", values="recovery_pct")
    rng = np.random.default_rng(SEED)
    boot = {}
    for group in GROUPS:
        measurements = wide.loc[group].reindex(columns=effect.index).to_numpy()
        indices = rng.integers(len(measurements), size=(2000, len(measurements)))
        boot[group] = np.nanmean(measurements[indices], axis=1)
    contrasts = boot["g1"] - (boot["g5"] + boot["g6"]) / 2
    differences = contrasts - contrasts[:, effect.index.get_loc("edgeR_binary_ungated"), None]
    summary = effect.to_frame().join(pd.DataFrame(np.quantile(differences, [.025, .975], axis=0).T,
                                                   index=effect.index, columns=["delta_lower", "delta_upper"]))
    summary.to_csv(output / "selection.csv")


def rarefaction(raw, displayed, output):
    registry, libraries = {}, {}
    for group, table in raw.items():
        indices = np.array([registry.setdefault(key, len(registry)) for key in table.index])
        values = [(sample, indices[column.gt(0)], column[column.gt(0)].to_numpy())
                  for sample, column in table.items()]
        libraries[group] = values
    genes, code = np.unique([key[1] for key in registry], return_inverse=True)
    outcomes, shortlists = [], {}
    for seed in RAREFACTION_SEEDS:
        sums = {}
        for group, values in libraries.items():
            accumulated = np.zeros(len(code), dtype=np.int64)
            for sample, ids, counts in values:
                rng = np.random.default_rng(_stable_sample_seed(seed, sample))
                accumulated[ids] += rng.multivariate_hypergeometric(counts, DEPTH) if counts.sum() > DEPTH else counts
            sums[group] = accumulated
        g1 = sums["g1"]
        present = g1 > 0
        remainder = present & ~np.logical_or.reduce([sums[g] > 0 for g in SUBTRACTION_GROUPS])
        private = present & (sums["g5"] == 0) & (sums["g6"] == 0)
        full = np.bincount(code[present], minlength=len(genes))
        keep = full > 0
        remaining = np.bincount(code[remainder], minlength=len(genes))[keep]
        share = full[keep] / present.sum()
        retention = remaining / full[keep]
        cutoffs = fit_plot_thresholds(100*share, 100*retention)
        eligible = cutoffs.includes(100*share, 100*retention)
        selected = genes[keep][eligible].tolist()
        shortlists[seed] = (selected, cutoffs)
        private_umi = np.bincount(code[private], weights=g1[private], minlength=len(genes))
        total_umi = np.bincount(code[present], weights=g1[present], minlength=len(genes))
        outcomes.append(pd.DataFrame({"v": genes[keep], "seed": seed,
                                      "only_g1_share": private_umi[keep] / total_umi[keep]}))
    base = displayed.loc[displayed["selected"]].sort_values("score_s", ascending=False)["v"].tolist()
    assert set(base) == set(shortlists[DISPLAY_SEED][0])
    compare_shortlists((base, shortlists[DISPLAY_SEED][1]), shortlists).to_csv(output / "rarefaction.csv", index=False)
    matrix, summary = regional_tables(pd.concat(outcomes), base, DISPLAY_SEED)
    matrix.to_csv(output / "rarefaction_regions.csv")
    summary.to_csv(output / "rarefaction_regions_summary.csv")
    return representative_seed(shortlists)[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--chain", choices=("TRA", "TRB"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--abundance-only", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(args.cache / "mouse_metadata.csv", dtype={"mouse_id": str})
    sample_ids = metadata.loc[metadata["chain"].eq(args.chain) & metadata["_tissue"].eq("spleen"), "sample_id"]
    raw = load_snapshot_tables(args.cache / "raw", args.chain, metadata, sample_ids=sample_ids)
    if args.abundance_only:
        matched, _ = matched_internal(raw, metadata, "pooled")
        tables, _ = downsample_group_tables(matched, depth=DEPTH, seed=DISPLAY_SEED)
        abundance(matched, tables, metadata, args.output, args.chain)
        return
    profiles, tables, matched = prepared_profiles(raw, metadata, args.chain, args.output)
    del raw
    print("Profiles and cutoffs", flush=True)

    depth = depth_check(pd.concat([matched[g] for g in ("g1", "g2")], axis=1).fillna(0).astype(np.int64), _stable_sample_seed)
    depth.to_csv(args.output / "depth.csv", index=False)
    mds = mds_points(pd.concat([tables["pooled"][g] for g in ("g1", "g2")], axis=1).fillna(0), metadata)
    mds.to_csv(args.output / "mds.csv")
    pd.DataFrame({"n": [len(tables["pooled"][g]) for g in tables["pooled"]]},
                 index=tables["pooled"].keys()).to_csv(args.output / "group_sizes.csv")
    graphic = plot_group_overlap(tables["pooled"], args.chain)
    graphic._suptitle.set_visible(False)
    with plt.rc_context({"font.family": "P052", "font.size": 10.5, "pdf.fonttype": 3}):
        for label in graphic.axes[0].texts:
            label.set_fontfamily("P052")
        graphic.savefig(args.output / "overlap.pdf", bbox_inches="tight", pad_inches=.07)
    plt.close(graphic)
    print("Depth and MDS", flush=True)

    common, correlation = abundance(matched, tables["pooled"], metadata, args.output, args.chain)
    print("Abundance", common, correlation, flush=True)
    selection_benchmark(matched, metadata, args.output, args.chain)
    print("Selection benchmark", flush=True)
    center = rarefaction(matched, profiles["pooled"], args.output)
    print("Rarefaction", center, flush=True)


if __name__ == "__main__":
    main()
