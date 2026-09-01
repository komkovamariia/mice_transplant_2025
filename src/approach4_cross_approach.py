"""Direction-aware comparison of V-segment rankings across methods and strata."""

from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .figures import save_figure
from .reporting import format_gene_list, save_conclusion
from .strata import STRATUM_LABELS, repository_root

APPROACH = "04_cross_approach"
SOURCES = {
    "set_count": {
        "folder": "01_set_count",
        "effect_column": "mean_edger_log2fc",
        "rank_column": "stratum_rank",
        "relative_path": "{stratum}/trav_ranking.csv",
    },
    "sequence_embedding": {
        "folder": "02_sequence_embedding",
        "effect_column": "median_effect_g1_vs_allogeneic",
        "rank_column": "rank",
        "ranking_suffix": "v_gene_ranking",
    },
    "clone_alloreactivity": {
        "folder": "03_clone_alloreactivity",
        "effect_column": "mean_effect_g1_vs_allogeneic",
        "rank_column": "rank",
        "ranking_suffix": "v_gene_ranking",
    },
}


def _output_dir() -> Path:
    path = repository_root() / "results" / APPROACH
    path.mkdir(parents=True, exist_ok=True)
    return path


def _ranking_path(folder: str, stratum: str, suffix: str) -> Path:
    return repository_root() / "results" / folder / f"{stratum}_{suffix}.csv"


def load_rankings(stratum: str) -> dict[str, pd.DataFrame]:
    rankings = {}
    missing = []
    for name, specification in SOURCES.items():
        if "relative_path" in specification:
            path = (
                repository_root()
                / "results"
                / specification["folder"]
                / specification["relative_path"].format(stratum=stratum)
            )
        else:
            path = _ranking_path(
                specification["folder"],
                stratum,
                specification["ranking_suffix"],
            )
        if not path.exists():
            missing.append(str(path))
            continue

        frame = pd.read_csv(path)
        effect_column = specification["effect_column"]
        rank_column = specification["rank_column"]
        required = {"v_gene", rank_column, effect_column}
        if not required.issubset(frame.columns):
            raise ValueError(f"{path} must contain {sorted(required)}.")
        rankings[name] = (
            frame[["v_gene", rank_column, effect_column]]
            .rename(columns={rank_column: "rank"})
            .dropna(subset=["v_gene", "rank"])
            .drop_duplicates("v_gene")
            .rename(columns={effect_column: (f"{name}_effect_g1_vs_allogeneic")})
        )

    if missing:
        raise FileNotFoundError(
            "Cross-approach comparison requires all three upstream "
            "rankings. Missing: " + "; ".join(missing)
        )
    return rankings


def _direction_agreement(
    effect_a: pd.Series,
    effect_b: pd.Series,
) -> tuple[int, float]:
    valid = effect_a.notna() & effect_b.notna()
    signs_a = np.sign(effect_a[valid].to_numpy(float))
    signs_b = np.sign(effect_b[valid].to_numpy(float))
    nonzero = (signs_a != 0) & (signs_b != 0)
    if not nonzero.any():
        return 0, np.nan
    return int(nonzero.sum()), float(np.mean(signs_a[nonzero] == signs_b[nonzero]))


def _pairwise_metrics(
    rankings: dict[str, pd.DataFrame],
    stratum: str,
) -> pd.DataFrame:
    rows = []
    for approach_a, approach_b in combinations(rankings, 2):
        merged = rankings[approach_a].merge(
            rankings[approach_b],
            on="v_gene",
            suffixes=("_a", "_b"),
        )
        if len(merged) >= 3:
            rank_rho, rank_p = spearmanr(
                merged["rank_a"],
                merged["rank_b"],
            )
        else:
            rank_rho, rank_p = np.nan, np.nan

        effect_a = merged[f"{approach_a}_effect_g1_vs_allogeneic"]
        effect_b = merged[f"{approach_b}_effect_g1_vs_allogeneic"]
        valid_effect = effect_a.notna() & effect_b.notna()
        if valid_effect.sum() >= 3:
            effect_rho, effect_p = spearmanr(
                effect_a[valid_effect],
                effect_b[valid_effect],
            )
        else:
            effect_rho, effect_p = np.nan, np.nan
        n_directional, directional_agreement = _direction_agreement(
            effect_a,
            effect_b,
        )

        top10_a = set(rankings[approach_a].nsmallest(10, "rank")["v_gene"])
        top10_b = set(rankings[approach_b].nsmallest(10, "rank")["v_gene"])
        rows.append(
            {
                "stratum": stratum,
                "approach_a": approach_a,
                "approach_b": approach_b,
                "n_shared_v_genes": len(merged),
                "rank_spearman_rho": rank_rho,
                "rank_spearman_p": rank_p,
                "effect_spearman_rho": effect_rho,
                "effect_spearman_p": effect_p,
                "n_directionally_comparable": n_directional,
                "direction_agreement_fraction": directional_agreement,
                "top10_overlap": len(top10_a & top10_b),
            }
        )
    return pd.DataFrame(rows)


def _consensus(
    rankings: dict[str, pd.DataFrame],
    stratum: str,
) -> pd.DataFrame:
    all_v_genes = sorted(
        set().union(*(set(frame["v_gene"]) for frame in rankings.values()))
    )
    table = pd.DataFrame({"v_gene": all_v_genes})

    for name, frame in rankings.items():
        maximum_rank = max(float(frame["rank"].max()), 1.0)
        renamed = frame.rename(columns={"rank": f"{name}_rank"}).copy()
        renamed[f"{name}_rank_fraction"] = renamed[f"{name}_rank"] / maximum_rank
        table = table.merge(
            renamed,
            on="v_gene",
            how="left",
        )

    rank_fraction_columns = [
        column for column in table if column.endswith("_rank_fraction")
    ]
    effect_columns = [
        column for column in table if column.endswith("_effect_g1_vs_allogeneic")
    ]
    table["n_approaches"] = table[rank_fraction_columns].notna().sum(axis=1)
    table["mean_rank_fraction"] = table[rank_fraction_columns].mean(axis=1, skipna=True)
    table["n_directional_effects"] = table[effect_columns].notna().sum(axis=1)
    table["n_g1_enriched_effects"] = table[effect_columns].gt(0).sum(axis=1)
    table["n_allogeneic_enriched_effects"] = table[effect_columns].lt(0).sum(axis=1)
    directional_maximum = table[
        ["n_g1_enriched_effects", "n_allogeneic_enriched_effects"]
    ].max(axis=1)
    table["direction_agreement_fraction"] = directional_maximum / table[
        "n_directional_effects"
    ].replace(0, np.nan)

    table["consensus_direction"] = "undetermined"
    table.loc[
        table["n_g1_enriched_effects"] > table["n_allogeneic_enriched_effects"],
        "consensus_direction",
    ] = "g1_enriched"
    table.loc[
        table["n_allogeneic_enriched_effects"] > table["n_g1_enriched_effects"],
        "consensus_direction",
    ] = "allogeneic_enriched"
    table.loc[
        table["n_g1_enriched_effects"].eq(table["n_allogeneic_enriched_effects"])
        & table["n_directional_effects"].gt(0),
        "consensus_direction",
    ] = "mixed"

    table = table.sort_values(
        [
            "n_approaches",
            "direction_agreement_fraction",
            "mean_rank_fraction",
        ],
        ascending=[False, False, True],
        na_position="last",
    ).reset_index(drop=True)
    table["consensus_rank"] = np.arange(1, len(table) + 1)
    table["stratum"] = stratum
    return table


def _plot_consensus(
    consensus: pd.DataFrame,
    stratum: str,
) -> None:
    top = (
        consensus[consensus["n_approaches"].eq(3)]
        .head(15)
        .sort_values("mean_rank_fraction", ascending=False)
    )
    if top.empty:
        return

    colors = top["consensus_direction"].map(
        {
            "g1_enriched": "#3B6FB6",
            "allogeneic_enriched": "#C65D3A",
            "mixed": "#8A8A8A",
            "undetermined": "#B8B8B8",
        }
    )
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(
        top["v_gene"],
        1 - top["mean_rank_fraction"],
        color=colors,
    )
    ax.set_xlabel("Cross-approach rank concordance")
    ax.set_ylabel("V segment")
    ax.set_title(f"Direction-aware consensus ranking: {STRATUM_LABELS[stratum]}")
    save_figure(
        fig,
        APPROACH,
        stratum,
        "consensus_v_segment_ranking",
    )


def run_stratum(stratum: str) -> dict:
    rankings = load_rankings(stratum)
    metrics = _pairwise_metrics(rankings, stratum)
    consensus = _consensus(rankings, stratum)
    output = _output_dir()
    metrics.to_csv(
        output / f"{stratum}_pairwise_metrics.csv",
        index=False,
    )
    consensus.to_csv(
        output / f"{stratum}_consensus_ranking.csv",
        index=False,
    )
    _plot_consensus(consensus, stratum)

    top = consensus[consensus["n_approaches"].eq(3)].head(5)
    top_genes = top["v_gene"].tolist()
    payload = {
        "stratum": stratum,
        "label": STRATUM_LABELS[stratum],
        "top_consensus_v_genes": top_genes,
        "top_consensus_directions": dict(
            zip(top["v_gene"], top["consensus_direction"])
        ),
        "pairwise_rank_spearman": {
            f"{row.approach_a}_vs_{row.approach_b}": (
                None if pd.isna(row.rank_spearman_rho) else float(row.rank_spearman_rho)
            )
            for row in metrics.itertuples()
        },
        "pairwise_top10_overlap": {
            f"{row.approach_a}_vs_{row.approach_b}": int(row.top10_overlap)
            for row in metrics.itertuples()
        },
        "pairwise_direction_agreement": {
            f"{row.approach_a}_vs_{row.approach_b}": (
                None
                if pd.isna(row.direction_agreement_fraction)
                else float(row.direction_agreement_fraction)
            )
            for row in metrics.itertuples()
        },
    }
    return save_conclusion(
        output,
        stratum,
        f"Cross-approach conclusion: {STRATUM_LABELS[stratum]}",
        payload,
        [
            (
                "All three approaches produced a V-segment ranking for "
                f"{STRATUM_LABELS[stratum]}."
            ),
            (f"The leading cross-method segments were {format_gene_list(top_genes)}."),
            (
                "Rank overlap and signed g1-versus-allogeneic effects are "
                "reported separately, so a shared segment with opposing "
                "directions remains visible rather than being forced into "
                "agreement."
            ),
        ],
    )


def compile_summary(
    strata: list[str] | tuple[str, ...],
) -> pd.DataFrame:
    output = _output_dir()
    rows = []
    for stratum in strata:
        path = output / f"{stratum}_conclusion.json"
        if path.exists():
            rows.append(json.loads(path.read_text(encoding="utf-8")))

    summary = pd.DataFrame(rows)
    summary.to_json(
        output / "cross_stratum_conclusions.json",
        orient="records",
        indent=2,
    )

    top_rows = []
    for stratum in strata:
        path = output / f"{stratum}_consensus_ranking.csv"
        if not path.exists():
            continue
        frame = pd.read_csv(path)
        frame = frame[frame["n_approaches"].eq(3)].head(5)
        top_rows.extend(
            {
                "stratum": stratum,
                "v_gene": row.v_gene,
                "consensus_direction": row.consensus_direction,
            }
            for row in frame.itertuples()
        )

    frequency = pd.DataFrame(top_rows)
    if frequency.empty:
        overall = pd.DataFrame(
            columns=[
                "v_gene",
                "n_strata_top5",
                "n_g1_enriched_strata",
                "n_allogeneic_enriched_strata",
                "n_mixed_strata",
            ]
        )
    else:
        overall = (
            frequency.assign(
                g1_enriched=lambda frame: frame["consensus_direction"].eq(
                    "g1_enriched"
                ),
                allogeneic_enriched=lambda frame: frame["consensus_direction"].eq(
                    "allogeneic_enriched"
                ),
                mixed=lambda frame: frame["consensus_direction"].eq("mixed"),
            )
            .groupby("v_gene", as_index=False)
            .agg(
                n_strata_top5=("stratum", "nunique"),
                n_g1_enriched_strata=("g1_enriched", "sum"),
                n_allogeneic_enriched_strata=(
                    "allogeneic_enriched",
                    "sum",
                ),
                n_mixed_strata=("mixed", "sum"),
            )
            .sort_values(
                ["n_strata_top5", "v_gene"],
                ascending=[False, True],
            )
        )
    overall.to_csv(
        output / "overall_consensus_v_genes.csv",
        index=False,
    )

    if not overall.empty:
        top = overall.head(15).sort_values("n_strata_top5")
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.barh(
            top["v_gene"],
            top["n_strata_top5"],
            color="#4F6D7A",
        )
        ax.set_xlabel("Strata with top-five consensus support")
        ax.set_ylabel("V segment")
        ax.set_title("Cross-stratum support across all three analytical approaches")
        save_figure(
            fig,
            APPROACH,
            None,
            "overall_cross_stratum_consensus",
        )

    leading = overall.head(10)["v_gene"].tolist()
    markdown = (
        "# Cross-stratum conclusion\n\n"
        f"The V segments most frequently supported across the selected "
        f"biological strata were {format_gene_list(leading, limit=10)}. "
        "The accompanying table retains the consensus direction in each "
        "stratum, allowing compartment-specific and directionally mixed "
        "signals to be distinguished.\n"
    )
    (output / "cross_stratum_conclusion.md").write_text(
        markdown,
        encoding="utf-8",
    )
    return summary
