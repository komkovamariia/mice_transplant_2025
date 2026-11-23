"""Add fixed aaV spike-ins to count matrices and summarize recovery."""

from __future__ import annotations

import ast
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

SUBTRACTION_GROUPS = ("g2", "g4", "g5", "g6")
DEFAULT_FRACTIONS = (0.0001, 0.001, 0.01)
ARTICLE_FRACTIONS = (
    0.0000001, 0.0000002, 0.0000005, 0.000001, 0.000002, 0.000005, 0.00001, 0.00002,
    0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.10,
)
EDGER_SUCCESS = "edgeR quasi-likelihood model fitted successfully."


def validate_case_name(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,99}", value):
        raise ValueError("Run and case names must be simple identifiers without slashes.")
    return value


def aa_v_key(value) -> tuple[str, str]:
    parsed = ast.literal_eval(value) if isinstance(value, str) else value
    if not isinstance(parsed, (tuple, list)) or len(parsed) != 2:
        raise ValueError(f"Expected an exact (CDR3 amino-acid sequence, V gene) key: {value!r}")
    return str(parsed[0]), str(parsed[1])


def validate_counts(table: pd.DataFrame) -> pd.DataFrame:
    out = table.copy().fillna(0)
    keys = [aa_v_key(value) for value in out.index]
    out.index = pd.Index(keys, tupleize_cols=False)
    if not out.index.is_unique or not out.columns.is_unique:
        raise ValueError("Count matrices require unique aaV keys and sample columns.")
    if all(pd.api.types.is_integer_dtype(dtype) for dtype in out.dtypes):
        values = out.to_numpy(copy=False)
        invalid = (values < 0).any() or (values > np.iinfo(np.int64).max).any()
    else:
        values = out.to_numpy(dtype=float)
        invalid = (not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any())
    if invalid:
        raise ValueError("Spike-in requires finite, non-negative integer UMI counts.")
    if not out.dtypes.eq(np.dtype("int64")).all():
        out = out.astype(np.int64)
    # exact-set analysis treats any retained row as present
    return out.loc[out.sum(axis=1).gt(0)]


def combine_count_tables(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Join count matrices by aaV while keeping sample columns separate."""
    checked = [validate_counts(table) for table in tables.values()]
    columns = [column for table in checked for column in table.columns]
    if len(columns) != len(set(columns)):
        raise ValueError("A sample column occurs in more than one count table.")
    return pd.concat(checked, axis=1).fillna(0).astype(np.int64)


def _metadata_digest(metadata: pd.DataFrame) -> str:
    fields = ["sample_id", "group_no", "mouse_id", "source", "subtype"]
    text = metadata[fields].astype(str).sort_values(fields).to_csv(index=False)
    return hashlib.sha256(text.encode()).hexdigest()


def snapshot_tables(tables, directory: Path, chain: str, metadata) -> None:
    folder = directory / "count_snapshots"
    folder.mkdir(parents=True, exist_ok=True)
    files = {}
    for group, table in tables.items():
        path = folder / f"{chain}_{group}.csv.gz"
        checked = validate_counts(table)
        checked.index = checked.index.map(repr)
        checked.to_csv(path, index_label="aaV")
        files[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {"metadata_sha256": _metadata_digest(metadata), "files": files}
    (folder / f"{chain}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def load_snapshot_tables(directory: Path, chain: str, metadata=None, *, sample_ids=None) -> dict:
    folder = directory / "count_snapshots"
    manifest = json.loads((folder / f"{chain}_manifest.json").read_text())
    if metadata is not None and _metadata_digest(metadata) != manifest["metadata_sha256"]:
        raise ValueError("Sample metadata changed after the spike-in baseline was frozen.")
    requested = set(sample_ids) if sample_ids is not None else None
    tables = {}
    for filename, expected_digest in manifest["files"].items():
        path = folder / filename
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected_digest:
            raise ValueError(f"Baseline count snapshot changed: {path}")
        group = filename.removeprefix(f"{chain}_").removesuffix(".csv.gz")
        if requested is None:
            frame = pd.read_csv(path, index_col=0)
        else:
            chunks = pd.read_csv(path, index_col=0, chunksize=50_000,
                                 usecols=lambda name: name == "aaV" or name in requested)
            frame = pd.concat([part.loc[part.ne(0).any(axis=1)] for part in chunks])
        tables[group] = validate_counts(frame)
    return tables


def select_targets(
    tables: dict,
    retention: pd.DataFrame,
    ranking: pd.DataFrame,
    *,
    v_gene: str | None = None,
    chain: str = "TRA",
    n_clones: int = 10,
    max_g1_fraction: float = 1e-5,
    seed: int = 1031,
) -> dict:
    """Select observed, control-absent aaV from a baseline noncandidate V segment."""
    if chain not in {"TRA", "TRB"}:
        raise ValueError("The spike-in chain must be TRA or TRB.")
    v_prefix = chain + "V"
    if n_clones < 1 or not np.isfinite(max_g1_fraction) or not 0 <= max_g1_fraction < 1:
        raise ValueError("Use n_clones >= 1 and a finite rarity fraction in [0, 1).")
    checked = {group: validate_counts(table) for group, table in tables.items()}
    g1 = checked["g1"]
    libraries = g1.sum(axis=0)
    if g1.empty or libraries.le(0).any():
        raise ValueError("Every g1 sample needs a nonzero baseline UMI library.")
    forbidden = set().union(
        *(set(checked[group].index) for group in SUBTRACTION_GROUPS)
    )
    observed = set().union(*(set(table.index) for table in checked.values()))
    frequencies = g1.div(libraries, axis=1).max(axis=1)

    graphical = set(retention.loc[
        retention["frequency_in_full_g1"].ge(0.006)
        & retention["retained_cdr3_share"].ge(0.66), "v"
    ])
    inferential = set(ranking.loc[
        ranking["edgeR_g1_enriched"].eq(True) | ranking["fisher_g1_enriched"].eq(True),
        "v_gene",
    ])
    excluded_genes = graphical | inferential
    eligible = {}
    for key in sorted(observed - forbidden):
        if not key[1].startswith(v_prefix) or key[1] in excluded_genes:
            continue
        if float(frequencies.get(key, 0)) <= max_g1_fraction:
            eligible.setdefault(key[1], []).append(key)
    enough = {gene: keys for gene, keys in eligible.items() if len(keys) >= n_clones}
    if v_gene is not None and v_gene not in enough:
        raise ValueError(
            f"{v_gene} is a baseline candidate or has fewer than {n_clones} eligible aaV. "
            "Inspect selection criteria; controls and rarity thresholds are never relaxed automatically."
        )
    if not enough:
        raise ValueError(
            f"No noncandidate {v_prefix} has {n_clones} observed, rare-in-g1 aaV absent from "
            "g2/g4/g5/g6. Try an explicitly smaller --spike-clones value or review "
            "--spike-max-g1-fraction. No synthetic sequences were invented."
        )
    rng = np.random.default_rng(seed)
    chosen_gene = v_gene or str(rng.choice(sorted(enough)))
    candidates = enough[chosen_gene]
    # prefer g1-absent clonotypes so the spike can change exact-set retention
    absent = [key for key in candidates if key not in g1.index]
    rare = [key for key in candidates if key in g1.index]
    selected = []
    for pool in (absent, rare):
        if pool:
            order = rng.permutation(len(pool))
            selected.extend(
                pool[position]
                for position in order[: n_clones - len(selected)]
            )
        if len(selected) == n_clones:
            break
    return {
        "v_gene": chosen_gene,
        "chain": chain,
        "clonotypes": [list(key) for key in selected],
        "n_clones": n_clones,
        "seed": seed,
        "max_g1_fraction": max_g1_fraction,
        "baseline_absent_in_g1": sum(key not in g1.index for key in selected),
        "baseline_rare_in_g1": sum(key in g1.index for key in selected),
        "eligible_genes": {gene: len(keys) for gene, keys in sorted(enough.items())},
        "excluded_graphical_candidates": sorted(graphical),
        "excluded_inferential_candidates": sorted(inferential),
        "relatedness": f"same exact {v_prefix} annotation; observed CDR3aa sequences",
        "dose_definition": "total added family UMI / original UMI library in each g1 sample",
        "allocation": "equal shares with deterministic largest-remainder rounding",
    }


def apply_spike_to_tables(tables: dict, selection: dict, fraction: float):
    if not np.isfinite(fraction) or not 0 <= fraction < 1:
        raise ValueError("Added UMI fraction must be finite and in [0, 1).")
    result = {group: validate_counts(table) for group, table in tables.items()}
    keys = [aa_v_key(key) for key in selection["clonotypes"]]
    if not keys or len(keys) != len(set(keys)):
        raise ValueError("Select at least one unique aaV clonotype.")
    if any(key[1] != selection["v_gene"] for key in keys):
        raise ValueError("Every selected aaV must use the selected V segment.")
    controls = set().union(
        *(set(result[group].index) for group in SUBTRACTION_GROUPS)
    )
    if set(keys) & controls:
        raise ValueError("Selected aaV occur in a subtraction control.")
    g1 = result["g1"]
    libraries = g1.sum(axis=0)
    if not len(libraries) or libraries.le(0).any():
        raise ValueError("Every g1 sample needs a nonzero baseline UMI library.")
    added_keys = [key for key in keys if key not in g1.index]
    g1 = g1.reindex(pd.Index(list(g1.index) + added_keys, tupleize_cols=False), fill_value=0)
    doses = []
    for sample, library in libraries.items():
        budget = int(np.floor(float(library) * fraction + 0.5))
        quotient, remainder = divmod(budget, len(keys))
        for position, key in enumerate(keys):
            amount = quotient + int(position < remainder)
            g1.at[key, sample] += amount
            doses.append({
                "sample_id": sample, "cdr3aa": key[0], "v_gene": key[1],
                "baseline_library_umi": int(library), "requested_fraction": fraction,
                "added_umi": amount, "sample_added_umi": budget,
                "realized_added_fraction": budget / int(library),
                "realized_added_fraction_post_spike": budget / (int(library) + budget),
            })
    result["g1"] = validate_counts(g1)
    return result, pd.DataFrame(doses)


def measure_case(selection, fraction, g1, remaining, retention, edger, ranking):
    keys = [aa_v_key(key) for key in selection["clonotypes"]]
    g1 = validate_counts(g1)
    key_set = set(keys)
    present = key_set & set(g1.index)
    recovered = key_set & {aa_v_key(key) for key in remaining}
    tested = edger.copy()
    tested["aaV"] = tested["feature_id"].map(aa_v_key)
    tested = tested[tested["aaV"].isin(key_set)]
    positive = tested[tested["logFC"].gt(0) & tested["FDR"].lt(0.05)]
    gene = selection["v_gene"]
    vrow = retention.loc[retention["v"].eq(gene)]
    rrow = ranking.loc[ranking["v_gene"].eq(gene)]
    def retention_value(column):
        return float(vrow.iloc[0][column]) if len(vrow) else 0.0
    gene_umi = float(g1.loc[[key for key in g1.index if key[1] == gene]].to_numpy().sum())
    metrics = {
        "fraction": fraction, "v_gene": gene, "n_selected": len(keys),
        "n_present_g1": len(present), "n_recovered_remainder": len(recovered),
        "recovery_fraction": len(recovered) / len(keys),
        "all_selected_in_remainder": len(recovered) == len(keys),
        "n_edger_tested_selected": len(tested),
        "n_edger_positive_fdr_selected": len(positive),
        "edger_positive_recovery_fraction": len(positive) / len(keys),
        "v_unique_g1": int(retention_value("cdr3_in_full_g1")),
        "v_remaining_g1": int(retention_value("cdr3_remaining_in_g1")),
        "v_retention": retention_value("retained_cdr3_share"),
        "v_unique_share_g1": retention_value("frequency_in_full_g1"),
        "v_umi_share_g1": gene_umi / float(g1.to_numpy().sum()),
        "v_article_score": retention_value("frequency_in_full_g1") * retention_value("retained_cdr3_share"),
        "v_final_rank": int(rrow.iloc[0]["stratum_rank"]) if len(rrow) else None,
    }
    rows = []
    lookup = {
        aa_v: row
        for aa_v, (_, row) in zip(tested["aaV"], tested.iterrows())
    }
    for key in keys:
        row = lookup.get(key)
        rows.append({
            "cdr3aa": key[0], "v_gene": key[1], "present_g1": key in present,
            "in_g1_remainder": key in recovered, "edger_tested": row is not None,
            "logFC": float(row["logFC"]) if row is not None else None,
            "FDR": float(row["FDR"]) if row is not None else None,
            "positive_fdr": key in set(positive["aaV"]),
        })
    return metrics, pd.DataFrame(rows)


def write_case_metrics(directory, selection, fraction, g1, remaining, retention, edger, ranking):
    metrics, recovery = measure_case(selection, fraction, g1, remaining, retention, edger, ranking)
    (directory / "spike_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    recovery.to_csv(directory / "spike_clonotype_recovery.csv", index=False)
    return metrics


def write_baseline_metrics(directory: Path, selection: dict):
    chain = selection.get("chain", "TRA")
    tables = load_snapshot_tables(directory, chain)
    remainder = pd.read_csv(directory / "spike_baseline_remainder.csv", dtype=str, keep_default_na=False)
    edger_name = "edger_aav_results.csv" if chain == "TRA" else "trb_edger_aav_results.csv"
    return write_case_metrics(
        directory, selection, 0.0, tables["g1"], set(remainder.itertuples(index=False, name=None)),
        pd.read_csv(directory / f"{chain.lower()}_v_retention.csv"),
        pd.read_csv(directory / edger_name),
        pd.read_csv(directory / f"{chain.lower()}v_ranking.csv"),
    )


def sensitivity_summary(frame: pd.DataFrame) -> dict:
    """Summarize pass/fail transitions across the tested dose grid."""
    doses = frame.loc[frame["fraction"].gt(0)].sort_values("fraction")
    observations = []
    for metric in ("all_selected_in_remainder", "any_edger_positive_fdr", "all_edger_positive_fdr"):
        values = doses[metric].astype(bool).to_numpy()
        fractions = doses["fraction"].to_numpy(dtype=float)
        successes = fractions[values]
        first = float(successes[0]) if len(successes) else None
        below = fractions[(fractions < first) & ~values] if first is not None else []
        transitions = [
            {
                "from_fraction": float(fractions[position - 1]),
                "to_fraction": float(fractions[position]),
                "from_pass": bool(values[position - 1]),
                "to_pass": bool(values[position]),
            }
            for position in range(1, len(values))
            if values[position] != values[position - 1]
        ]
        observations.append({
            "criterion": metric, "lowest_tested_passing_fraction": first,
            "largest_tested_failure_below_first_pass": float(max(below)) if len(below) else None,
            "all_tested_doses_pass": bool(values.all()) if len(values) else False,
            "no_tested_dose_passes": not bool(values.any()),
            "observed_transitions": transitions,
            "upper_failure_observed": any(
                transition["from_pass"] and not transition["to_pass"]
                for transition in transitions
            ),
        })
    return {
        "tested_fractions": doses["fraction"].tolist(), "criteria": observations,
        "interpretation": (
            "These are observed brackets for one fixed family and uniform sample dosing. "
            "An upper sensitivity limit is not established unless recovery is lost at a "
            "higher tested dose. UMI rounding, filtering and biological replication affect detection."
        ),
    }


def summarize_experiment(
    result_root: Path, figure_root: Path, cases: list[str], *, chain: str = "TRA",
) -> pd.DataFrame:
    rows = [json.loads((result_root / case / "spike_metrics.json").read_text()) for case in cases]
    frame = pd.DataFrame(rows).sort_values("fraction").reset_index(drop=True)
    frame["v_final_rank"] = pd.to_numeric(frame["v_final_rank"], errors="coerce")
    baseline = frame.loc[frame["fraction"].eq(0)].iloc[0]
    for column in ("v_remaining_g1", "v_retention", "v_unique_share_g1", "v_umi_share_g1", "v_article_score"):
        frame[f"delta_{column}"] = frame[column] - baseline[column]
    frame["rank_improvement"] = baseline["v_final_rank"] - frame["v_final_rank"]
    frame["entered_ranking"] = pd.isna(baseline["v_final_rank"]) & frame["v_final_rank"].notna()
    frame["remainder_increased"] = frame["delta_v_remaining_g1"].gt(0)
    frame["retention_increased"] = frame["delta_v_retention"].gt(0)
    frame["unique_share_increased"] = frame["delta_v_unique_share_g1"].gt(0)
    frame["umi_share_increased"] = frame["delta_v_umi_share_g1"].gt(0)
    frame["rank_improved"] = frame["rank_improvement"].gt(0) | frame["entered_ranking"]
    frame["any_edger_positive_fdr"] = frame["n_edger_positive_fdr_selected"].gt(0)
    frame["all_edger_positive_fdr"] = frame["n_edger_positive_fdr_selected"].eq(frame["n_selected"])
    frame.to_csv(result_root / "spike_in_summary.csv", index=False)
    bounds = sensitivity_summary(frame)
    (result_root / "sensitivity_bounds.json").write_text(json.dumps(bounds, indent=2) + "\n")

    import matplotlib.pyplot as plt
    from src.article_figures import (
        ARTICLE_FONTS, article_dpi, article_name, panel_bbox, plot_spike_in,
    )

    output = figure_root / chain
    output.mkdir(parents=True, exist_ok=True)
    filename = article_name(chain, "spike_in")
    surplus = output / "surplus"
    surplus.mkdir(exist_ok=True)
    with plt.rc_context({"font.family": "sans-serif", "font.sans-serif": ARTICLE_FONTS}):
        fig = plot_spike_in(frame, chain)
        fig.savefig(output / filename, dpi=article_dpi(filename))
        for ax, panel in zip(fig.axes, ("exact_recovery", "edger_recovery", "umi_share", "final_rank"), strict=True):
            fig.savefig(surplus / f"spike_{panel}.png", dpi=300, bbox_inches=panel_bbox(ax))
        plt.close(fig)

    report = [
        "# Spike-in sensitivity experiment",
        "",
        (
            f"Selected V segment: {frame.iloc[0]['v_gene']}. "
            f"Family size: {int(frame.iloc[0]['n_selected'])} aaV."
        ),
        "",
        bounds["interpretation"],
        "",
        "The CSV reports each hypothesis separately.",
        "Exact-set metrics count unique aaV and can plateau after full recovery.",
        (
            "See spike_in_summary.csv, sensitivity_bounds.json and each case's "
            "spike_clonotype_recovery.csv."
        ),
        "",
    ]
    (result_root / "conclusion.md").write_text("\n".join(report))
    return frame

