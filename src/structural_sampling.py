"""Deterministic per-library rarefaction for structural aaV analysis."""

import hashlib

import numpy as np
import pandas as pd


def _stable_sample_seed(base_seed, sample_id):
    digest = hashlib.sha256(f"{base_seed}|{sample_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "little") % (2**32 - 1)


def downsample_count_table(
    table,
    *,
    depth,
    seed,
    eligible_samples=None,
):
    numeric = table.apply(pd.to_numeric, errors="raise").fillna(0).astype(np.int64)
    eligible = None if eligible_samples is None else {str(x) for x in eligible_samples}
    sampled_columns = {}
    audit_rows = []
    for sample_id in numeric.columns:
        counts = numeric[sample_id].to_numpy(dtype=np.int64)
        total = int(counts.sum())
        kept = total >= depth and (eligible is None or str(sample_id) in eligible)
        audit_rows.append(
            {
                "sample_id": str(sample_id),
                "functional_umi_before": total,
                "target_umi": int(depth),
                "kept_for_structural_analysis": bool(kept),
                "functional_umi_after": int(depth if kept else 0),
            }
        )
        if not kept:
            continue

        positive_positions = np.flatnonzero(counts)
        positive_counts = counts[positive_positions]
        if total == depth:
            sampled_positive = positive_counts.copy()
        else:
            rng = np.random.default_rng(_stable_sample_seed(seed, sample_id))
            sampled_positive = rng.multivariate_hypergeometric(positive_counts, depth)
        retained = sampled_positive > 0
        sampled_columns[str(sample_id)] = pd.Series(
            sampled_positive[retained].astype(np.int64),
            index=numeric.index[positive_positions[retained]],
            dtype=np.int64,
        )

    if not sampled_columns:
        raise ValueError(
            f"No sample reached the structural downsampling depth of {depth:,} functional UMI."
        )
    sampled_table = pd.concat(sampled_columns, axis=1).fillna(0).astype(np.int64)
    sampled_table.index = pd.Index(sampled_table.index, tupleize_cols=False)
    return sampled_table, pd.DataFrame(audit_rows)


def downsample_group_tables(
    tables,
    *,
    depth,
    seed,
    eligible_by_group=None,
):
    sampled = {}
    audits = []
    for group, table in tables.items():
        eligible = None if eligible_by_group is None else eligible_by_group.get(group, set())
        sampled[group], audit = downsample_count_table(
            table,
            depth=depth,
            seed=seed,
            eligible_samples=eligible,
        )
        audit.insert(0, "group_no", group)
        audits.append(audit)
    return sampled, pd.concat(audits, ignore_index=True)
