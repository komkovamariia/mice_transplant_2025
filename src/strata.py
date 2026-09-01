"""Data loading and biological stratification shared by all analyses."""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd

STRATA = (
    "all_combined",
    "cd4_thymus",
    "cd8_thymus",
    "cd4_spleen",
    "cd8_spleen",
    "cd4_combined",
    "cd8_combined",
    "thymus_combined",
    "spleen_combined",
)

STRATUM_LABELS = {
    "all_combined": "CD4 + CD8 T cells, thymus + spleen pooled within mouse",
    "cd4_thymus": "CD4 T cells, thymus",
    "cd8_thymus": "CD8 T cells, thymus",
    "cd4_spleen": "CD4 T cells, spleen",
    "cd8_spleen": "CD8 T cells, spleen",
    "cd4_combined": "CD4 T cells, mouse-level thymus + spleen pool",
    "cd8_combined": "CD8 T cells, mouse-level thymus + spleen pool",
    "thymus_combined": "Thymus, mouse-level CD4 + CD8 pool",
    "spleen_combined": "Spleen, mouse-level CD4 + CD8 pool",
}

STRATUM_DEFINITIONS = {
    "all_combined": {
        "cell_subsets": ("cd4", "cd8"),
        "tissues": ("thymus", "spleen"),
        "pool_within_mouse": True,
    },
    "cd4_thymus": {
        "cell_subsets": ("cd4",),
        "tissues": ("thymus",),
        "pool_within_mouse": False,
    },
    "cd8_thymus": {
        "cell_subsets": ("cd8",),
        "tissues": ("thymus",),
        "pool_within_mouse": False,
    },
    "cd4_spleen": {
        "cell_subsets": ("cd4",),
        "tissues": ("spleen",),
        "pool_within_mouse": False,
    },
    "cd8_spleen": {
        "cell_subsets": ("cd8",),
        "tissues": ("spleen",),
        "pool_within_mouse": False,
    },
    "cd4_combined": {
        "cell_subsets": ("cd4",),
        "tissues": ("thymus", "spleen"),
        "pool_within_mouse": True,
    },
    "cd8_combined": {
        "cell_subsets": ("cd8",),
        "tissues": ("thymus", "spleen"),
        "pool_within_mouse": True,
    },
    "thymus_combined": {
        "cell_subsets": ("cd4", "cd8"),
        "tissues": ("thymus",),
        "pool_within_mouse": True,
    },
    "spleen_combined": {
        "cell_subsets": ("cd4", "cd8"),
        "tissues": ("spleen",),
        "pool_within_mouse": True,
    },
}

REQUIRED_COLUMNS = {
    "cdr3",
    "v_gene",
    "umi",
    "group",
    "sample_id",
    "mouse_id",
    "source",
    "subtype",
}

_CANONICAL_AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def repository_root() -> Path:
    """Return the repository root, including during notebook execution."""
    configured = os.environ.get("MICE_TCR_REPO")
    if configured:
        return Path(configured).expanduser().resolve()

    here = Path.cwd().resolve()
    for candidate in (here, *here.parents):
        if (candidate / "environment.yml").exists() and (
            candidate / "notebooks"
        ).exists():
            return candidate
    return here


def data_dir() -> Path:
    return (
        Path(os.environ.get("MICE_TCR_DATA_DIR", repository_root() / "data"))
        .expanduser()
        .resolve()
    )


def clean_repertoire_path() -> Path:
    explicit = os.environ.get("MICE_TCR_CLEAN_PARQUET")
    if explicit:
        return Path(explicit).expanduser().resolve()
    return data_dir() / "clean_clonotypes_aaV.parquet"


def _classify(
    primary: pd.Series,
    fallback: pd.Series,
    patterns: dict[str, str],
) -> pd.Series:
    """Classify from the explicit field first and use the sample name only as fallback."""
    primary = primary.fillna("").astype(str).str.lower()
    fallback = fallback.fillna("").astype(str).str.lower()
    result = pd.Series(pd.NA, index=primary.index, dtype="object")

    primary_hits = {
        label: primary.str.contains(pattern, regex=True)
        for label, pattern in patterns.items()
    }
    primary_count = sum(primary_hits.values())
    for label, hit in primary_hits.items():
        result.loc[hit & primary_count.eq(1)] = label

    use_fallback = primary_count.eq(0)
    fallback_hits = {
        label: fallback.str.contains(pattern, regex=True)
        for label, pattern in patterns.items()
    }
    fallback_count = sum(fallback_hits.values())
    for label, hit in fallback_hits.items():
        result.loc[use_fallback & hit & fallback_count.eq(1)] = label
    return result


def _infer_cell_subset(df: pd.DataFrame) -> pd.Series:
    return _classify(
        df["subtype"],
        df["sample_id"],
        {
            "cd4": r"(?:^|[^a-z0-9])cd4(?:[^a-z0-9]|$)",
            "cd8": r"(?:^|[^a-z0-9])cd8(?:[^a-z0-9]|$)",
        },
    )


def _infer_tissue(df: pd.DataFrame) -> pd.Series:
    return _classify(
        df["source"],
        df["sample_id"],
        {
            "spleen": r"spleen",
            "thymus": r"thym|allothymus",
        },
    )


def _validate_sample_assignments(df: pd.DataFrame) -> None:
    unresolved = df[df["_cell_subset"].isna() | df["_tissue"].isna()][
        "sample_id"
    ].unique()
    if len(unresolved):
        preview = ", ".join(map(str, unresolved[:8]))
        raise ValueError(
            "Every sample must map unambiguously to CD4/CD8 and thymus/spleen. "
            f"Unresolved sample(s): {preview}"
        )

    fields = ["mouse_id", "group", "_cell_subset", "_tissue"]
    inconsistent = []
    for field in fields:
        counts = df.groupby("sample_id", observed=True)[field].nunique(dropna=False)
        inconsistent.extend(
            (sample_id, field) for sample_id in counts[counts > 1].index
        )
    if inconsistent:
        preview = ", ".join(f"{sample}:{field}" for sample, field in inconsistent[:8])
        raise ValueError(f"Inconsistent sample metadata detected: {preview}")


def normalize_repertoire(
    df: pd.DataFrame,
    *,
    source_label: str = "repertoire table",
) -> pd.DataFrame:
    """Validate and annotate a sample-resolved aaV repertoire table."""
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{source_label} is missing required columns: {missing}")

    out = df.copy()
    out["cdr3"] = out["cdr3"].fillna("").astype(str).str.strip().str.upper()
    out["v_gene"] = out["v_gene"].fillna("").astype(str).str.strip()
    numeric_umi = pd.to_numeric(out["umi"], errors="coerce")
    invalid_umi = numeric_umi.isna() | numeric_umi.lt(0)
    if invalid_umi.any():
        examples = out.loc[invalid_umi, ["sample_id", "umi"]].head(5).to_dict("records")
        raise ValueError(
            f"UMI counts must be numeric and non-negative. Examples: {examples}"
        )
    out["umi"] = numeric_umi.astype(float)

    invalid_sequence = ~out["cdr3"].str.match(_CANONICAL_AA)
    invalid_v = out["v_gene"].isin(["", "nan", "None"])
    if invalid_sequence.any() or invalid_v.any():
        examples = out.loc[
            invalid_sequence | invalid_v,
            ["sample_id", "cdr3", "v_gene"],
        ].head(5)
        raise ValueError(
            "The canonical table contains empty or non-functional aaV clonotypes. "
            f"Examples: {examples.to_dict('records')}"
        )

    out = out[out["umi"].gt(0)].copy()
    if out.empty:
        raise ValueError("No positive-UMI clonotypes remain after input validation.")

    out["group"] = out["group"].astype(str).str.lower().str.strip()
    out["sample_id"] = out["sample_id"].astype(str).str.strip()
    out["mouse_id"] = out["mouse_id"].astype(str).str.strip()
    out["source"] = out["source"].astype(str).str.strip()
    out["subtype"] = out["subtype"].astype(str).str.strip()
    for identifier in ("group", "sample_id", "mouse_id"):
        invalid_identifier = out[identifier].str.lower().isin({"", "nan", "none"})
        if invalid_identifier.any():
            raise ValueError(f"Column '{identifier}' contains missing identifiers.")
    out["_cell_subset"] = _infer_cell_subset(out)
    out["_tissue"] = _infer_tissue(out)
    _validate_sample_assignments(out)

    if "chain" in out.columns:
        out["chain"] = out["chain"].astype(str).str.upper().str.strip()
        out["chain"] = out["chain"].replace({"ALPHA": "TRA", "TCRA": "TRA"})
        unsupported_chains = sorted(set(out["chain"].unique()).difference({"TRA"}))
        if unsupported_chains:
            raise ValueError(
                "The active embedding workflow is TRA-specific. "
                f"Unsupported chain values: {unsupported_chains}"
            )
        out["ckey"] = out["chain"] + "|" + out["cdr3"] + "|" + out["v_gene"]
    else:
        out["ckey"] = out["cdr3"] + "|" + out["v_gene"]
    return out


def load_repertoire(path: str | Path | None = None) -> pd.DataFrame:
    """Load the later derived sample-resolved aaV Parquet interface."""
    repertoire_path = Path(path) if path else clean_repertoire_path()
    if not repertoire_path.exists():
        raise FileNotFoundError(
            f"Canonical repertoire table not found: {repertoire_path}. "
            "Set MICE_TCR_DATA_DIR or MICE_TCR_CLEAN_PARQUET."
        )

    return normalize_repertoire(
        pd.read_parquet(repertoire_path),
        source_label=str(repertoire_path),
    )


def select_stratum(df: pd.DataFrame, stratum: str) -> pd.DataFrame:
    """Select one stratum and define its independent analysis unit."""
    if stratum not in STRATA:
        raise ValueError(f"Unknown stratum '{stratum}'. Expected one of {STRATA}.")

    definition = STRATUM_DEFINITIONS[stratum]
    mask = df["_cell_subset"].isin(definition["cell_subsets"])
    mask &= df["_tissue"].isin(definition["tissues"])

    out = df.loc[mask].copy()
    if out.empty:
        raise ValueError(
            f"No clonotypes were assigned to {STRATUM_LABELS[stratum]}. "
            "Check subtype, source, and sample naming."
        )

    if definition["pool_within_mouse"]:
        out["analysis_unit"] = out["group"] + "|" + out["mouse_id"] + f"|{stratum}"
    else:
        out["analysis_unit"] = out["sample_id"]
    return out


def validate_strata(
    df: pd.DataFrame,
    strata: Sequence[str] = STRATA,
) -> pd.DataFrame:
    rows = []
    for stratum in strata:
        sub = select_stratum(df, stratum)
        rows.append(
            {
                "stratum": stratum,
                "label": STRATUM_LABELS[stratum],
                "n_rows": len(sub),
                "n_samples": sub["sample_id"].nunique(),
                "n_analysis_units": sub["analysis_unit"].nunique(),
                "n_mice": sub["mouse_id"].nunique(),
                "groups": ",".join(sorted(sub["group"].dropna().unique())),
            }
        )
    return pd.DataFrame(rows)


def requested_strata() -> tuple[str, ...]:
    raw = os.environ.get("MICE_TCR_STRATA", "").strip()
    if not raw or raw.lower() == "all":
        return STRATA

    chosen = tuple(item.strip().lower() for item in raw.split(",") if item.strip())
    unsupported = [item for item in chosen if item not in STRATA]
    if unsupported:
        raise ValueError(
            f"Unsupported strata: {unsupported}. Expected any of {STRATA}."
        )
    return chosen


def analysis_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Return one metadata row per independent analysis unit."""
    fields = ["analysis_unit", "mouse_id", "group"]
    metadata = df[fields].drop_duplicates()
    conflicts = metadata.groupby("analysis_unit").agg(
        n_mice=("mouse_id", "nunique"),
        n_groups=("group", "nunique"),
    )
    bad = conflicts[(conflicts > 1).any(axis=1)]
    if not bad.empty:
        raise ValueError(
            f"Analysis-unit metadata are inconsistent: {bad.index[:8].tolist()}"
        )
    return (
        metadata.drop_duplicates("analysis_unit")
        .set_index("analysis_unit")
        .sort_index()
    )


def ensure_groups(
    df: pd.DataFrame,
    required: Iterable[str],
    context: str,
) -> None:
    present = set(df["group"].dropna().astype(str))
    missing = [group for group in required if group not in present]
    if missing:
        raise ValueError(
            f"{context}: missing biological groups {missing}; "
            f"present={sorted(present)}."
        )
