"""Shared data loading and biological stratification for the article analyses."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import pandas as pd

STRATA = (
    "cd4_thymus",
    "cd4_spleen",
    "cd8_thymus",
    "cd8_spleen",
    "cd4_combined",
    "cd8_combined",
)

STRATUM_LABELS = {
    "cd4_thymus": "CD4 T cells, thymus",
    "cd4_spleen": "CD4 T cells, spleen",
    "cd8_thymus": "CD8 T cells, thymus",
    "cd8_spleen": "CD8 T cells, spleen",
    "cd4_combined": "CD4 T cells, thymus + spleen",
    "cd8_combined": "CD8 T cells, thymus + spleen",
}

REQUIRED_COLUMNS = {
    "cdr3", "v_gene", "umi", "group", "sample_id", "mouse_id", "source", "subtype"
}


def repository_root() -> Path:
    env = os.environ.get("MICE_TCR_REPO")
    if env:
        return Path(env).expanduser().resolve()
    here = Path.cwd().resolve()
    for candidate in (here, *here.parents):
        if (candidate / "environment.yml").exists() and (candidate / "approaches").exists():
            return candidate
    return here


def data_dir() -> Path:
    return Path(os.environ.get("MICE_TCR_DATA_DIR", repository_root() / "data")).expanduser().resolve()


def clean_repertoire_path() -> Path:
    explicit = os.environ.get("MICE_TCR_CLEAN_PARQUET")
    return Path(explicit).expanduser().resolve() if explicit else data_dir() / "clean_clonotypes_aaV.parquet"


def load_repertoire(path: str | Path | None = None) -> pd.DataFrame:
    path = Path(path) if path else clean_repertoire_path()
    if not path.exists():
        raise FileNotFoundError(
            f"Canonical repertoire table not found: {path}. "
            "Set MICE_TCR_DATA_DIR or MICE_TCR_CLEAN_PARQUET."
        )
    df = pd.read_parquet(path)
    missing = sorted(REQUIRED_COLUMNS.difference(df.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")

    out = df.copy()
    out["cdr3"] = out["cdr3"].astype(str).str.strip()
    out["v_gene"] = out["v_gene"].astype(str).str.strip()
    out["umi"] = pd.to_numeric(out["umi"], errors="coerce").fillna(0).astype(float)
    out["group"] = out["group"].astype(str).str.lower().str.strip()
    out["sample_id"] = out["sample_id"].astype(str)
    out["mouse_id"] = out["mouse_id"].astype(str)
    out["source"] = out["source"].astype(str)
    out["subtype"] = out["subtype"].astype(str)

    out["_cell_subset"] = _infer_cell_subset(out)
    out["_tissue"] = _infer_tissue(out)
    out["ckey"] = out["cdr3"] + "|" + out["v_gene"]
    return out


def _infer_cell_subset(df: pd.DataFrame) -> pd.Series:
    text = (
        df["subtype"].fillna("").astype(str)
        + " "
        + df["sample_id"].fillna("").astype(str)
    ).str.lower()
    out = pd.Series(pd.NA, index=df.index, dtype="object")
    out[text.str.contains(r"(^|[^a-z0-9])cd4([^a-z0-9]|$)", regex=True)] = "cd4"
    out[text.str.contains(r"(^|[^a-z0-9])cd8([^a-z0-9]|$)", regex=True)] = "cd8"
    return out


def _infer_tissue(df: pd.DataFrame) -> pd.Series:
    text = (
        df["source"].fillna("").astype(str)
        + " "
        + df["sample_id"].fillna("").astype(str)
    ).str.lower()
    out = pd.Series(pd.NA, index=df.index, dtype="object")
    out[text.str.contains("spleen", regex=False)] = "spleen"
    out[text.str.contains("thym", regex=False)] = "thymus"
    return out


def select_stratum(df: pd.DataFrame, stratum: str) -> pd.DataFrame:
    if stratum not in STRATA:
        raise ValueError(f"Unknown stratum '{stratum}'. Expected one of {STRATA}.")
    cell, compartment = stratum.split("_", 1)
    mask = df["_cell_subset"].eq(cell)
    if compartment == "combined":
        mask &= df["_tissue"].isin(["thymus", "spleen"])
    else:
        mask &= df["_tissue"].eq(compartment)
    out = df.loc[mask].copy()
    if out.empty:
        raise ValueError(
            f"No clonotypes were assigned to {STRATUM_LABELS[stratum]}. "
            "Check 'subtype', 'source', and sample naming."
        )
    return out


def validate_strata(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for stratum in STRATA:
        sub = select_stratum(df, stratum)
        rows.append({
            "stratum": stratum,
            "label": STRATUM_LABELS[stratum],
            "n_rows": len(sub),
            "n_samples": sub["sample_id"].nunique(),
            "n_mice": sub["mouse_id"].nunique(),
            "groups": ",".join(sorted(sub["group"].dropna().unique())),
        })
    return pd.DataFrame(rows)


def requested_strata() -> tuple[str, ...]:
    raw = os.environ.get("MICE_TCR_STRATA", "").strip()
    if not raw or raw.lower() == "all":
        return STRATA
    chosen = tuple(x.strip().lower() for x in raw.split(",") if x.strip())
    bad = [x for x in chosen if x not in STRATA]
    if bad:
        raise ValueError(f"Unsupported strata: {bad}. Expected any of {STRATA}.")
    return chosen


def sample_metadata(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["sample_id", "mouse_id", "group", "_cell_subset", "_tissue"]
    optional = [c for c in ("source", "subtype", "treatment") if c in df.columns]
    return (
        df[cols + optional]
        .drop_duplicates(subset=["sample_id"])
        .set_index("sample_id")
        .sort_index()
    )


def ensure_groups(df: pd.DataFrame, required: Iterable[str], context: str) -> None:
    present = set(df["group"].dropna().astype(str))
    missing = [g for g in required if g not in present]
    if missing:
        raise ValueError(f"{context}: missing biological groups {missing}; present={sorted(present)}.")
