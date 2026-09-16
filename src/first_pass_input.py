"""Verified CSV and MiXCR input route for the count-based first approach."""

from __future__ import annotations

import os
import re
from pathlib import Path

import pandas as pd

from .strata import normalize_repertoire

DEFAULT_METADATA_CSV = Path(
    "/projects/mice_transplant_2025/metadata_mice_transplant.csv"
)
DEFAULT_CLONOSET_INDEX = Path(
    "/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv"
)
HISTORICAL_EXCLUSIONS = {
    "g1_m3_thymus_cd8_80_alpha",
    "g5_m1_spleen_cd4_93_alpha",
    "g5_m4_thymus_cd8_100_alpha",
}
_CANONICAL_AA = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def _identifier(value) -> str:
    if pd.isna(value):
        return ""
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value).strip()
    if numeric.is_integer():
        return str(int(numeric))
    return str(value).strip()


def first_pass_paths(
    metadata_csv: str | Path | None = None,
    clonoset_index: str | Path | None = None,
) -> tuple[Path, Path]:
    metadata = Path(
        metadata_csv or os.environ.get("MICE_TCR_METADATA_CSV", DEFAULT_METADATA_CSV)
    ).expanduser()
    index = Path(
        clonoset_index
        or os.environ.get("MICE_TCR_CLONOSET_INDEX", DEFAULT_CLONOSET_INDEX)
    ).expanduser()
    return metadata.resolve(), index.resolve()


def _resolve_clonoset_file(value, index_path: Path) -> Path:
    candidate = Path(str(value)).expanduser()
    candidates = [candidate]
    if not candidate.is_absolute():
        candidates.extend(
            [
                index_path.parent / candidate,
                index_path.parent / "mixcr" / candidate,
            ]
        )
    for path in candidates:
        if path.exists():
            return path.resolve()
    return candidates[-1].resolve()


def prepare_sample_index(metadata_path: Path, index_path: Path) -> pd.DataFrame:
    """Reproduce the sample-ID merge used by the original first-pass notebook."""
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Metadata CSV not found: {metadata_path}")
    if not index_path.is_file():
        raise FileNotFoundError(f"Clonoset-index CSV not found: {index_path}")

    metadata = pd.read_csv(metadata_path)
    required_metadata = {
        "chain",
        "sample_no",
        "sample_id",
        "group_no",
        "mouse_no",
        "source",
        "subtype",
    }
    missing_metadata = sorted(required_metadata.difference(metadata.columns))
    if missing_metadata:
        raise ValueError(
            f"{metadata_path} is missing metadata columns: {missing_metadata}"
        )
    metadata = metadata.copy()
    # the historical index keys use the metadata prefixes verbatim, for example
    # alpha-100 and beta-100. Chain standardization belongs after this merge
    metadata["sample_id_old"] = (
        metadata["chain"].astype(str).str.strip()
        + "-"
        + metadata["sample_no"].map(_identifier)
    )
    metadata = metadata.drop(columns="chain")

    clonosets = pd.read_csv(index_path)
    required_index = {"sample_id", "filename", "chain"}
    missing_index = sorted(required_index.difference(clonosets.columns))
    if missing_index:
        raise ValueError(
            f"{index_path} is missing clonoset-index columns: {missing_index}"
        )
    clonosets = clonosets.rename(columns={"sample_id": "sample_id_old"})
    merged = clonosets.merge(
        metadata,
        on="sample_id_old",
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    unmatched = merged.loc[merged["_merge"].ne("both"), "sample_id_old"].tolist()
    if unmatched:
        raise ValueError(
            "Clonoset rows could not be matched to metadata by chain and sample number: "
            + ", ".join(map(str, unmatched[:10]))
        )
    merged = merged.drop(columns=["_merge", "sample_id_old"])
    merged["chain"] = merged["chain"].astype(str).str.upper().str.strip()
    merged["group_no"] = merged["group_no"].astype(str).str.lower().str.strip()
    merged = merged[
        merged["chain"].eq("TRA")
        & merged["group_no"].isin(["g1", "g2", "g5", "g6"])
        & ~merged["sample_id"].isin(HISTORICAL_EXCLUSIONS)
    ].copy()
    if merged.empty:
        raise ValueError(
            "No TRA samples remained after the first-pass metadata selection."
        )

    merged["filename"] = merged["filename"].map(
        lambda value: str(_resolve_clonoset_file(value, index_path))
    )
    missing_files = [path for path in merged["filename"] if not Path(path).is_file()]
    if missing_files:
        raise FileNotFoundError(
            "MiXCR clonotype export(s) referenced by the index were not found: "
            + "; ".join(missing_files[:10])
        )
    if merged["sample_id"].duplicated().any():
        duplicated = merged.loc[merged["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"Duplicate TRA sample identifiers: {duplicated[:10]}")
    return merged.sort_values("sample_id").reset_index(drop=True)


def count_table_to_repertoire(
    count_table: pd.DataFrame,
    sample_index: pd.DataFrame,
) -> pd.DataFrame:
    """Convert a repseq aaV count table to the shared sample-resolved interface."""
    features = [tuple(feature) for feature in count_table.index]
    invalid_features = [feature for feature in features if len(feature) != 2]
    if invalid_features:
        raise ValueError(
            "The repseq count table must use aaV feature tuples: "
            f"{invalid_features[:5]}"
        )
    matrix = count_table.copy()
    matrix.index = pd.MultiIndex.from_tuples(features, names=["cdr3", "v_gene"])
    long = matrix.reset_index().melt(
        id_vars=["cdr3", "v_gene"],
        var_name="sample_id",
        value_name="umi",
    )
    long["umi"] = pd.to_numeric(long["umi"], errors="coerce")
    long = long[long["umi"].gt(0)].copy()
    functional = long["cdr3"].fillna("").astype(str).str.match(_CANONICAL_AA)
    valid_v = long["v_gene"].fillna("").astype(str).str.startswith("TRAV")
    removed = int((~(functional & valid_v)).sum())
    long = long[functional & valid_v].copy()
    if long.empty:
        raise ValueError("No functional TRA aaV clonotypes remained after validation.")

    metadata_columns = [
        "sample_id",
        "group_no",
        "mouse_no",
        "source",
        "subtype",
        "chain",
    ]
    long = long.merge(
        sample_index[metadata_columns].drop_duplicates("sample_id"),
        on="sample_id",
        how="left",
        validate="many_to_one",
    )
    long = long.rename(columns={"group_no": "group", "mouse_no": "mouse_id"})
    long["mouse_id"] = long["mouse_id"].map(_identifier)
    long["abundance_measure"] = "repseq_count"
    if removed:
        print(
            f"Excluded {removed:,} non-functional or non-TRAV aaV rows; "
            "no minimum abundance threshold was applied."
        )
    return normalize_repertoire(
        long,
        source_label="repseq aaV count table generated from MiXCR exports",
    )


def load_first_pass_repertoire(
    metadata_csv: str | Path | None = None,
    clonoset_index: str | Path | None = None,
) -> pd.DataFrame:
    """Load Approach 1 directly from the verified CSV and MiXCR data route."""
    try:
        from repseq import intersections
    except ImportError as error:
        raise ImportError(
            "Approach 1 requires the pinned repseq package from environment.yml. "
            "Update the Conda environment before execution."
        ) from error

    metadata_path, index_path = first_pass_paths(metadata_csv, clonoset_index)
    sample_index = prepare_sample_index(metadata_path, index_path)
    count_table = intersections.count_table(
        sample_index,
        overlap_type="aaV",
        mismatches=0,
    )
    return count_table_to_repertoire(count_table, sample_index)
