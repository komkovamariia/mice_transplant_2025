"""Utilities for concise, result-specific stratum conclusions."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path


def format_gene_list(genes: Iterable[str], limit: int = 5) -> str:
    selected = [str(gene) for gene in genes if str(gene)]
    if not selected:
        return "no V segments passed the ranking criteria"
    return ", ".join(selected[:limit])


def save_conclusion(
    output_dir: Path,
    stratum: str,
    title: str,
    payload: dict,
    paragraphs: list[str],
) -> dict:
    """Write the machine-readable result and its publication-oriented interpretation."""
    interpretation = " ".join(part.strip() for part in paragraphs if part.strip())
    result = dict(payload)
    result["interpretation"] = interpretation

    json_path = output_dir / f"{stratum}_conclusion.json"
    json_path.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    markdown = f"# {title}\n\n{interpretation}\n"
    (output_dir / f"{stratum}_conclusion.md").write_text(
        markdown,
        encoding="utf-8",
    )
    return result
