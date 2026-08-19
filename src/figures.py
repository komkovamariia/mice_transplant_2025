"""Consistent publication-oriented figure output."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt

from .strata import repository_root


def _apply_style() -> None:
    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "figure.dpi": 120,
            "font.size": 10,
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def figure_dir(
    approach: str,
    stratum: str | None = None,
) -> Path:
    root = repository_root() / "figures" / approach
    path = root / stratum if stratum else root
    path.mkdir(parents=True, exist_ok=True)
    return path


def slugify(text: str) -> str:
    normalized = (
        re.sub(
            r"[^A-Za-z0-9]+",
            "_",
            str(text),
        )
        .strip("_")
        .lower()
    )
    return normalized or "figure"


def _title(fig) -> str:
    figure_title = fig._suptitle.get_text() if fig._suptitle else ""
    axis_titles = [
        axis.get_title().strip() for axis in fig.axes if axis.get_title().strip()
    ]
    return figure_title.strip() or "; ".join(axis_titles) or "TCR analysis figure"


def save_figure(
    fig,
    approach: str,
    stratum: str | None,
    name: str,
    *,
    close: bool = True,
    dpi: int = 300,
):
    """Save every figure as a raster and editable vector file."""
    _apply_style()
    output = figure_dir(approach, stratum)
    stem = slugify(name)
    png_path = output / f"{stem}.png"
    pdf_path = output / f"{stem}.pdf"
    title = _title(fig)

    fig.savefig(
        png_path,
        dpi=dpi,
        bbox_inches="tight",
        metadata={"Title": title, "Software": "mice_transplant_2025"},
    )
    fig.savefig(
        pdf_path,
        bbox_inches="tight",
        metadata={"Title": title, "Creator": "mice_transplant_2025"},
    )
    if close:
        plt.close(fig)
    return png_path, pdf_path


def write_figure_inventory() -> Path:
    """Write one auditable inventory for every persisted publication figure."""
    root = repository_root()
    figure_root = root / "figures"
    output = root / "results" / "figure_inventory.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for png_path in sorted(figure_root.rglob("*.png")):
        relative = png_path.relative_to(figure_root)
        parts = relative.parts
        approach = parts[0] if len(parts) >= 2 else ""
        stratum = parts[1] if len(parts) >= 3 else "cross_stratum"
        pdf_path = png_path.with_suffix(".pdf")
        rows.append(
            {
                "approach": approach,
                "stratum": stratum,
                "figure": png_path.stem,
                "png": png_path.relative_to(root).as_posix(),
                "pdf": (
                    pdf_path.relative_to(root).as_posix() if pdf_path.exists() else ""
                ),
            }
        )
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["approach", "stratum", "figure", "png", "pdf"],
        )
        writer.writeheader()
        writer.writerows(rows)
    return output
