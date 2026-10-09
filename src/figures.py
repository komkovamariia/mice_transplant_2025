"""Publication figure helpers used across the repertoire analyses."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt

from .strata import repository_root


CORAL = "#FF785C"
PURPLE = "#792080"
MAGENTA = "#BF1674"
LINK_PINK = "#D47BAE"
THRESHOLD_BLUE = "#5B9BD5"


def apply_article_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 14,
            "axes.titleweight": "normal",
            "axes.labelsize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 120,
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white",
        }
    )


def figure_dir(
    approach: str,
    stratum: str | None = None,
    chain: str | None = None,
) -> Path:
    root = repository_root() / "figures" / approach
    path = root / stratum if stratum else root
    if chain:
        normalized = str(chain).upper()
        if normalized not in {"TRA", "TRB"}:
            raise ValueError("chain must be TRA or TRB")
        path = path / normalized
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
    chain: str | None = None,
    close: bool = True,
    dpi: int = 300,
) -> Path:
    """Save one publication PNG in the analysis figure directory."""
    apply_article_style()
    output = figure_dir(approach, stratum, chain)
    png_path = output / f"{slugify(name)}.png"
    fig.savefig(
        png_path,
        dpi=dpi,
        bbox_inches="tight",
        metadata={"Title": _title(fig), "Software": "mice_transplant_2025"},
    )
    if close:
        plt.close(fig)
    return png_path


def write_figure_inventory() -> Path:
    """Write one inventory row for every persisted PNG figure."""
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
        chain = parts[2] if len(parts) >= 4 and parts[2] in {"TRA", "TRB"} else ""
        rows.append(
            {
                "approach": approach,
                "stratum": stratum,
                "chain": chain,
                "figure": png_path.stem,
                "png": png_path.relative_to(root).as_posix(),
            }
        )
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["approach", "stratum", "chain", "figure", "png"],
        )
        writer.writeheader()
        writer.writerows(rows)
    return output
