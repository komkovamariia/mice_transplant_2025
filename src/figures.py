"""Consistent publication-oriented figure output."""

from __future__ import annotations

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
