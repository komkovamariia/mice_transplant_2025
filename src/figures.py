"""Publication-oriented figure output utilities."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt

from .strata import repository_root


def figure_dir(approach: str, stratum: str | None = None) -> Path:
    root = repository_root() / "figures" / approach
    path = root / stratum if stratum else root
    path.mkdir(parents=True, exist_ok=True)
    return path


def slugify(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_").lower()
    return text or "figure"


def save_figure(fig, approach: str, stratum: str | None, name: str, *, close: bool = True, dpi: int = 300):
    out = figure_dir(approach, stratum)
    stem = slugify(name)
    png = out / f"{stem}.png"
    pdf = out / f"{stem}.pdf"
    fig.savefig(png, dpi=dpi, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    if close:
        plt.close(fig)
    return png, pdf
