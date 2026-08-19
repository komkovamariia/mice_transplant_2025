#!/usr/bin/env python3
"""Static repository validation for the article analysis layout."""

from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CYRILLIC = re.compile(r"[\u0400-\u04FF]")
EXPECTED_NOTEBOOKS = [
    ROOT / "approaches/01_set_count/set_count_analysis.ipynb",
    ROOT / "approaches/02_sequence_embedding/sequence_embedding_analysis.ipynb",
    ROOT / "approaches/03_clone_alloreactivity/clone_alloreactivity_analysis.ipynb",
    ROOT / "approaches/04_cross_approach/cross_approach_comparison.ipynb",
]
EXPECTED_STRATA = {"cd4_thymus", "cd4_spleen", "cd8_thymus", "cd8_spleen", "cd4_combined", "cd8_combined"}


def fail(message: str) -> None:
    raise SystemExit(f"VALIDATION FAILED: {message}")


def main() -> int:
    required = [ROOT / "README.md", ROOT / "RUN_GUIDE_RU.md", ROOT / "scripts/run_analysis.py",
                ROOT / "src/strata.py", ROOT / "src/runtime.py", *EXPECTED_NOTEBOOKS]
    missing = [str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:
        fail(f"missing required files: {missing}")

    extensions = {".md", ".py", ".ipynb", ".yml", ".yaml", ".txt"}
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in extensions:
            continue
        rel = path.relative_to(ROOT)
        if rel.as_posix() == "RUN_GUIDE_RU.md":
            continue
        if any(part in {".git", "outputs"} for part in rel.parts):
            continue
        text = path.read_text(errors="ignore")
        if CYRILLIC.search(text):
            fail(f"Cyrillic text remains in active English file: {rel}")

    for path in EXPECTED_NOTEBOOKS:
        nb = json.loads(path.read_text())
        code = [c for c in nb["cells"] if c["cell_type"] == "code"]
        if len(code) < 7:
            fail(f"{path.relative_to(ROOT)} has too few code cells for stratum-level progress.")
        joined = "\n".join("".join(c.get("source", [])) for c in nb["cells"])
        absent = EXPECTED_STRATA.difference(x for x in EXPECTED_STRATA if x in joined)
        if absent:
            fail(f"{path.relative_to(ROOT)} does not expose all six strata: {sorted(absent)}")
        if "DESeq2" in joined or "deseq" in joined.lower():
            fail(f"stale DESeq2 reference in {path.relative_to(ROOT)}")

    analysis_modules = [ROOT / "src/approach1_set_count.py", ROOT / "src/approach2_sequence_embedding.py",
                        ROOT / "src/approach3_clone_alloreactivity.py", ROOT / "src/approach4_cross_approach.py"]
    for path in analysis_modules:
        text = path.read_text()
        if "plt.show(" in text:
            fail(f"inline-only plotting remains in {path.relative_to(ROOT)}")
        if text.count("plt.subplots(") != text.count("save_figure("):
            fail(f"every generated matplotlib figure must be persisted in {path.relative_to(ROOT)} "
                 f"(subplots={text.count('plt.subplots(')}, saves={text.count('save_figure(')})")

    readme = (ROOT / "README.md").read_text()
    for token in ("--approach 1", "--approach 2", "--approach 3", "--approach all", "--resume"):
        if token not in readme:
            fail(f"README is missing launch command token: {token}")

    env = (ROOT / "environment.yml").read_text()
    if "bioconductor-deseq2" in env.lower():
        fail("DESeq2 remains in environment.yml")
    if "nbclient" not in env.lower():
        fail("nbclient is required by scripts/run_analysis.py but missing from environment.yml")

    print("Repository validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
