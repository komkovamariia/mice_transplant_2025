#!/usr/bin/env python3
"""Static validation of the publication-oriented repository contract."""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CYRILLIC = re.compile(r"[\u0400-\u04FF]")
EXPECTED_STRATA = {
    "cd4_thymus",
    "cd8_thymus",
    "cd4_spleen",
    "cd8_spleen",
    "cd4_combined",
    "cd8_combined",
    "thymus_combined",
    "spleen_combined",
}
EXPECTED_NOTEBOOKS = [
    ROOT / "notebooks/01_set_count_analysis.ipynb",
    ROOT / "notebooks/02_sequence_embedding_analysis.ipynb",
    ROOT / "notebooks/03_clone_alloreactivity_analysis.ipynb",
    ROOT / "notebooks/04_cross_approach_comparison.ipynb",
]
ANALYSIS_MODULES = [
    ROOT / "src/approach1_set_count.py",
    ROOT / "src/approach2_sequence_embedding.py",
    ROOT / "src/approach3_clone_alloreactivity.py",
    ROOT / "src/approach4_cross_approach.py",
]
ACTIVE_TEXT_EXTENSIONS = {
    ".md",
    ".py",
    ".ipynb",
    ".yml",
    ".yaml",
    ".txt",
    ".env",
}
FIGURE_EXTENSIONS = {
    ".png",
    ".pdf",
    ".svg",
    ".jpg",
    ".jpeg",
    ".tif",
    ".tiff",
}


def fail(message: str) -> None:
    raise SystemExit(f"VALIDATION FAILED: {message}")


def _active_files():
    excluded_roots = {
        ".git",
        "data",
        "outputs",
        "results",
        "project_sources",
    }
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if relative.parts and relative.parts[0] in excluded_roots:
            continue
        yield path, relative


def _validate_required_files() -> None:
    required = [
        ROOT / "README.md",
        ROOT / "RUN_GUIDE.md",
        ROOT / "scripts/run_analysis.py",
        ROOT / "scripts/validate_repository.py",
        ROOT / "src/strata.py",
        ROOT / "src/runtime.py",
        ROOT / "src/reporting.py",
        ROOT / "src/first_pass_input.py",
        ROOT / ".github/workflows/validation.yml",
        *EXPECTED_NOTEBOOKS,
        *ANALYSIS_MODULES,
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        fail(f"missing required files: {missing}")
    if (ROOT / "RUN_GUIDE_RU.md").exists():
        fail("legacy Russian execution guide remains in the active tree")
    if (ROOT / "approaches").exists():
        fail(
            "legacy approaches/ directory remains; source notebooks belong in notebooks/"
        )


def _validate_language_and_paths() -> None:
    figure_roots = [
        path
        for path in ROOT.rglob("figures")
        if path.is_dir() and path != ROOT / "figures"
    ]
    if figure_roots:
        fail(
            "multiple figure roots remain: "
            + ", ".join(str(path.relative_to(ROOT)) for path in figure_roots)
        )
    for path, relative in _active_files():
        relative_text = relative.as_posix()
        if CYRILLIC.search(relative_text):
            fail(f"Cyrillic remains in an active path: {relative}")
        if any(character.isspace() for character in relative.name):
            fail(f"whitespace remains in an active filename: {relative}")

        if path.suffix.lower() not in ACTIVE_TEXT_EXTENSIONS:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if CYRILLIC.search(text):
            fail(f"Cyrillic text remains in active file: {relative}")

    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in FIGURE_EXTENSIONS:
            continue
        relative = path.relative_to(ROOT)
        if not relative.parts or relative.parts[0] != "figures":
            fail(f"figure file is stored outside figures/: {relative}")


def _validate_python() -> None:
    for path in [
        *ANALYSIS_MODULES,
        ROOT / "src/strata.py",
        ROOT / "src/runtime.py",
        ROOT / "src/reporting.py",
        ROOT / "src/first_pass_input.py",
        ROOT / "src/figures.py",
        ROOT / "scripts/run_analysis.py",
    ]:
        try:
            ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError as error:
            fail(f"Python syntax error in {path.relative_to(ROOT)}: {error}")

    for path in ANALYSIS_MODULES:
        text = path.read_text(encoding="utf-8")
        if "plt.show(" in text:
            fail(f"inline-only plotting remains in {path.relative_to(ROOT)}")
        subplot_count = text.count("plt.subplots(")
        save_count = text.count("save_figure(")
        if subplot_count != save_count:
            fail(
                "every matplotlib figure must be persisted in "
                f"{path.relative_to(ROOT)} "
                f"(subplots={subplot_count}, saves={save_count})"
            )

    required_effect = "effect_g1_vs_allogeneic"
    for path in ANALYSIS_MODULES[:3]:
        if required_effect not in path.read_text(encoding="utf-8"):
            fail(f"{path.relative_to(ROOT)} lacks the standardized signed-effect field")


def _validate_notebooks() -> None:
    for path in EXPECTED_NOTEBOOKS:
        notebook = json.loads(path.read_text(encoding="utf-8"))
        cells = notebook.get("cells", [])
        code_cells = [cell for cell in cells if cell.get("cell_type") == "code"]
        if len(code_cells) < 7:
            fail(
                f"{path.relative_to(ROOT)} has too few code cells "
                "for stratum-level progress"
            )

        identifiers = [cell.get("id") for cell in cells if cell.get("id")]
        if len(identifiers) != len(set(identifiers)):
            fail(f"duplicate cell IDs in {path.relative_to(ROOT)}")

        joined = "\n".join("".join(cell.get("source", [])) for cell in cells)
        absent = {stratum for stratum in EXPECTED_STRATA if stratum not in joined}
        if absent:
            fail(
                f"{path.relative_to(ROOT)} does not expose all eight "
                f"strata: {sorted(absent)}"
            )
        if "DESeq2" in joined or "deseq" in joined.lower():
            fail(f"stale DESeq2 reference in {path.relative_to(ROOT)}")
        if any(cell.get("outputs") for cell in code_cells):
            fail(
                f"source notebook contains committed outputs: {path.relative_to(ROOT)}"
            )
        if any(cell.get("execution_count") is not None for cell in code_cells):
            fail(f"source notebook contains execution counts: {path.relative_to(ROOT)}")
        plotting_tokens = (
            "plt.show(",
            "plt.subplots(",
            "sns.",
            ".plot(",
        )
        if any(token in joined for token in plotting_tokens):
            fail(
                f"plotting remains in notebook instead of src/: "
                f"{path.relative_to(ROOT)}"
            )


def _validate_documentation_and_runtime() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    guide = (ROOT / "RUN_GUIDE.md").read_text(encoding="utf-8")
    for token in (
        "--approach 1",
        "--approach 2",
        "--approach 3",
        "--approach 4",
        "--approach all",
        "--strata",
        "--resume",
        "pipeline.log",
    ):
        if token not in readme or token not in guide:
            fail(f"README.md or RUN_GUIDE.md is missing command token: {token}")

    runner = (ROOT / "scripts/run_analysis.py").read_text(encoding="utf-8")
    for token in (
        "[CELL ",
        "remaining",
        "DONE",
        "pipeline.log",
        "id=",
        "step=",
    ):
        if token not in runner:
            fail(f"runner lacks progress-log token: {token}")

    strata = (ROOT / "src/strata.py").read_text(encoding="utf-8")
    for token in (
        "analysis_unit",
        "mouse-level thymus + spleen pool",
        *EXPECTED_STRATA,
    ):
        if token not in strata:
            fail(f"stratification code lacks required token: {token}")

    approach1 = (ROOT / "src/approach1_set_count.py").read_text(encoding="utf-8")
    for token in (
        "01_exact_clonotype_overlap",
        "02_g1_exclusive_trav_usage",
        "03_edger_volcano",
        "04_edger_volcano_labeled",
        "05_edger_ma",
        "06_edger_ma_labeled",
        "07_top_g1_enriched_clonotypes",
        "08_trav_mean_effect",
        "09_trav_significant_feature_count",
        "10_fisher_edger_concordance",
        "eight_stratum_summary.csv",
        "distinctive_trav_summary.csv",
    ):
        if token not in approach1:
            fail(f"Approach 1 lacks required output token: {token}")

    first_notebook = EXPECTED_NOTEBOOKS[0].read_text(encoding="utf-8")
    if "load_first_pass_repertoire" not in first_notebook:
        fail("Approach 1 notebook does not use the verified CSV and MiXCR loader")
    if "load_repertoire()" in first_notebook:
        fail("Approach 1 notebook still invokes the derived Parquet loader")

    environment = (ROOT / "environment.yml").read_text(encoding="utf-8")
    requirements = (ROOT / "requirements.txt").read_text(encoding="utf-8")
    if "bioconductor-deseq2" in environment.lower():
        fail("DESeq2 remains in environment.yml")
    for token in ("nbclient", "pytest"):
        if token not in environment.lower():
            fail(f"environment.yml is missing {token}")
        if token not in requirements.lower():
            fail(f"requirements.txt is missing {token}")


def main() -> int:
    _validate_required_files()
    _validate_language_and_paths()
    _validate_python()
    _validate_notebooks()
    _validate_documentation_and_runtime()
    print("Repository validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
