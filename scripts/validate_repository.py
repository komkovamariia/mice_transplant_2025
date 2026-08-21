#!/usr/bin/env python3
"""Validate the restored notebook-first repository contract."""

from __future__ import annotations

import ast
import json
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CYRILLIC = re.compile(r"[\u0400-\u04FF]")
EXPECTED_STRATA = (
    "cd4_thymus",
    "cd8_thymus",
    "cd4_spleen",
    "cd8_spleen",
    "cd4_combined",
    "cd8_combined",
    "thymus_combined",
    "spleen_combined",
)
EXPECTED_AUDIT_STEMS = (
    "01_cd4_thymus",
    "02_cd8_thymus",
    "03_cd4_spleen",
    "04_cd8_spleen",
    "05_cd4_thymus_spleen",
    "06_cd8_thymus_spleen",
    "07_thymus_cd4_cd8",
    "08_spleen_cd4_cd8",
)
PRIMARY_NOTEBOOK = ROOT / "venn_original.ipynb"
LATER_NOTEBOOKS = [
    ROOT / "notebooks/02_sequence_embedding_analysis.ipynb",
    ROOT / "notebooks/03_clone_alloreactivity_analysis.ipynb",
    ROOT / "notebooks/04_cross_approach_comparison.ipynb",
]
ANALYSIS_MODULES = [
    ROOT / "src/approach2_sequence_embedding.py",
    ROOT / "src/approach3_clone_alloreactivity.py",
    ROOT / "src/approach4_cross_approach.py",
]
PYTHON_FILES = [
    *ANALYSIS_MODULES,
    ROOT / "src/strata.py",
    ROOT / "src/runtime.py",
    ROOT / "src/reporting.py",
    ROOT / "src/first_pass_input.py",
    ROOT / "src/figures.py",
    ROOT / "scripts/run_analysis.py",
    ROOT / "scripts/validate_repository.py",
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
    try:
        completed = subprocess.run(
            ["git", "-C", str(ROOT), "ls-files", "-z"],
            check=True,
            capture_output=True,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        candidates = ROOT.rglob("*")
    else:
        candidates = (
            ROOT / relative
            for relative in completed.stdout.decode("utf-8").split("\0")
            if relative
        )

    excluded_roots = {
        ".git",
        "audit_runs",
        "data",
        "logs",
        "outputs",
        "project_sources",
        "results",
    }
    excluded_parts = {
        ".ipynb_checkpoints",
        ".pytest_cache",
        ".ruff_cache",
        ".venv",
        "__pycache__",
    }
    for path in candidates:
        if not path.is_file():
            continue
        relative = path.relative_to(ROOT)
        if relative.parts and relative.parts[0] in excluded_roots:
            continue
        if any(part in excluded_parts for part in relative.parts):
            continue
        yield path, relative


def _validate_required_files() -> None:
    required = [
        ROOT / "README.md",
        ROOT / "RUN_GUIDE.md",
        PRIMARY_NOTEBOOK,
        ROOT / "scripts/run_analysis.py",
        ROOT / "scripts/validate_repository.py",
        ROOT / "src/strata.py",
        ROOT / "src/runtime.py",
        ROOT / "src/reporting.py",
        ROOT / "src/first_pass_input.py",
        ROOT / ".github/workflows/validation.yml",
        *LATER_NOTEBOOKS,
        *ANALYSIS_MODULES,
    ]
    missing = [str(path.relative_to(ROOT)) for path in required if not path.exists()]
    if missing:
        fail(f"missing required files: {missing}")

    retired = [
        ROOT / "notebooks/01_set_count_analysis.ipynb",
        ROOT / "src/approach1_set_count.py",
        ROOT / "RUN_GUIDE_RU.md",
    ]
    remaining = [str(path.relative_to(ROOT)) for path in retired if path.exists()]
    if remaining:
        fail(f"retired refactor files remain active: {remaining}")

    legacy_root = ROOT / "approaches"
    if legacy_root.exists() and any(path.is_file() for path in legacy_root.rglob("*")):
        fail("legacy approaches/ directory remains")


def _validate_language_and_paths() -> None:
    active_files = list(_active_files())
    for path, relative in active_files:
        if CYRILLIC.search(relative.as_posix()):
            fail(f"Cyrillic remains in an active path: {relative}")
        if any(character.isspace() for character in relative.name):
            fail(f"whitespace remains in an active filename: {relative}")
        if path.suffix.lower() in ACTIVE_TEXT_EXTENSIONS:
            text = path.read_text(encoding="utf-8", errors="ignore")
            if CYRILLIC.search(text):
                fail(f"Cyrillic text remains in active file: {relative}")

    for path, relative in active_files:
        if path.suffix.lower() in FIGURE_EXTENSIONS:
            if not relative.parts or relative.parts[0] != "figures":
                fail(f"figure file is stored outside figures/: {relative}")
        if "figures" in relative.parts[1:]:
            fail(f"multiple figure roots remain: {relative}")


def _compile_notebook_cells(path: Path, notebook: dict) -> None:
    for index, cell in enumerate(notebook.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        try:
            compile(source, f"{path.name}:cell-{index}", "exec")
        except SyntaxError as error:
            fail(
                f"Python syntax error in {path.relative_to(ROOT)}, "
                f"cell {index}: {error}"
            )


def _validate_primary_notebook() -> None:
    notebook = json.loads(PRIMARY_NOTEBOOK.read_text(encoding="utf-8"))
    cells = notebook.get("cells", [])
    code_cells = [cell for cell in cells if cell.get("cell_type") == "code"]
    if len(cells) < 40 or len(code_cells) < 30:
        fail(
            "venn_original.ipynb is missing substantial historical analysis content "
            f"(cells={len(cells)}, code_cells={len(code_cells)})"
        )

    identifiers = [cell.get("id") for cell in cells if cell.get("id")]
    if len(identifiers) != len(cells) or len(identifiers) != len(set(identifiers)):
        fail("venn_original.ipynb requires one unique stable ID per cell")
    if any(cell.get("outputs") for cell in code_cells):
        fail("source venn_original.ipynb contains committed outputs")
    if any(cell.get("execution_count") is not None for cell in code_cells):
        fail("source venn_original.ipynb contains execution counts")

    joined = "\n".join("".join(cell.get("source", [])) for cell in cells)
    required_tokens = (
        "6678766a0e8913ec2feacb3efb2daebd6a27eda4",
        "plot_double_triple_venn_like",
        "compare_full_vs_remaining_after_subtractions",
        "make_clone_intersection_table_venn3",
        "plot_v_landscape",
        "compare_tra_trb_from_count_tables",
        "compute_pairwise_clone_umi_correlations",
        "fisher_exact",
        "glmQLFit",
        "MICE_TCR_STRATUM",
        "_persist_open_figures",
        "figure_manifest.csv",
        "eight_stratum_distinctive_trav_summary.csv",
        *EXPECTED_STRATA,
    )
    missing = [token for token in required_tokens if token not in joined]
    if missing:
        fail(f"venn_original.ipynb lacks restored analysis token(s): {missing}")
    if re.search(r"deseq", joined, re.IGNORECASE):
        fail(
            "an excluded differential-expression implementation remains "
            "in venn_original.ipynb"
        )
    if joined.count("plt.show(") < 10:
        fail("venn_original.ipynb is missing the retained article figure set")
    required_figure_scope = (
        'region_order = ["100", "110", "101", "111"]',
        "np.log2(len(target_region_order))",
        "retention_score",
        "TRA_g1_top10_bubble_genes_heatmap.png",
        "TRB_g1_top10_bubble_genes_heatmap.png",
    )
    missing_scope = [token for token in required_figure_scope if token not in joined]
    if missing_scope:
        fail(f"venn_original.ipynb lacks article figure contract: {missing_scope}")
    forbidden_figure_scope = (
        "impact_score",
        "specific_score",
        "double_triple_g2",
        "_heatmaps.png",
        "_scatter.png",
        "_boxplot.png",
        "mouse_venn_panel",
        "distinctive_trav_across_eight_strata",
    )
    remaining_scope = [token for token in forbidden_figure_scope if token in joined]
    if remaining_scope:
        fail(f"venn_original.ipynb retains excluded figure scope: {remaining_scope}")
    if (
        "plt.savefig = _article_savefig" not in joined
        or "plt.show = _article_show" not in joined
    ):
        fail("venn_original.ipynb does not persist explicit and displayed figures")
    _compile_notebook_cells(PRIMARY_NOTEBOOK, notebook)


def _validate_later_notebooks() -> None:
    for path in LATER_NOTEBOOKS:
        notebook = json.loads(path.read_text(encoding="utf-8"))
        cells = notebook.get("cells", [])
        code_cells = [cell for cell in cells if cell.get("cell_type") == "code"]
        if len(code_cells) < 7:
            fail(
                f"{path.relative_to(ROOT)} has too few code cells "
                "for progress logging"
            )
        identifiers = [cell.get("id") for cell in cells if cell.get("id")]
        if len(identifiers) != len(set(identifiers)):
            fail(f"duplicate cell IDs in {path.relative_to(ROOT)}")
        if any(cell.get("outputs") for cell in code_cells):
            fail(
                f"source notebook contains committed outputs: "
                f"{path.relative_to(ROOT)}"
            )
        if any(cell.get("execution_count") is not None for cell in code_cells):
            fail(
                f"source notebook contains execution counts: "
                f"{path.relative_to(ROOT)}"
            )
        joined = "\n".join("".join(cell.get("source", [])) for cell in cells)
        absent = [stratum for stratum in EXPECTED_STRATA if stratum not in joined]
        if absent:
            fail(
                f"{path.relative_to(ROOT)} does not expose all eight strata: {absent}"
            )
        _compile_notebook_cells(path, notebook)


def _validate_python() -> None:
    for path in PYTHON_FILES:
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
                "every later-approach figure must be persisted in "
                f"{path.relative_to(ROOT)} "
                f"(subplots={subplot_count}, saves={save_count})"
            )


def _validate_runner_and_documentation() -> None:
    runner = (ROOT / "scripts/run_analysis.py").read_text(encoding="utf-8")
    for token in (
        "audit_runs",
        "logs",
        "[CELL ",
        "code cell(s) remain",
        "verify_first_approach_outputs",
        "prepare_first_approach_output_directories",
        "archive_first_approach_figures",
        "01_set_count.zip",
        "figure_manifest.csv",
        *EXPECTED_STRATA,
        *EXPECTED_AUDIT_STEMS,
    ):
        if token not in runner:
            fail(f"runner lacks required token: {token}")

    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    guide = (ROOT / "RUN_GUIDE.md").read_text(encoding="utf-8")
    for token in (
        "venn_original.ipynb",
        "--approach 1",
        "--approach 2",
        "--approach 3",
        "--approach 4",
        "--approach all",
        "--strata",
        "--resume",
        "audit_runs/01_cd4_thymus.executed.ipynb",
        "logs/01_cd4_thymus.log",
        "Parquet",
    ):
        if token not in readme or token not in guide:
            fail(
                f"README.md or RUN_GUIDE.md is missing command/output token: {token}"
            )

    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
    for token in ("audit_runs/", "logs/", ".ipynb_checkpoints/"):
        if token not in gitignore:
            fail(f".gitignore lacks generated path: {token}")


def main() -> int:
    _validate_required_files()
    _validate_language_and_paths()
    _validate_primary_notebook()
    _validate_later_notebooks()
    _validate_python()
    _validate_runner_and_documentation()
    print("Repository validation passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
