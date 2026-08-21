#!/usr/bin/env python3
"""Execute the restored notebook workflow with durable cell-level logging."""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import shutil
import threading
import time
import zipfile
from datetime import datetime
from pathlib import Path

import nbformat
import pandas as pd
from nbclient import NotebookClient


STRATUM_RUNS = (
    ("01_cd4_thymus", "cd4_thymus"),
    ("02_cd8_thymus", "cd8_thymus"),
    ("03_cd4_spleen", "cd4_spleen"),
    ("04_cd8_spleen", "cd8_spleen"),
    ("05_cd4_thymus_spleen", "cd4_combined"),
    ("06_cd8_thymus_spleen", "cd8_combined"),
    ("07_thymus_cd4_cd8", "thymus_combined"),
    ("08_spleen_cd4_cd8", "spleen_combined"),
)
STRATUM_BY_KEY = {key: stem for stem, key in STRATUM_RUNS}

APPROACHES = {
    "1": ("01_set_count", "venn_original.ipynb"),
    "2": ("02_sequence_embedding", "notebooks/02_sequence_embedding_analysis.ipynb"),
    "3": ("03_clone_alloreactivity", "notebooks/03_clone_alloreactivity_analysis.ipynb"),
    "4": ("04_cross_approach", "notebooks/04_cross_approach_comparison.ipynb"),
}
ALIASES = {
    "set-count": "1",
    "venn": "1",
    "embedding": "2",
    "clone": "3",
    "compare": "4",
    "comparison": "4",
    "cross": "4",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


class Logger:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self._lock = threading.Lock()
        self._handle = path.open("a", buffering=1, encoding="utf-8")

    def write(self, message: str) -> None:
        timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
        line = f"{timestamp} | {message}"
        with self._lock:
            print(line, flush=True)
            self._handle.write(line + "\n")
            self._handle.flush()

    def close(self) -> None:
        self._handle.close()


@contextlib.contextmanager
def temporary_env(updates: dict[str, str]):
    previous = {key: os.environ.get(key) for key in updates}
    os.environ.update(updates)
    try:
        yield
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _heartbeat(
    logger: Logger,
    stop: threading.Event,
    label: str,
    started: float,
) -> None:
    while not stop.wait(30):
        logger.write(f"{label} RUNNING | elapsed={time.monotonic() - started:.0f}s")


def _cell_description(cell) -> str:
    configured = str(cell.metadata.get("analysis_step", "")).strip()
    if configured:
        return configured[:100]
    source = "".join(cell.get("source", []))
    first_line = next(
        (line.strip() for line in source.splitlines() if line.strip()),
        "code cell",
    )
    return first_line[:100]


def _save_checkpoint(notebook, path: Path, completed_ids: set[str]) -> None:
    notebook.metadata.setdefault("article_runner", {})
    notebook.metadata["article_runner"]["completed_cell_ids"] = sorted(completed_ids)
    notebook.metadata["article_runner"]["updated_at"] = (
        datetime.now().astimezone().isoformat()
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(notebook, path)


def _restore_checkpoint_outputs(
    source_notebook,
    checkpoint_notebook,
    completed_ids: set[str],
) -> None:
    old_cells = {
        cell.id: cell for cell in checkpoint_notebook.cells if cell.cell_type == "code"
    }
    for cell in source_notebook.cells:
        if cell.cell_type != "code" or cell.id not in completed_ids:
            continue
        old = old_cells.get(cell.id)
        if old is None:
            continue
        cell["outputs"] = old.get("outputs", [])
        cell["execution_count"] = old.get("execution_count")


def execute_notebook(
    notebook_path: Path,
    artifact_stem: str,
    *,
    environment: dict[str, str],
    resume: bool,
) -> Path:
    root = repo_root()
    notebook = nbformat.read(notebook_path, as_version=4)
    output = root / "audit_runs" / f"{artifact_stem}.executed.ipynb"
    log_path = root / "logs" / f"{artifact_stem}.log"
    logger = Logger(log_path)

    code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
    current_ids = {cell.id for cell in code_cells}
    completed: set[str] = set()
    if resume and output.exists():
        checkpoint = nbformat.read(output, as_version=4)
        completed = set(
            checkpoint.metadata.get("article_runner", {}).get(
                "completed_cell_ids", []
            )
        )
        completed.intersection_update(current_ids)
        _restore_checkpoint_outputs(notebook, checkpoint, completed)
        logger.write(
            "RESUME requested | notebook state will be replayed from the first cell | "
            f"{len(completed)} prior completion marker(s) found"
        )

    run_environment = {"MICE_TCR_REPO": str(root), **environment}
    total = len(code_cells)
    ordinal = {cell.id: index + 1 for index, cell in enumerate(code_cells)}
    logger.write(
        f"START {artifact_stem} | notebook={notebook_path.relative_to(root)} | "
        f"code_cells={total} | stratum={environment.get('MICE_TCR_STRATUM', 'multiple')}"
    )

    with temporary_env(run_environment):
        client = NotebookClient(
            notebook,
            timeout=None,
            kernel_name="python3",
            resources={"metadata": {"path": str(root)}},
            allow_errors=False,
        )
        try:
            with client.setup_kernel(cwd=str(root)):
                for notebook_index, cell in enumerate(notebook.cells):
                    if cell.cell_type != "code":
                        continue

                    position = ordinal[cell.id]
                    remaining_after = total - position
                    description = _cell_description(cell)
                    label = (
                        f"[CELL {position:03d}/{total:03d}] "
                        f"id={cell.id} | step={description}"
                    )
                    phase = "REPLAY" if resume and cell.id in completed else "START"
                    started = time.monotonic()
                    logger.write(
                        f"{label} | {phase} | "
                        f"{remaining_after} code cell(s) remain after this cell"
                    )
                    stop = threading.Event()
                    heartbeat = threading.Thread(
                        target=_heartbeat,
                        args=(logger, stop, label, started),
                        daemon=True,
                    )
                    heartbeat.start()
                    try:
                        client.execute_cell(cell, notebook_index, store_history=True)
                    except BaseException:
                        stop.set()
                        heartbeat.join(timeout=1)
                        logger.write(
                            f"{label} | FAILED after {time.monotonic() - started:.1f}s"
                        )
                        _save_checkpoint(notebook, output, completed)
                        raise
                    else:
                        stop.set()
                        heartbeat.join(timeout=1)
                        completed.add(cell.id)
                        remaining_uncheckpointed = len(
                            current_ids.difference(completed)
                        )
                        logger.write(
                            f"{label} | DONE in {time.monotonic() - started:.1f}s | "
                            f"{remaining_uncheckpointed} code cell(s) not yet checkpointed"
                        )
                        _save_checkpoint(notebook, output, completed)
        finally:
            logger.write(
                f"CHECKPOINT {output.relative_to(root)} | completed={len(completed)}/{total}"
            )
            logger.close()
    return output


def resolve_sequence(value: str) -> list[str]:
    normalized = ALIASES.get(value.lower(), value.lower())
    if normalized == "all":
        return ["1", "2", "3", "4"]
    selected = [
        ALIASES.get(item.strip().lower(), item.strip().lower())
        for item in normalized.split(",")
        if item.strip()
    ]
    unsupported = [item for item in selected if item not in APPROACHES]
    if unsupported:
        raise SystemExit(f"Unknown approach value(s): {unsupported}. Use 1,2,3,4,all.")
    return list(dict.fromkeys(selected))


def resolve_strata(value: str) -> list[str]:
    if value.strip().lower() == "all":
        return [key for _, key in STRATUM_RUNS]
    selected = [item.strip().lower() for item in value.split(",") if item.strip()]
    unsupported = [item for item in selected if item not in STRATUM_BY_KEY]
    if unsupported:
        raise SystemExit(
            f"Unknown stratum value(s): {unsupported}. Expected: "
            + ", ".join(STRATUM_BY_KEY)
        )
    return list(dict.fromkeys(selected))


def _base_environment(arguments) -> dict[str, str]:
    environment = {}
    if arguments.data_dir is not None:
        environment["MICE_TCR_DATA_DIR"] = str(arguments.data_dir.resolve())
    if arguments.metadata_csv is not None:
        environment["MICE_TCR_METADATA_CSV"] = str(
            arguments.metadata_csv.resolve()
        )
    if arguments.clonoset_index is not None:
        environment["MICE_TCR_CLONOSET_INDEX"] = str(
            arguments.clonoset_index.resolve()
        )
    if arguments.working_dir is not None:
        environment["MICE_TCR_WORKING_DIR"] = str(arguments.working_dir.resolve())
    return environment


def prepare_first_approach_output_directories(selected_strata: list[str]) -> None:
    """Remove stale generated figures and tables before Approach 1 execution."""
    root = repo_root()
    figure_root = root / "figures" / "01_set_count"
    result_root = root / "results" / "01_set_count"
    archive_path = root / "figures" / "01_set_count.zip"
    if archive_path.is_file():
        archive_path.unlink()

    if set(selected_strata) == set(STRATUM_BY_KEY):
        if figure_root.is_dir():
            shutil.rmtree(figure_root)
        if result_root.is_dir():
            shutil.rmtree(result_root)
        return

    for stratum in selected_strata:
        for output_root in (figure_root, result_root):
            stratum_directory = output_root / stratum
            if stratum_directory.is_dir():
                shutil.rmtree(stratum_directory)


def archive_first_approach_figures() -> Path:
    """Create one ZIP containing the complete figures/01_set_count directory."""
    root = repo_root()
    figure_root = root / "figures" / "01_set_count"
    if not figure_root.is_dir():
        raise RuntimeError(
            "Cannot archive Approach 1 figures because figures/01_set_count is missing."
        )

    archive_path = root / "figures" / "01_set_count.zip"
    temporary_path = archive_path.with_suffix(".zip.tmp")
    if temporary_path.exists():
        temporary_path.unlink()
    with zipfile.ZipFile(
        temporary_path,
        mode="w",
        compression=zipfile.ZIP_DEFLATED,
    ) as archive:
        for path in sorted(figure_root.rglob("*")):
            if path.is_file():
                archive.write(path, path.relative_to(figure_root.parent))
    temporary_path.replace(archive_path)
    return archive_path


def verify_first_approach_outputs(selected_strata: list[str]) -> None:
    """Fail immediately when an executed notebook, table, or displayed figure is absent."""
    root = repo_root()
    missing = []
    for stratum in selected_strata:
        stem = STRATUM_BY_KEY[stratum]
        required = [
            root / "audit_runs" / f"{stem}.executed.ipynb",
            root / "logs" / f"{stem}.log",
            root / "results" / "01_set_count" / stratum / "run_summary.json",
            root / "results" / "01_set_count" / stratum / "trav_ranking.csv",
            root / "results" / "01_set_count" / stratum / "distinctive_trav.csv",
            root / "results" / "01_set_count" / stratum / "figure_manifest.csv",
        ]
        missing.extend(path for path in required if not path.is_file())

        summary_path = required[2]
        if summary_path.is_file():
            run_summary = json.loads(summary_path.read_text(encoding="utf-8"))
            edger_status = str(run_summary.get("edger_status", "missing"))
            if edger_status != "edgeR quasi-likelihood model fitted successfully.":
                raise RuntimeError(
                    f"edgeR did not complete for {stratum}: {edger_status}"
                )

        manifest_path = required[-1]
        if not manifest_path.is_file():
            continue
        manifest = pd.read_csv(manifest_path)
        if manifest.empty:
            raise RuntimeError(f"Figure manifest is empty for {stratum}.")
        heatmap_names = {
            Path(path).name
            for path in manifest["png"].dropna().astype(str)
            if "heatmap" in Path(path).name.lower()
        }
        expected_heatmap_names = {
            "TRA_g1_top10_bubble_genes_heatmap.png",
            "TRB_g1_top10_bubble_genes_heatmap.png",
        }
        if heatmap_names != expected_heatmap_names:
            raise RuntimeError(
                f"Unexpected heatmap output for {stratum}: {sorted(heatmap_names)}"
            )
        for column in ("png", "pdf"):
            for relative in manifest[column].dropna().astype(str):
                if not (root / relative).is_file():
                    missing.append(root / relative)

    if missing:
        rendered = "; ".join(str(path.relative_to(root)) for path in missing[:20])
        raise RuntimeError(f"Approach 1 output validation failed; missing: {rendered}")

    if set(selected_strata) == set(STRATUM_BY_KEY):
        complete_run_files = [
            root
            / "results"
            / "01_set_count"
            / "eight_stratum_distinctive_trav_summary.csv",
        ]
        missing_complete = [path for path in complete_run_files if not path.is_file()]
        if missing_complete:
            raise RuntimeError(
                "The complete eight-stratum summary is missing: "
                + "; ".join(str(path.relative_to(root)) for path in missing_complete)
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run the restored notebook-first analysis, selected later approaches, "
            "or the complete publication workflow."
        )
    )
    parser.add_argument(
        "--approach",
        default="all",
        help="1, 2, 3, 4/compare, comma-separated values, or all (default).",
    )
    parser.add_argument(
        "--strata",
        default="all",
        help="all or comma-separated canonical biological strata.",
    )
    parser.add_argument("--data-dir", type=Path)
    parser.add_argument("--metadata-csv", type=Path)
    parser.add_argument("--clonoset-index", type=Path)
    parser.add_argument("--working-dir", type=Path)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Replay the checkpointed notebook state and continue with durable logs.",
    )
    arguments = parser.parse_args()

    root = repo_root()
    selected_approaches = resolve_sequence(arguments.approach)
    selected_strata = resolve_strata(arguments.strata)
    base_environment = _base_environment(arguments)
    pipeline_logger = Logger(root / "logs" / "pipeline.log")
    pipeline_logger.write(
        "PIPELINE START | approaches="
        + ",".join(selected_approaches)
        + " | strata="
        + ",".join(selected_strata)
    )
    try:
        for approach_position, key in enumerate(selected_approaches, start=1):
            name, relative_path = APPROACHES[key]
            notebook_path = root / relative_path
            if key == "1":
                prepare_first_approach_output_directories(selected_strata)
                for stratum_position, stratum in enumerate(selected_strata, start=1):
                    stem = STRATUM_BY_KEY[stratum]
                    final_complete_run = (
                        set(selected_strata) == set(STRATUM_BY_KEY)
                        and stratum_position == len(selected_strata)
                    )
                    pipeline_logger.write(
                        f"STRATUM {stratum_position}/{len(selected_strata)} START | "
                        f"{stratum} | output=audit_runs/{stem}.executed.ipynb"
                    )
                    environment = {
                        **base_environment,
                        "MICE_TCR_STRATUM": stratum,
                        "MICE_TCR_FINALIZE_SUMMARY": (
                            "1" if final_complete_run else "0"
                        ),
                    }
                    execute_notebook(
                        notebook_path,
                        stem,
                        environment=environment,
                        resume=arguments.resume,
                    )
                    pipeline_logger.write(
                        f"STRATUM {stratum_position}/{len(selected_strata)} DONE | {stratum}"
                    )
                verify_first_approach_outputs(selected_strata)
                archive_path = archive_first_approach_figures()
                pipeline_logger.write(
                    f"FIGURE ARCHIVE CREATED | {archive_path.relative_to(root)}"
                )
                pipeline_logger.write(
                    f"APPROACH {approach_position}/{len(selected_approaches)} DONE | "
                    f"{name} | verified_strata={len(selected_strata)}"
                )
                continue

            remaining = len(selected_approaches) - approach_position
            pipeline_logger.write(
                f"APPROACH {approach_position}/{len(selected_approaches)} START | "
                f"{name} | {remaining} approach(es) remain after this"
            )
            environment = {
                **base_environment,
                "MICE_TCR_STRATA": ",".join(selected_strata),
            }
            output = execute_notebook(
                notebook_path,
                name,
                environment=environment,
                resume=arguments.resume,
            )
            pipeline_logger.write(
                f"APPROACH {approach_position}/{len(selected_approaches)} DONE | "
                f"{name} | output={output.relative_to(root)}"
            )
    finally:
        pipeline_logger.write("PIPELINE END")
        pipeline_logger.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
