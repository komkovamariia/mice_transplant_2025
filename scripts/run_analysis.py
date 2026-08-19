#!/usr/bin/env python3
"""Execute one or more analysis notebooks with cell-level durable progress."""

from __future__ import annotations

import argparse
import contextlib
import os
import threading
import time
from datetime import datetime
from pathlib import Path

import nbformat
from nbclient import NotebookClient

APPROACHES = {
    "1": (
        "01_set_count",
        "approaches/01_set_count/set_count_analysis.ipynb",
    ),
    "2": (
        "02_sequence_embedding",
        "approaches/02_sequence_embedding/sequence_embedding_analysis.ipynb",
    ),
    "3": (
        "03_clone_alloreactivity",
        "approaches/03_clone_alloreactivity/clone_alloreactivity_analysis.ipynb",
    ),
    "4": (
        "04_cross_approach",
        "approaches/04_cross_approach/cross_approach_comparison.ipynb",
    ),
}
ALIASES = {
    "set-count": "1",
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
        self._handle = path.open(
            "a",
            buffering=1,
            encoding="utf-8",
        )

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


def _save_checkpoint(
    notebook,
    path: Path,
    completed_ids: set[str],
) -> None:
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
    approach_name: str,
    *,
    strata: str,
    data_dir: Path | None,
    resume: bool,
) -> Path:
    root = repo_root()
    notebook = nbformat.read(notebook_path, as_version=4)
    output = root / "outputs" / "notebooks" / f"{approach_name}.executed.ipynb"
    log_path = root / "outputs" / "logs" / f"{approach_name}.log"
    logger = Logger(log_path)

    code_cells = [cell for cell in notebook.cells if cell.cell_type == "code"]
    current_ids = {cell.id for cell in code_cells}
    completed: set[str] = set()
    if resume and output.exists():
        checkpoint = nbformat.read(output, as_version=4)
        completed = set(
            checkpoint.metadata.get(
                "article_runner",
                {},
            ).get("completed_cell_ids", [])
        )
        completed.intersection_update(current_ids)
        _restore_checkpoint_outputs(
            notebook,
            checkpoint,
            completed,
        )
        logger.write(
            f"RESUME requested | {len(completed)} valid completed cell(s) restored"
        )

    environment = {
        "MICE_TCR_REPO": str(root),
        "MICE_TCR_STRATA": strata,
    }
    if data_dir is not None:
        environment["MICE_TCR_DATA_DIR"] = str(data_dir.resolve())

    total = len(code_cells)
    ordinal = {cell.id: index + 1 for index, cell in enumerate(code_cells)}
    logger.write(
        f"START {approach_name} | "
        f"notebook={notebook_path.relative_to(root)} | "
        f"code_cells={total} | strata={strata}"
    )

    with temporary_env(environment):
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
                    bootstrap = "bootstrap" in set(cell.metadata.get("tags", []))
                    description = _cell_description(cell)
                    label = (
                        f"[CELL {position:03d}/{total:03d}] "
                        f"id={cell.id} | step={description}"
                    )
                    if resume and cell.id in completed and not bootstrap:
                        logger.write(
                            f"{label} | SKIP completed checkpoint | "
                            f"{remaining_after} position(s) remain"
                        )
                        continue

                    started = time.monotonic()
                    logger.write(
                        f"{label} | START | "
                        f"{remaining_after} position(s) remain after this cell"
                    )
                    stop = threading.Event()
                    heartbeat = threading.Thread(
                        target=_heartbeat,
                        args=(logger, stop, label, started),
                        daemon=True,
                    )
                    heartbeat.start()
                    try:
                        client.execute_cell(
                            cell,
                            notebook_index,
                            store_history=True,
                        )
                    except BaseException:
                        stop.set()
                        heartbeat.join(timeout=1)
                        logger.write(
                            f"{label} | FAILED after {time.monotonic() - started:.1f}s"
                        )
                        _save_checkpoint(
                            notebook,
                            output,
                            completed,
                        )
                        raise
                    else:
                        stop.set()
                        heartbeat.join(timeout=1)
                        completed.add(cell.id)
                        remaining_uncheckpointed = len(
                            current_ids.difference(completed)
                        )
                        logger.write(
                            f"{label} | DONE in "
                            f"{time.monotonic() - started:.1f}s | "
                            f"{remaining_uncheckpointed} code cell(s) "
                            "not yet checkpointed"
                        )
                        _save_checkpoint(
                            notebook,
                            output,
                            completed,
                        )
        finally:
            logger.write(
                f"CHECKPOINT {output.relative_to(root)} | "
                f"completed={len(completed)}/{total}"
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


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run one analytical approach, selected approaches, "
            "or the complete publication pipeline."
        )
    )
    parser.add_argument(
        "--approach",
        default="all",
        help=("1, 2, 3, 4/compare, comma-separated values, or all (default)."),
    )
    parser.add_argument(
        "--strata",
        default="all",
        help="all or comma-separated canonical biological strata.",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Directory containing clean_clonotypes_aaV.parquet.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Restore checkpoint outputs, re-run bootstrap cells, "
            "and skip completed analytical cells."
        ),
    )
    arguments = parser.parse_args()

    root = repo_root()
    selected = resolve_sequence(arguments.approach)
    pipeline_logger = Logger(root / "outputs" / "logs" / "pipeline.log")
    pipeline_logger.write(
        f"PIPELINE START | approaches={','.join(selected)} | strata={arguments.strata}"
    )
    try:
        for position, key in enumerate(selected, start=1):
            name, relative_path = APPROACHES[key]
            remaining = len(selected) - position
            pipeline_logger.write(
                f"APPROACH {position}/{len(selected)} START | "
                f"{name} | {remaining} approach(es) remain after this"
            )
            try:
                output = execute_notebook(
                    root / relative_path,
                    name,
                    strata=arguments.strata,
                    data_dir=arguments.data_dir,
                    resume=arguments.resume,
                )
            except BaseException:
                pipeline_logger.write(
                    f"APPROACH {position}/{len(selected)} FAILED | {name}"
                )
                raise
            pipeline_logger.write(
                f"APPROACH {position}/{len(selected)} DONE | "
                f"{name} | output={output.relative_to(root)}"
            )
    finally:
        pipeline_logger.write("PIPELINE END")
        pipeline_logger.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
