#!/usr/bin/env python3
"""Execute article notebooks with cell-level progress, checkpoints, and optional resume."""

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
    "1": ("01_set_count", "approaches/01_set_count/set_count_analysis.ipynb"),
    "2": ("02_sequence_embedding", "approaches/02_sequence_embedding/sequence_embedding_analysis.ipynb"),
    "3": ("03_clone_alloreactivity", "approaches/03_clone_alloreactivity/clone_alloreactivity_analysis.ipynb"),
    "4": ("04_cross_approach", "approaches/04_cross_approach/cross_approach_comparison.ipynb"),
}
ALIASES = {"set-count": "1", "embedding": "2", "clone": "3", "compare": "4", "cross": "4"}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


class Logger:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self._lock = threading.Lock()
        self._fh = path.open("a", buffering=1)

    def write(self, message: str) -> None:
        stamp = datetime.now().astimezone().isoformat(timespec="seconds")
        line = f"{stamp} | {message}"
        with self._lock:
            print(line, flush=True)
            self._fh.write(line + "\n")
            self._fh.flush()

    def close(self) -> None:
        self._fh.close()


@contextlib.contextmanager
def temporary_env(updates: dict[str, str]):
    old = {k: os.environ.get(k) for k in updates}
    os.environ.update(updates)
    try:
        yield
    finally:
        for key, value in old.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


def _heartbeat(logger: Logger, stop: threading.Event, label: str, started: float) -> None:
    while not stop.wait(30):
        logger.write(f"{label} RUNNING | elapsed={time.monotonic() - started:.0f}s")


def _save_checkpoint(nb, path: Path, completed_ids: set[str]) -> None:
    nb.metadata.setdefault("article_runner", {})
    nb.metadata["article_runner"]["completed_cell_ids"] = sorted(completed_ids)
    nb.metadata["article_runner"]["updated_at"] = datetime.now().astimezone().isoformat()
    path.parent.mkdir(parents=True, exist_ok=True)
    nbformat.write(nb, path)


def execute_notebook(notebook_path: Path, approach_name: str, *, strata: str, data_dir: Path | None, resume: bool) -> Path:
    root = repo_root()
    source_nb = nbformat.read(notebook_path, as_version=4)
    output = root / "outputs" / "notebooks" / f"{approach_name}.executed.ipynb"
    log_path = root / "outputs" / "logs" / f"{approach_name}.log"
    logger = Logger(log_path)

    completed: set[str] = set()
    if resume and output.exists():
        old_nb = nbformat.read(output, as_version=4)
        completed.update(old_nb.metadata.get("article_runner", {}).get("completed_cell_ids", []))
        logger.write(f"RESUME requested | {len(completed)} previously completed cell(s) found")

    env = {"MICE_TCR_REPO": str(root), "MICE_TCR_STRATA": strata}
    if data_dir is not None:
        env["MICE_TCR_DATA_DIR"] = str(data_dir.resolve())

    code_cells = [c for c in source_nb.cells if c.cell_type == "code"]
    total = len(code_cells)
    ordinal = {c.id: i + 1 for i, c in enumerate(code_cells)}
    logger.write(f"START {approach_name} | notebook={notebook_path.relative_to(root)} | code_cells={total} | strata={strata}")

    with temporary_env(env):
        client = NotebookClient(source_nb, timeout=None, kernel_name="python3",
                                resources={"metadata": {"path": str(root)}}, allow_errors=False)
        try:
            with client.setup_kernel(cwd=str(root)):
                for idx, cell in enumerate(source_nb.cells):
                    if cell.cell_type != "code":
                        continue
                    pos = ordinal[cell.id]
                    remaining = total - pos
                    bootstrap = "bootstrap" in set(cell.metadata.get("tags", []))
                    if resume and cell.id in completed and not bootstrap:
                        logger.write(f"[CELL {pos:03d}/{total:03d}] SKIP completed checkpoint | {remaining} remaining")
                        continue

                    label = f"[CELL {pos:03d}/{total:03d}]"
                    started = time.monotonic()
                    logger.write(f"{label} START | {remaining} remaining after this cell")
                    stop = threading.Event()
                    beat = threading.Thread(target=_heartbeat, args=(logger, stop, label, started), daemon=True)
                    beat.start()
                    try:
                        client.execute_cell(cell, idx, store_history=True)
                    except BaseException:
                        stop.set(); beat.join(timeout=1)
                        logger.write(f"{label} FAILED after {time.monotonic() - started:.1f}s")
                        _save_checkpoint(source_nb, output, completed)
                        raise
                    else:
                        stop.set(); beat.join(timeout=1)
                        completed.add(cell.id)
                        logger.write(f"{label} DONE in {time.monotonic() - started:.1f}s | {total - len(completed)} code cell(s) not yet checkpointed")
                        _save_checkpoint(source_nb, output, completed)
        finally:
            logger.write(f"CHECKPOINT {output.relative_to(root)} | completed={len(completed)}/{total}")
            logger.close()
    return output


def resolve_sequence(value: str) -> list[str]:
    value = ALIASES.get(value.lower(), value.lower())
    if value == "all":
        return ["1", "2", "3", "4"]
    parts = [ALIASES.get(x.strip().lower(), x.strip().lower()) for x in value.split(",")]
    bad = [x for x in parts if x not in APPROACHES]
    if bad:
        raise SystemExit(f"Unknown approach value(s): {bad}. Use 1,2,3,4,all.")
    return parts


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one analytical approach or the complete article pipeline.")
    parser.add_argument("--approach", default="all", help="1, 2, 3, 4/compare, comma-separated values, or all (default).")
    parser.add_argument("--strata", default="all", help="all or comma-separated canonical strata.")
    parser.add_argument("--data-dir", type=Path, help="Directory containing clean_clonotypes_aaV.parquet.")
    parser.add_argument("--resume", action="store_true",
                        help="Re-run bootstrap cells, then skip completed non-bootstrap cells. Tables and figures remain durable checkpoints.")
    args = parser.parse_args()

    root = repo_root()
    for key in resolve_sequence(args.approach):
        name, relative = APPROACHES[key]
        execute_notebook(root / relative, name, strata=args.strata, data_dir=args.data_dir, resume=args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
