"""Scheduler-aware CPU limits and parallel helpers."""

from __future__ import annotations

import importlib
import os
import re
from collections.abc import Iterable

_DEFAULT_WORKERS: int | None = None
_REPSEQ_DEFAULT_WORKERS: int | None = None


def _positive_int(value: str | None) -> int | None:
    if not value:
        return None
    match = re.search(r"\d+", str(value))
    if not match:
        return None
    number = int(match.group())
    return number if number > 0 else None


def cpu_limit() -> int:
    """Return the CPU ceiling visible to this process."""
    limits: list[int] = []

    if hasattr(os, "sched_getaffinity"):
        try:
            affinity = len(os.sched_getaffinity(0))
            if affinity > 0:
                limits.append(affinity)
        except (OSError, AttributeError):
            pass

    for key in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        limit = _positive_int(os.environ.get(key))
        if limit is not None:
            limits.append(limit)

    limits.append(max(1, os.cpu_count() or 1))
    # a slurm job without a per-task limit falls back to one worker
    if os.environ.get("SLURM_JOB_ID") and not _positive_int(os.environ.get("SLURM_CPUS_PER_TASK")):
        limits.append(1)
    return max(1, min(limits))


def available_cpus() -> int:
    """Return the requested worker count within the current CPU limit."""
    requested = _positive_int(os.environ.get("ARTICLE_N_JOBS"))
    if requested is None:
        requested = _positive_int(os.environ.get("SLURM_CPUS_PER_TASK")) or 1
    return min(requested, cpu_limit())


def configure_runtime(n_jobs: int | None = None, *, verbose: bool = True) -> int:
    """Set thread limits before numerical libraries are imported."""
    global _DEFAULT_WORKERS

    workers = int(n_jobs) if n_jobs is not None else available_cpus()
    if workers < 1:
        raise ValueError("n_jobs must be >= 1")
    workers = min(workers, cpu_limit())
    _DEFAULT_WORKERS = workers

    thread_env = {
        "OMP_NUM_THREADS": workers,
        "OPENBLAS_NUM_THREADS": workers,
        "MKL_NUM_THREADS": workers,
        "NUMEXPR_NUM_THREADS": workers,
        "VECLIB_MAXIMUM_THREADS": workers,
        "POLARS_MAX_THREADS": workers,
    }
    for key, value in thread_env.items():
        os.environ[key] = str(value)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["MKL_DYNAMIC"] = "FALSE"
    os.environ["OMP_MAX_ACTIVE_LEVELS"] = "1"

    if verbose:
        affinity = None
        if hasattr(os, "sched_getaffinity"):
            try:
                affinity = len(os.sched_getaffinity(0))
            except OSError:
                pass
        print(
            f"Parallel runtime: {workers} CPU worker(s)"
            + (f" | affinity={affinity}" if affinity is not None else "")
            + (
                f" | SLURM_CPUS_PER_TASK={os.environ.get('SLURM_CPUS_PER_TASK')}"
                if os.environ.get("SLURM_CPUS_PER_TASK")
                else ""
            )
        )
    return workers


def _resolve_jobs(n_jobs: int | None, n_tasks: int | None = None) -> int:
    workers = int(n_jobs) if n_jobs is not None else (
        _DEFAULT_WORKERS or available_cpus()
    )
    if workers < 1:
        raise ValueError("n_jobs must be >= 1")
    workers = min(workers, cpu_limit())
    if n_tasks is not None:
        workers = min(workers, max(1, n_tasks))
    return max(1, workers)


def _repseq_parallel_runner(
    function,
    tasks,
    program_name,
    object_name="tasks",
    verbose=True,
    cpu=None,
):
    """Run repseq sample tasks within the configured CPU limit."""
    tasks = list(tasks)
    n_tasks = len(tasks)
    if n_tasks == 0:
        return []

    workers = _resolve_jobs(cpu or _REPSEQ_DEFAULT_WORKERS, n_tasks)
    if verbose:
        print(
            f"{program_name}: {n_tasks} {object_name}; "
            f"using {workers} process worker(s)"
        )

    if workers == 1:
        return [function(task) for task in tasks]

    from joblib import Parallel, delayed, parallel_config

    with parallel_config(backend="loky", n_jobs=workers, inner_max_num_threads=1):
        return Parallel()(delayed(function)(task) for task in tasks)


def configure_repseq_parallelism(
    n_jobs: int | None = None, *, verbose: bool = True
) -> int:
    """Apply the configured worker limit to repseq modules used here."""
    global _REPSEQ_DEFAULT_WORKERS
    workers = _resolve_jobs(n_jobs)
    _REPSEQ_DEFAULT_WORKERS = workers

    module_names = (
        "repseq.common_functions",
        "repseq.stats",
        "repseq.intersections",
        "repseq.clustering",
        "repseq.diffexp",
        "repseq.pgen_calculation",
    )
    patched: list[str] = []
    for module_name in module_names:
        try:
            module = importlib.import_module(module_name)
        except ImportError:
            continue
        if hasattr(module, "run_parallel_calculation"):
            module.run_parallel_calculation = _repseq_parallel_runner
            patched.append(module_name)

    if verbose:
        print(f"repseq parallelism: {workers} worker(s); patched {', '.join(patched)}")
    return workers


def parallel_map(
    function,
    items: Iterable,
    *,
    n_jobs: int | None = None,
    prefer: str = "threads",
    batch_size="auto",
):
    """Map independent tasks in input order."""
    items = list(items)
    if not items:
        return []
    workers = _resolve_jobs(n_jobs, len(items))
    if workers == 1:
        return [function(item) for item in items]

    from joblib import Parallel, delayed, parallel_config

    if prefer == "processes":
        with parallel_config(backend="loky", n_jobs=workers, inner_max_num_threads=1):
            return Parallel(batch_size=batch_size)(
                delayed(function)(item) for item in items
            )

    if prefer != "threads":
        raise ValueError("prefer must be 'threads' or 'processes'")

    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=1):
        return Parallel(n_jobs=workers, prefer="threads", batch_size=batch_size)(
            delayed(function)(item) for item in items
        )


def permanova_parallel(
    distance_matrix,
    labels,
    n_perm: int = 9999,
    seed: int = 0,
    n_jobs: int | None = None,
):
    """Run one-factor PERMANOVA with a fixed permutation schedule."""
    import numpy as np

    distance_matrix = np.asarray(distance_matrix, dtype=np.float64)
    labels = np.asarray(labels)
    if (
        distance_matrix.ndim != 2
        or distance_matrix.shape[0] != distance_matrix.shape[1]
    ):
        raise ValueError("distance_matrix must be square")
    if distance_matrix.shape[0] != labels.shape[0]:
        raise ValueError("labels length must match distance_matrix")

    sample_count = len(labels)
    groups, group_codes = np.unique(labels, return_inverse=True)
    group_count = len(groups)
    if group_count < 2 or sample_count <= group_count:
        raise ValueError(
            "PERMANOVA requires at least two groups and residual degrees of freedom"
        )

    squared_distances = distance_matrix * distance_matrix
    upper_triangle = np.triu_indices(sample_count, 1)
    total_sum_squares = (
        squared_distances[upper_triangle].sum() / sample_count
    )

    def pseudo_f(codes_for_groups) -> float:
        within_sum_squares = 0.0
        for group_index in range(group_count):
            member_indices = np.flatnonzero(codes_for_groups == group_index)
            if len(member_indices) < 2:
                continue
            within_group = squared_distances[
                np.ix_(member_indices, member_indices)
            ]
            upper_group = np.triu_indices(len(member_indices), 1)
            within_sum_squares += (
                within_group[upper_group].sum() / len(member_indices)
            )

        between_sum_squares = total_sum_squares - within_sum_squares
        return (
            between_sum_squares / (group_count - 1)
        ) / (
            within_sum_squares / (sample_count - group_count)
        )

    observed_f = pseudo_f(group_codes)
    if n_perm <= 0:
        return observed_f, 1.0

    # generate permutations before dispatch so worker count cannot change the rng stream
    random_generator = np.random.default_rng(seed)
    permutation_codes = np.empty(
        (n_perm, sample_count),
        dtype=np.int16 if group_count < 32768 else np.int32,
    )
    for permutation_index in range(n_perm):
        permutation_codes[permutation_index] = random_generator.permutation(
            group_codes
        )

    workers = _resolve_jobs(n_jobs, n_perm)
    batch_count = min(n_perm, max(workers * 4, 1))
    batches = [
        batch
        for batch in np.array_split(permutation_codes, batch_count)
        if len(batch)
    ]

    def evaluate_batch(batch) -> int:
        return sum(pseudo_f(row) >= observed_f for row in batch)

    if workers == 1:
        exceedances = sum(evaluate_batch(batch) for batch in batches)
    else:
        from joblib import Parallel, delayed, parallel_config

        with parallel_config(
            backend="loky",
            n_jobs=workers,
            inner_max_num_threads=1,
        ):
            exceedances = sum(
                Parallel()(
                    delayed(evaluate_batch)(batch)
                    for batch in batches
                )
            )

    return observed_f, (exceedances + 1) / (n_perm + 1)
