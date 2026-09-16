"""Runtime helpers for reproducible, scheduler-aware parallel execution.

The article notebooks run both libraries that parallelize internally (mirpy/seqtree,
SciPy cKDTree, BLAS) and repseq routines that parallelize across samples.  This module
uses the requested allocation while avoiding nested process x BLAS oversubscription.
Outside a scheduler it defaults to one worker; explicit requests remain affinity-bound.

Scientific invariants:
- scheduling never changes clonotype definitions, filters, thresholds or contrasts;
- ordered maps preserve input order;
- PERMANOVA permutations are generated deterministically in the parent process before
  they are evaluated in parallel;
- exact mirpy density (kdtree) remains exact; this module does not switch to ANN.
"""

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
    m = re.search(r"\d+", str(value))
    if not m:
        return None
    n = int(m.group())
    return n if n > 0 else None


def cpu_limit() -> int:
    """Hard ceiling from affinity, scheduler allocation and host CPU count."""
    limits: list[int] = []

    if hasattr(os, "sched_getaffinity"):
        try:
            affinity = len(os.sched_getaffinity(0))
            if affinity > 0:
                limits.append(affinity)
        except (OSError, AttributeError):
            pass

    for key in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        n = _positive_int(os.environ.get(key))
        if n is not None:
            limits.append(n)

    limits.append(max(1, os.cpu_count() or 1))
    # A job without a per-task allocation must not use the entire shared node
    if os.environ.get("SLURM_JOB_ID") and not _positive_int(os.environ.get("SLURM_CPUS_PER_TASK")):
        limits.append(1)
    return max(1, min(limits))


def available_cpus() -> int:
    """Return the requested workers, always capped by the current allocation."""
    requested = _positive_int(os.environ.get("ARTICLE_N_JOBS"))
    if requested is None:
        requested = _positive_int(os.environ.get("SLURM_CPUS_PER_TASK")) or 1
    return min(requested, cpu_limit())


def configure_runtime(n_jobs: int | None = None, *, verbose: bool = True) -> int:
    """Configure native thread pools before NumPy/Polars/mirpy are imported.

    Call this at the very beginning of a fresh kernel.  The variables cover BLAS,
    OpenMP, NumExpr and Polars.  ``TCREmp`` still receives the same returned value
    explicitly, so seqtree uses exactly the CPUs granted by the scheduler.
    """
    global _DEFAULT_WORKERS

    n = int(n_jobs) if n_jobs is not None else available_cpus()
    if n < 1:
        raise ValueError("n_jobs must be >= 1")
    n = min(n, cpu_limit())
    _DEFAULT_WORKERS = n

    thread_env = {
        "OMP_NUM_THREADS": n,
        "OPENBLAS_NUM_THREADS": n,
        "MKL_NUM_THREADS": n,
        "NUMEXPR_NUM_THREADS": n,
        "VECLIB_MAXIMUM_THREADS": n,
        "POLARS_MAX_THREADS": n,
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
            f"Parallel runtime: {n} CPU worker(s)"
            + (f" | affinity={affinity}" if affinity is not None else "")
            + (
                f" | SLURM_CPUS_PER_TASK={os.environ.get('SLURM_CPUS_PER_TASK')}"
                if os.environ.get("SLURM_CPUS_PER_TASK")
                else ""
            )
        )
    return n


def _resolve_jobs(n_jobs: int | None, n_tasks: int | None = None) -> int:
    n = int(n_jobs) if n_jobs is not None else (_DEFAULT_WORKERS or available_cpus())
    if n < 1:
        raise ValueError("n_jobs must be >= 1")
    n = min(n, cpu_limit())
    if n_tasks is not None:
        n = min(n, max(1, n_tasks))
    return max(1, n)


def _repseq_parallel_runner(
    function,
    tasks,
    program_name,
    object_name="tasks",
    verbose=True,
    cpu=None,
):
    """Drop-in replacement for repseq.common_functions.run_parallel_calculation.

    repseq's pinned implementation uses ``ProcessPoolExecutor(max_workers=None)``.
    On Python 3.12 that default may not use every allocated core and nested numerical
    libraries can oversubscribe each worker.  joblib/loky lets us request the exact
    scheduler allocation and cap each worker's native BLAS/OpenMP pool to one thread.
    Results are returned in input order, matching ``executor.map`` semantics.
    """
    tasks = list(tasks)
    n_tasks = len(tasks)
    if n_tasks == 0:
        return []

    n = _resolve_jobs(cpu or _REPSEQ_DEFAULT_WORKERS, n_tasks)
    if verbose:
        print(f"{program_name}: {n_tasks} {object_name}; using {n} process worker(s)")

    if n == 1:
        return [function(task) for task in tasks]

    from joblib import Parallel, delayed, parallel_config

    with parallel_config(backend="loky", n_jobs=n, inner_max_num_threads=1):
        return Parallel()(delayed(function)(task) for task in tasks)


def configure_repseq_parallelism(
    n_jobs: int | None = None, *, verbose: bool = True
) -> int:
    """Make repseq batch/sample calculations use all allocated CPUs safely.

    Several repseq modules import ``run_parallel_calculation`` into module scope, so
    patching only ``repseq.common_functions`` is insufficient.  We patch all modules
    used by this repository.  This changes scheduling only, not calculations.
    """
    global _REPSEQ_DEFAULT_WORKERS
    n = _resolve_jobs(n_jobs)
    _REPSEQ_DEFAULT_WORKERS = n

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
        print(f"repseq parallelism: {n} worker(s); patched {', '.join(patched)}")
    return n


def parallel_map(
    function,
    items: Iterable,
    *,
    n_jobs: int | None = None,
    prefer: str = "threads",
    batch_size="auto",
):
    """Ordered parallel map for independent tasks.

    ``prefer='threads'`` is appropriate for independent file I/O and NumPy kernels
    that release the GIL.  ``prefer='processes'`` uses loky with one native thread per
    process to prevent oversubscription.  The returned list has the same order as
    ``items``.
    """
    items = list(items)
    if not items:
        return []
    n = _resolve_jobs(n_jobs, len(items))
    if n == 1:
        return [function(item) for item in items]

    from joblib import Parallel, delayed, parallel_config

    if prefer == "processes":
        with parallel_config(backend="loky", n_jobs=n, inner_max_num_threads=1):
            return Parallel(batch_size=batch_size)(
                delayed(function)(item) for item in items
            )

    if prefer != "threads":
        raise ValueError("prefer must be 'threads' or 'processes'")

    from threadpoolctl import threadpool_limits

    with threadpool_limits(limits=1):
        return Parallel(n_jobs=n, prefer="threads", batch_size=batch_size)(
            delayed(function)(item) for item in items
        )


def permanova_parallel(
    D, labels, n_perm: int = 9999, seed: int = 0, n_jobs: int | None = None
):
    """Deterministic parallel one-factor PERMANOVA on a distance matrix.

    The permutation schedule is generated *once* in the parent process using the given
    seed and then split into batches.  Therefore worker completion order cannot change
    the p-value.  This also fixes a bug in the historical notebook implementation where
    the helper ignored its permuted-label argument and accidentally recomputed the
    observed statistic for every permutation.
    """
    import numpy as np

    D = np.asarray(D, dtype=np.float64)
    labels = np.asarray(labels)
    if D.ndim != 2 or D.shape[0] != D.shape[1]:
        raise ValueError("D must be a square distance matrix")
    if D.shape[0] != labels.shape[0]:
        raise ValueError("labels length must match D")

    n = len(labels)
    uniq, codes = np.unique(labels, return_inverse=True)
    a = len(uniq)
    if a < 2 or n <= a:
        raise ValueError(
            "PERMANOVA requires at least two groups and residual degrees of freedom"
        )

    D2 = D * D
    upper = np.triu_indices(n, 1)
    sst = D2[upper].sum() / n

    def pseudo_f(group_codes) -> float:
        ssw = 0.0
        for g in range(a):
            idx = np.flatnonzero(group_codes == g)
            if len(idx) < 2:
                continue
            sub = D2[np.ix_(idx, idx)]
            ssw += sub[np.triu_indices(len(idx), 1)].sum() / len(idx)
        ssb = sst - ssw
        return (ssb / (a - 1)) / (ssw / (n - a))

    F_obs = pseudo_f(codes)
    if n_perm <= 0:
        return F_obs, 1.0

    # generate permutations serially so the RNG stream is independent of worker count
    rng = np.random.default_rng(seed)
    perm_codes = np.empty((n_perm, n), dtype=np.int16 if a < 32768 else np.int32)
    for i in range(n_perm):
        perm_codes[i] = rng.permutation(codes)

    n_workers = _resolve_jobs(n_jobs, n_perm)
    n_batches = min(n_perm, max(n_workers * 4, 1))
    batches = [b for b in np.array_split(perm_codes, n_batches) if len(b)]

    def eval_batch(batch) -> int:
        return sum(pseudo_f(row) >= F_obs for row in batch)

    if n_workers == 1:
        exceed = sum(eval_batch(batch) for batch in batches)
    else:
        from joblib import Parallel, delayed, parallel_config

        with parallel_config(backend="loky", n_jobs=n_workers, inner_max_num_threads=1):
            exceed = sum(Parallel()(delayed(eval_batch)(batch) for batch in batches))

    return F_obs, (exceed + 1) / (n_perm + 1)
