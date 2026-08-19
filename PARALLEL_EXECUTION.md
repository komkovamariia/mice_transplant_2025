# Parallel execution for article notebooks

The article notebooks use `article_runtime.py` to consume all CPUs assigned to the current process without oversubscribing nested numerical libraries.

## CPU selection

At notebook startup `configure_runtime()` chooses the usable worker count in this order:

1. explicit `ARTICLE_N_JOBS`, when set;
2. Linux CPU affinity (`sched_getaffinity`), which reflects cpuset/cgroup restrictions on most HPC systems;
3. Slurm `SLURM_CPUS_PER_TASK` / `SLURM_CPUS_ON_NODE`;
4. `os.cpu_count()` only when no scheduler/affinity limit is visible.

The most restrictive detected scheduler/affinity value is used. Thus the code does not consume CPUs outside the allocation.

To verify the allocation from a shell:

```bash
python - <<'PY'
import os
from article_runtime import available_cpus
print('available_cpus =', available_cpus())
print('affinity       =', len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else 'n/a')
print('SLURM_CPUS_PER_TASK =', os.environ.get('SLURM_CPUS_PER_TASK'))
PY
```

`ARTICLE_N_JOBS` may be used to make a run deliberately smaller, for example:

```bash
ARTICLE_N_JOBS=8 jupyter nbconvert --to notebook --execute mirpy_analysis.ipynb \
  --ExecutePreprocessor.kernel_name=python3 --ExecutePreprocessor.timeout=-1
```

Do not set `ARTICLE_N_JOBS` above the scheduler allocation.

## What is parallelized

### `venn_original.ipynb`

`repseq` sample/batch calculations are dispatched through a fixed-size joblib/loky pool using all allocated CPUs. Each process is restricted to one inner BLAS/OpenMP thread, preventing `N processes × N BLAS threads` oversubscription. Input/result order is preserved.

### `mirpy_analysis.ipynb`

- `TCREmp` receives `threads=N_WORKERS`, so sequence embedding uses the full CPU allocation internally.
- BLAS/OpenMP/NumExpr/Polars thread pools receive the same scheduler-aware limit.
- Exact `neighbor_enrichment(..., backend='kdtree')` is retained. No approximate-nearest-neighbour shortcut is introduced.
- Large outer embedding/projection chunks remain sequential intentionally: the operation inside each chunk is already parallel, while concurrent chunks would duplicate large embedding buffers and can exceed RAM.
- The global 60k fit, PCA basis and union coordinates are computed once and cached for downstream MMD, witness and contrast blocks rather than being fitted repeatedly.
- PERMANOVA permutations are generated deterministically from a fixed seed, split into batches and evaluated across all allocated workers. Parallel scheduling cannot alter the permutation schedule.

### `study2_alloreactivity/alloreactivity_study.ipynb`

- `repseq` batch calculations use the same scheduler-aware process pool.
- Independent Hill-number file calculations are parallelized with an ordered thread map.
- The stochastic diversity calculation uses a fixed seed for reproducibility.
- Large clustering jobs, when regenerated rather than loaded from saved tables, use the configured `repseq` worker pool.

## Scientific invariants

Parallelization does not change the clonotype definition, UMI thresholds, downsampling depth, group contrasts, exact-density backend, multiple-testing thresholds or ranking criteria. Large memory-bound stages are batched instead of naively duplicated across processes.

One separate correctness repair was made while parallelizing `mirpy_analysis.ipynb`: the historical PERMANOVA helper ignored the permuted-label argument and therefore recomputed the observed grouping on every permutation. The replacement performs the intended label permutation. Consequently, PERMANOVA p-values from the corrected notebook must be treated as newly calculated results rather than expected to reproduce the erroneous historical p-values.
