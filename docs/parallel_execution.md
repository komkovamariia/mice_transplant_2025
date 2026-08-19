# Parallel execution

The active pipeline uses `src/runtime.py` to consume the CPU allocation visible to the current process without oversubscribing nested numerical libraries.

## CPU selection

`configure_runtime()` resolves the worker count from an explicit `ARTICLE_N_JOBS` override, Linux CPU affinity, Slurm allocation variables, and finally `os.cpu_count()` when no scheduler restriction is visible. When multiple allocation signals are present, the most restrictive positive value is used.

Check the effective allocation:

```bash
python - <<'PY'
import os
from src.runtime import available_cpus
print("available_cpus:", available_cpus())
print("os.cpu_count:", os.cpu_count())
if hasattr(os, "sched_getaffinity"):
    print("affinity:", len(os.sched_getaffinity(0)))
print("SLURM_CPUS_PER_TASK:", os.environ.get("SLURM_CPUS_PER_TASK"))
print("SLURM_CPUS_ON_NODE:", os.environ.get("SLURM_CPUS_ON_NODE"))
PY
```

Use `ARTICLE_N_JOBS` only to reduce a run within an existing scheduler allocation:

```bash
ARTICLE_N_JOBS=8 python scripts/run_analysis.py --approach 2
```

Do not set it above the allocated CPU count.

## Approach-specific policy

Approach 1 delegates R/edgeR numerical work to the configured native thread environment. Tissue-specific count tables remain sample resolved; combined strata pool thymus and spleen within mouse before model fitting. The six strata are analyzed serially so their R model objects do not compete for memory.

Approach 2 gives TCREmp the complete scheduler-aware CPU count. The expensive global sequence embedding is chunked in the outer loop while the native embedding and linear-algebra kernels use the allocated cores internally. Running several embedding chunks concurrently would duplicate multi-gigabyte buffers and is intentionally avoided. Local density enrichment retains the exact `kdtree` backend.

Approach 2 PERMANOVA uses `permanova_parallel()`. The complete permutation schedule is generated from a fixed seed in the parent process before worker dispatch. Worker completion order therefore cannot change the p-value.

Approach 3 is dominated by pandas aggregation and CDR3-neighborhood indexing. Abundance is normalized within mouse, and the six biological strata are executed as separate notebook cells. This provides durable checkpoints without duplicating the canonical input table across multiple process pools.

## Memory

CPU parallelism is not used as a reason to multiply memory-bound outer tasks. On HPC systems, request enough memory for the global embedding stage and monitor the Slurm step rather than a single Python PID:

```bash
sstat -j <jobid>.0 --format=JobID,AveCPU,AveRSS,MaxRSS,MaxVMSize
```

For process-level monitoring:

```bash
watch -n 10 "ps -u \$USER -o pid,ppid,etime,time,%cpu,%mem,rss,stat,cmd --sort=-%cpu | head -n 25"
```

## Scientific invariants

Scheduling does not change the aaV clonotype definition, biological strata, g1 versus g5+g6 contrast, FDR threshold, exact density backend, or deterministic permutation schedule.

The corrected PERMANOVA implementation is a scientific repair rather than a performance optimization. The historical helper ignored its permuted-label argument. Old p-values from that implementation must not be reused; only regenerated values from the active pipeline are valid.
