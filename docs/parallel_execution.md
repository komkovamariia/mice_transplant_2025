# Parallel execution

`src/runtime.py` limits worker and native-thread counts to the CPU allocation visible to the current process.

## CPU selection

`configure_runtime()` checks, in order of restriction, `ARTICLE_N_JOBS`, Linux CPU affinity, Slurm allocation variables and the host CPU count. The smallest positive limit wins.

Check the current allocation with:

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

`ARTICLE_N_JOBS` can reduce the worker count inside an existing allocation:

```bash
ARTICLE_N_JOBS=8 python scripts/run_analysis.py --approach 2
```

Do not set it above the scheduler allocation.

## Per-approach behavior

Approach 1 runs the nine strata serially. R/edgeR and native numerical libraries use the configured thread limits. Pooled strata combine repeated samples within mouse before count-based inference.

Approach 2 gives TCREmp the available scheduler-aware CPU count. The global embedding is chunked in the outer loop, while the embedding and linear-algebra kernels use the allocated cores internally. Local density enrichment keeps the `kdtree` backend.

PERMANOVA uses `permanova_parallel()`. The parent process creates the full permutation schedule from a fixed seed before work is split across processes, so worker completion order cannot change the p-value.

Approach 3 is mainly pandas aggregation and CDR3-neighborhood indexing. Its strata are executed as separate notebook cells rather than several concurrent process pools.

## Memory monitoring

For Slurm steps:

```bash
sstat -j <jobid>.0 --format=JobID,AveCPU,AveRSS,MaxRSS,MaxVMSize
```

For local process monitoring:

```bash
watch -n 10 "ps -u \$USER -o pid,ppid,etime,time,%cpu,%mem,rss,stat,cmd --sort=-%cpu | head -n 25"
```

Avoid running several memory-heavy outer tasks just because additional CPUs are available.

## Scientific behavior

Scheduling leaves the aaV definition, biological strata, g1 versus g5+g6 contrast, FDR threshold, density backend and permutation schedule unchanged.

The current PERMANOVA implementation uses the permuted labels when computing the null distribution. Older results from the helper that ignored those labels should be regenerated.
