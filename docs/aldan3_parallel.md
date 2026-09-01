# Bounded parallel execution on Aldan-3

The combined sensitivity experiment uses a **24-CPU budget** in two non-overlapping
phases: one baseline with up to 24 workers, then at most three dose notebooks with
8 workers each. A single finalizer runs after every dose succeeds. Slurm assigns and
binds the CPUs; no hard-coded CPU IDs, GPU reservations or exclusive nodes are used.
The budget applies to this experiment, not to other jobs already running in your account.

## Cluster settings and first measurement

The supplied Aldan-3 documentation screenshots specify `short` (2 hours), `medium`
(16 hours) and `--constraint=hpc` for CPU/HPC resources. Their live web pages require
GitLab authentication in the development environment. Confirm current availability
from the cluster before submission:

```bash
cd ~/mice_transplant_env_test
conda activate mice-transplant-2025
python scripts/slurm_spike.py diagnose
```

This performs four read-only scheduler queries. It does not submit or cancel jobs.
Memory consumption of the real study has not been measured here. In the commands
below, **64G is an example reservation per job**, not a measured requirement. Three
concurrent dose jobs would reserve up to 192G in total. Use the existing baseline's
Slurm accounting, or a new pilot, to choose this reservation with headroom. Too small
a reservation can terminate the job; a large reservation may increase queue time.

Prepare a fresh experiment, without submitting anything:

```bash
python scripts/slurm_spike.py prepare \
  --run-id combined_24cpu_01 \
  --cpu-budget 24 --baseline-cpus 24 --cpus 8 --max-parallel 3 \
  --mem 64G --partition short --time 02:00:00
```

The default grid contains sixteen **family doses**:

```text
0.00001%, 0.00002%, 0.00005%, 0.0001%, 0.0002%, 0.0005%,
0.001%, 0.002%, 0.005%, 0.01%, 0.02%, 0.05%, 0.1%, 0.2%, 0.5%, 1%
```

`--fractions` accepts unit fractions, so `1e-7` means `0.00001%` and `0.01` means
`1%`. Ten observed aaV and seed 1031 are used by default. Each dose starts from the
same baseline. See [the protocol](spike_in.md) for selection, integer rounding,
hypotheses and limitations. Very small doses can round to zero UMI.

The plan freezes input-path overrides and the current Python interpreter. Activate
the intended Conda environment before preparation. Paths can also be supplied with
`--metadata-csv`, `--clonoset-index`, `--working-dir` and `--data-dir`. Source fingerprints
are checked at preparation, worker startup and completion; keep this checkout and its
software environment unchanged until the experiment finishes.

## Baseline, then dose array

Submit the baseline pilot:

```bash
python scripts/slurm_spike.py submit \
  logs/slurm/combined_24cpu_01/plan.json --stage baseline
```

Read its job ID from the printed output or `logs/slurm/combined_24cpu_01/jobs.json`.
Use that ID in the following commands:

```bash
squeue -u "$USER"
sstat -j JOB_ID.0 --format=JobID,AveCPU,MaxRSS,AveRSS
sacct -j JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
tail -f logs/spike_in/combined_24cpu_01/baseline.log
```

`sstat` is useful while running if accounting is enabled; `sacct` includes completed
steps. Inspect the `.0` srun step as well as `.batch`. With process workers, MaxRSS
may represent the largest task rather than a complete simultaneous sum; use it with
the cluster's available memory accounting and leave headroom. `TotalCPU / (Elapsed *
AllocCPUS)` estimates CPU utilization when those fields cover all task processes.
Cell timings identify whether input/count construction, Fisher, edgeR or figures
dominate. A serial cell can legitimately use one CPU out of the allocation.

After the pilot succeeds and its memory reservation is suitable, submit remaining jobs:

```bash
python scripts/slurm_spike.py submit \
  logs/slurm/combined_24cpu_01/plan.json --stage remaining
```

This reuses the completed baseline. It submits `--array=1-16%3` with
`--dependency=afterok:BASELINE_ID`, then a finalizer with `afterok:ARRAY_ID`. If resources
are already measured, `--stage all` submits all three scheduler jobs in one call with
the same dependencies. Repeating submission does not duplicate recorded job IDs.
A failed baseline prevents all doses; a failed array element prevents publication
of an incomplete summary. Slurm's invalid-dependency cancellation is enabled.

Monitor specific jobs without polling the scheduler rapidly:

```bash
squeue -r -u "$USER"
tail -f logs/spike_in/combined_24cpu_01/dose_01.log
```

To cancel this experiment, pass only its baseline, array and finalizer job IDs from
`jobs.json` to `scancel`. Do not use account-wide cancellation. Failed case directories
are preserved for diagnosis. Automatic requeue and overwrite are disabled. After a
code or experiment-setting change, prepare a fresh run ID. Older sequential runs are
not automatically adopted because they lack the staged completion/provenance markers.
Let a currently running old experiment finish or stop its specific job deliberately
before starting another 24-CPU experiment; the runner does not cancel it for you.

## What is parallel, and what preserves biological identity

| Level | Scheduling | Invariant |
|---|---|---|
| Baseline | One notebook, up to 24 workers | Unmodified combined TRA/TRB inputs |
| Dose experiments | Up to 3 independent processes/kernels, 8 workers each | One frozen family and baseline; independent g1 additions |
| repseq operations | Ordered process workers, capped at allocated CPUs and task count | Original aaV keys and sample order |
| Fisher tests | Ordered batches of independent exact two-sided tests | Same contingency tables; one global BH correction after gathering |
| edgeR | One model per notebook, no outer split across features/mice | Same TMM, filter, dispersion, QL fit, contrast and mouse-level pooling |
| Summary and ZIP | One finalizer after all successful doses | Canonical dose order, no concurrent archive writes |

Native BLAS/OpenMP threads in process workers are capped at one. Threaded maps cap
their inner numerical pools too. Explicit worker requests cannot exceed the allocation
or CPU affinity. Outside Slurm, automatic execution defaults to one worker; inside a
proper task allocation it uses `SLURM_CPUS_PER_TASK`. Requesting 24 CPUs does not make
every Python/R operation parallel or imply a 24-fold speedup. More process workers also
require more memory, so compare measured wall time and memory before increasing counts.

Each dose loads checksummed baseline snapshots; g1 receives its declared addition.
Controls and TRB are unchanged. Source/configuration/selection checks reject mismatches.
Exclusive case claims prevent duplicate jobs overwriting outputs. Every case has its
own executed notebook, log, result directory and two TRA/TRB V-segment heatmaps.
There is no pooling of doses, mice or statistical tests across notebooks.

The outputs remain under `results/01_spike_in/<run_id>/`,
`audit_runs/spike_in/<run_id>/`, `logs/spike_in/<run_id>/` and
`figures/01_spike_in/<run_id>/`. Scheduler plans, job IDs and stdout logs are under
`logs/slurm/<run_id>/`. The finalizer writes the hypothesis summary, sensitivity
bounds and `figures/01_spike_in/<run_id>.zip` only after all cases pass validation.

This scheduler entry point covers the combined baseline and independent spike-in
doses. The standard nine-stratum runner and dependent later approaches retain their
existing execution order; do not launch copies against shared output paths.

## Sources

- [Aldan-3 Slurm guide](https://aldan3.pages.itm-rsmu.ru/docs/software/slurm/), supported by the supplied screenshots.
- [Aldan-3 hardware](https://aldan3.pages.itm-rsmu.ru/docs/hardware/); hardware totals were not independently verified.
- [Slurm job arrays](https://slurm.schedmd.com/job_array.html): concurrency limits and whole-array dependencies.
- [Slurm CPU management](https://slurm.schedmd.com/cpu_management.html): allocation and binding.
- [Slurm sbatch](https://slurm.schedmd.com/sbatch.html): CPU, memory, time and dependency options.
