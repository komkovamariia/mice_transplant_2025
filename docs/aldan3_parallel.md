# Aldan-3 spike-in runs

The Slurm helper runs the combined spike-in experiment with a 24-CPU cap:

- baseline: up to 24 CPUs;
- dose array: up to three concurrent jobs with 8 CPUs each;
- finalizer: one CPU after all dose jobs succeed.

The helper does not request GPUs or exclusive nodes.

## Check the cluster

```bash
cd ~/mice_transplant_env_test
conda activate mice-transplant-2025
python scripts/slurm_spike.py diagnose
```

`diagnose` only reads Slurm state.

The examples below use the `short` partition, `--constraint=hpc` and 64G per analysis job. The 64G value is an example, not a measured requirement. Check `MaxRSS` from a baseline job before reusing it for the full array.

## Prepare a run

```bash
python scripts/slurm_spike.py prepare \
  --run-id combined_24cpu_01 \
  --cpu-budget 24 \
  --baseline-cpus 24 \
  --cpus 8 \
  --max-parallel 3 \
  --mem 64G \
  --partition short \
  --time 02:00:00
```

Preparation creates `logs/slurm/<run_id>/plan.json` and the stage scripts. It does not submit jobs.

The scheduler default contains 16 family fractions:

```text
0.00001%, 0.00002%, 0.00005%, 0.0001%, 0.0002%, 0.0005%,
0.001%, 0.002%, 0.005%, 0.01%, 0.02%, 0.05%, 0.1%, 0.2%, 0.5%, 1%
```

`--fractions` takes unit fractions, so `1e-7` is 0.00001% and `0.01` is 1%.

The plan stores the input-path overrides, Python interpreter and source fingerprint. Keep the checkout and environment unchanged after preparation. If code or settings change, prepare a new run ID.

## Submit the baseline

```bash
python scripts/slurm_spike.py submit \
  logs/slurm/combined_24cpu_01/plan.json \
  --stage baseline
```

The job ID is printed and stored in `logs/slurm/combined_24cpu_01/jobs.json`.

Useful checks:

```bash
squeue -u "$USER"
sstat -j JOB_ID.0 --format=JobID,AveCPU,MaxRSS,AveRSS
sacct -j JOB_ID --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,TotalCPU,MaxRSS,ReqMem
tail -f logs/spike_in/combined_24cpu_01/baseline.log
```

Use the baseline accounting to adjust the memory request if needed.

## Submit the dose array

After the baseline finishes successfully:

```bash
python scripts/slurm_spike.py submit \
  logs/slurm/combined_24cpu_01/plan.json \
  --stage remaining
```

This submits the dose array as `1-16%3` with an `afterok` dependency on the baseline. The finalizer depends on the full array.

If the resource settings are already established, all stages can be submitted at once:

```bash
python scripts/slurm_spike.py submit \
  logs/slurm/combined_24cpu_01/plan.json \
  --stage all
```

Submission is idempotent with respect to job IDs already recorded in `jobs.json`.

## Monitor or stop the run

```bash
squeue -r -u "$USER"
tail -f logs/spike_in/combined_24cpu_01/dose_01.log
```

To stop the experiment, pass the specific job IDs from `jobs.json` to `scancel`. Do not use account-wide cancellation.

Failed case directories are kept. Automatic requeue and output overwrite are disabled.

## Parallel layout

| Stage | Slurm layout | Analysis behavior |
|---|---|---|
| Baseline | one job, up to 24 CPUs | unmodified combined TRA/TRB baseline |
| Doses | up to 3 jobs at once, 8 CPUs each | one frozen target family, independent g1 additions |
| repseq work | process workers capped by the allocation | original aaV keys and sample order |
| Fisher | ordered batches | same contingency tables and one global BH correction |
| edgeR | one model per notebook | same TMM, filtering, dispersion, QL fit and mouse pooling |
| Finalizer | one job | ordered summary and ZIP creation |

Native BLAS/OpenMP threads inside process workers are capped to avoid nested oversubscription. `SLURM_CPUS_PER_TASK` is checked against the prepared plan before a worker starts.

Each dose loads the checksummed baseline matrices. Only g1 TRA receives the spike-in. Controls and TRB are unchanged.

## Output paths

```text
results/01_spike_in/<run_id>/
figures/01_spike_in/<run_id>/
audit_runs/spike_in/<run_id>/
logs/spike_in/<run_id>/
logs/slurm/<run_id>/
```

The finalizer writes the run summary, sensitivity bounds and `figures/01_spike_in/<run_id>.zip` after every case passes validation.

For the biological definitions and dose calculation, see [spike_in.md](spike_in.md).
