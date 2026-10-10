# Running the analysis

## Update an HPC checkout

```bash
cd ~/mice_transplant_env_test
git status --short
git switch main
git pull --ff-only origin main
conda env update -n mice-transplant-2025 -f environment.yml --prune
conda activate mice-transplant-2025
python scripts/validate_repository.py
```

Commit or move local source changes before switching branches. Generated outputs are ignored by Git, so old untracked files can remain after a pull.

## Approach 1

Run all nine strata:

```bash
python scripts/run_analysis.py --approach 1
```

Run selected strata:

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined

python scripts/run_analysis.py --approach 1 \
  --strata all_combined,cd4_combined,cd8_combined
```

Available keys:

```text
all_combined
cd4_thymus
cd8_thymus
cd4_spleen
cd8_spleen
cd4_combined
cd8_combined
thymus_combined
spleen_combined
```

The structural branch uses functional clonotypes and a 15,000-UMI library depth. Libraries below that depth are left out of structural set calculations. The representative structural seed for `all_combined`, `cd4_combined` and `cd8_combined` is chosen from 100 baseline rarefactions. Count-based inference keeps functional raw UMI counts and pools repeated samples within treatment group and mouse. The exact model and score definitions are in [docs/methods.md](docs/methods.md).

## Input paths

Default Approach 1 inputs:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
```

Override them on the command line:

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined \
  --metadata-csv /projects/mice_transplant_2025/metadata_mice_transplant.csv \
  --clonoset-index /projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv \
  --working-dir /projects/mice_transplant_2025/test_run
```

The metadata table must contain `chain`, `sample_no`, `sample_id`, `group_no`, `mouse_no`, `source` and `subtype`. The clonoset index must contain `sample_id`, `filename` and `chain`. Approach 1 reads the MiXCR exports referenced by the index and does not require the derived Parquet table.

See [config/example.env](config/example.env) for environment-variable overrides and [docs/notebook_input_provenance.md](docs/notebook_input_provenance.md) for the input map.

## Logs, resume and output checks

Follow a running notebook with:

```bash
tail -f logs/00_all_combined.log
tail -f logs/01_cd4_thymus.log
```

To rerun after interruption:

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined --resume
```

Resume starts a new kernel and replays cells from the beginning. Existing checkpoint outputs are kept for inspection, but Python and R objects are rebuilt.

A successful Approach 1 stratum contains at least:

```text
results/01_set_count/<stratum>/run_summary.json
results/01_set_count/<stratum>/trav_ranking.csv
results/01_set_count/<stratum>/distinctive_trav.csv
results/01_set_count/<stratum>/tra_v_retention.csv
results/01_set_count/<stratum>/figure_manifest.csv
```

The representative-seed strata also write:

```text
representative_downsampling_seed.json
tra_rarefaction_seed_summary.csv
tra_rarefaction_v_profiles.csv
tra_top10_only_g1_by_seed.csv
tra_top10_only_g1_summary.csv
tra_rarefaction_top10_stability.csv
trb_rarefaction_top10_stability.csv
tra_balanced_mouse_top10_stability.csv
trb_balanced_mouse_top10_stability.csv
```

`run_summary.json` should contain:

```text
edgeR quasi-likelihood model fitted successfully.
```

Figures are stored under `figures/01_set_count/<stratum>/TRA/` and `TRB/`. The main V-segment files include:

```text
TRA_g1_v_segment_representation.png
TRA_g1_v_segment_retention.png
TRA_g1_v_segment_localization.png
TRB_g1_v_segment_representation.png
TRB_g1_v_segment_retention.png
TRB_g1_v_segment_localization.png
```

For the representative-seed strata, the stability outputs also include `TRA_g1_rarefaction_stability.png`, `TRB_g1_rarefaction_stability.png`, `TRA_g1_balanced_mouse_stability.png`, `TRB_g1_balanced_mouse_stability.png` and `TRA_g1_top10_only_g1_stability.png`.

A complete nine-stratum run writes `results/01_set_count/all_strata_distinctive_trav_summary.csv` and refreshes `figures/01_set_count.zip`.

## Spike-in experiment

Sequential run:

```bash
python scripts/run_analysis.py --mode spike-in --spike-run-id sensitivity_01
```

The sequential default uses fractions `0.0001,0.001,0.01`, corresponding to 0.01%, 0.1% and 1% of the original g1 TRA UMI library.

Custom grid:

```bash
python scripts/run_analysis.py --mode spike-in \
  --spike-run-id sensitivity_expanded \
  --spike-clones 10 \
  --spike-fractions 0.00001,0.0001,0.001,0.01,0.1 \
  --spike-seed 1031
```

An exact V segment can be requested with `--spike-v-gene`. Selection still checks the baseline eligibility rules. Use a new run ID for a new experiment.

Results are written under `results/01_spike_in/<run_id>/`; executed notebooks and logs use `audit_runs/spike_in/<run_id>/` and `logs/spike_in/<run_id>/`. The protocol and output definitions are in [docs/spike_in.md](docs/spike_in.md).

## Aldan-3 spike-in jobs

For the bounded 24-CPU Slurm setup, use [docs/aldan3_parallel.md](docs/aldan3_parallel.md). Preparation writes the plan and shell scripts only. Submission is a separate command.

```bash
python scripts/slurm_spike.py diagnose

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

Use the Aldan-3 guide before submitting the baseline or dose array.

## Direct notebook execution

`venn_original.ipynb` defaults to `all_combined`. A compartment can also be executed directly with nbconvert:

```bash
mkdir -p audit_runs
MICE_TCR_STRATUM=cd4_thymus jupyter nbconvert \
  --to notebook --execute venn_original.ipynb \
  --ExecutePreprocessor.timeout=-1 \
  --output-dir audit_runs \
  --output 01_cd4_thymus.executed.ipynb
```

The runner is preferred for normal analysis because it also handles logs, checkpoints, output verification and figure archives.

## Approaches 2 to 4

```bash
python scripts/run_analysis.py --approach 2 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 3 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 4
python scripts/run_analysis.py --approach all --data-dir /absolute/path/to/derived_data
```

Approaches 2 and 3 use `clean_clonotypes_aaV.parquet`. Approach 4 reads their standardized ranking tables together with Approach 1 outputs.

CPU allocation outside the spike-in scheduler is described in [docs/parallel_execution.md](docs/parallel_execution.md).

## Validation

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```
