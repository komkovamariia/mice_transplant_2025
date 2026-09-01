# Running the analysis

## Update an existing HPC checkout

```bash
cd ~/mice_transplant_env_test
git status --short
```

Commit or back up any local source changes before switching branches. Generated output
folders are ignored by Git. Once the working tree is ready:

```bash
git switch main
git pull --ff-only origin main
conda env update -n mice-transplant-2025 -f environment.yml --prune
conda activate mice-transplant-2025
python scripts/validate_repository.py
```

If validation reports a legacy `approaches/` directory, inspect it and move it outside
the checkout under a new backup name. Old untracked files can survive a Git update.
Do not restore an older source notebook over the current `venn_original.ipynb`.

## Default and selected analyses

```bash
# Default: Approach 1, all nine strata, combined analysis first.
python scripts/run_analysis.py
python scripts/run_analysis.py --approach 1

# Only CD4 + CD8, thymus + spleen.
python scripts/run_analysis.py --approach 1 --strata all_combined

# One of the original compartments, or a specified sequence.
python scripts/run_analysis.py --approach 1 --strata cd4_thymus
python scripts/run_analysis.py --approach 1 --strata all_combined,cd4_thymus,cd8_thymus
```

The keys are `all_combined`, `cd4_thymus`, `cd8_thymus`, `cd4_spleen`, `cd8_spleen`,
`cd4_combined`, `cd8_combined`, `thymus_combined` and `spleen_combined`.

CD4 + CD8 denotes the union of these sampled subsets. The combined analysis includes
both organs; counts for repeated samples from one mouse are summed before edgeR.
Each distinct treatment-group/mouse pair remains an independent analysis unit.

## Input paths

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined \
  --metadata-csv /projects/mice_transplant_2025/metadata_mice_transplant.csv \
  --clonoset-index /projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv \
  --working-dir /projects/mice_transplant_2025/test_run
```

Metadata require `chain`, `sample_no`, `sample_id`, `group_no`, `mouse_no`, `source`
and `subtype`. The clonoset index requires `sample_id`, `filename` and `chain`.
Historical `alpha-N`/`beta-N` keys are matched before chain names are normalized.
All referenced MiXCR exports must be accessible on the analysis host.

Approach 1 does not read Parquet. See [config/example.env](config/example.env) for
persistent path overrides and [input provenance](docs/notebook_input_provenance.md)
for the distinction between source data and derived inputs.

## Progress, outputs and interruptions

```bash
tail -f logs/00_all_combined.log
tail -f logs/01_cd4_thymus.log
```

The logs report cell number, total cells, elapsed time and remaining cells. The combined
notebook is `audit_runs/00_all_combined.executed.ipynb`; the original CD4 thymus file is
`audit_runs/01_cd4_thymus.executed.ipynb`. Both retain completed cell outputs after failure.

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined --resume
```

`--resume` replays execution from the first cell in a new kernel. It preserves progress
visibility; it does not restore Python or R objects in memory. A standard run clears
old generated figures and tables for its selected strata. Copy outputs elsewhere first
if they need to be retained for comparison.

After successful execution, inspect `results/01_set_count/<stratum>/run_summary.json`,
`trav_ranking.csv`, `tra_v_retention.csv`, `tra_analysis_units.csv`, `provenance.json`
and `R_sessionInfo.txt`. `edger_status` must read:

```text
edgeR quasi-likelihood model fitted successfully.
```

Each figure listed in `figure_manifest.csv` must exist as PNG and PDF. Each stratum
has the two approved heatmaps, `TRA_g1_top10_bubble_genes_heatmap` and
`TRB_g1_top10_bubble_genes_heatmap`. These show at most ten threshold-eligible genes.

A complete nine-stratum run creates
`results/01_set_count/all_strata_distinctive_trav_summary.csv`.
Every successful Approach 1 invocation creates `figures/01_set_count.zip` from the
whole current figure folder. After a partial rerun, that ZIP also includes previously
generated figures for untouched strata. Run all nine strata together for a fresh,
consistent complete archive.

## Spike-in experiment

```bash
python scripts/run_analysis.py --mode spike-in --spike-run-id sensitivity_01
```

The default family has ten observed aaV clonotypes sharing a noncandidate TRAV.
The family receives 0.01%, 0.1% and 1% added UMI in each original g1 sample library.
The runner creates a fresh baseline and three full executed copies of the primary
notebook, all for `all_combined`. Controls and ordinary analysis outputs are preserved.

Custom family and dose grid:

```bash
python scripts/run_analysis.py --mode spike-in \
  --spike-run-id sensitivity_expanded \
  --spike-clones 10 \
  --spike-fractions 0.00001,0.0001,0.001,0.01,0.1 \
  --spike-seed 1031
```

An exact TRAV can be requested with `--spike-v-gene TRAV9-2`; this is an example,
not a claim that this segment is eligible in the study. Selection verifies eligibility.
`--spike-max-g1-fraction` defaults to `0.00001`, the maximum baseline frequency allowed
in any g1 sample. The runner never relaxes this threshold or invents sequences to
complete a family. If selection fails, inspect the baseline and explicitly choose a
smaller family or revised threshold.

Results are under `results/01_spike_in/<run_id>/`:

- `experiment.json` records the requested configuration.
- `selection.json` records the fixed TRAV, exact aaV keys and selection criteria.
- `baseline/count_snapshots/` contains hashed original TRA/TRB matrices.
- `baseline/` and `dose_01/`, `dose_02/`, ... contain full notebook outputs.
- Each dose has `spike_dose_by_sample.csv`, `spike_clonotype_recovery.csv` and `spike_metrics.json`.
- `spike_in_summary.csv` compares all hypotheses with baseline.
- `sensitivity_bounds.json` records observed pass/fail transitions.

Executed notebooks and logs use `audit_runs/spike_in/<run_id>/` and
`logs/spike_in/<run_id>/`. Figures use `figures/01_spike_in/<run_id>/` and a separate ZIP.
Use a fresh run identifier to rerun the experiment. `--resume` is supported for ordinary
runs only. See [the spike-in protocol](docs/spike_in.md) before interpreting limits.

## Direct notebook execution

Open `venn_original.ipynb` in the configured environment and run all cells from a fresh
kernel. Its direct default is `all_combined`. To execute a compartment via nbconvert:

```bash
mkdir -p audit_runs
MICE_TCR_STRATUM=cd4_thymus jupyter nbconvert \
  --to notebook --execute venn_original.ipynb \
  --ExecutePreprocessor.timeout=-1 \
  --output-dir audit_runs --output 01_cd4_thymus.executed.ipynb
```

The runner additionally manages cell logs, output replacement, the complete summary and
ZIP creation. Use the runner for the standard article output contract and spike-in mode.

## Additional approaches and validation

```bash
python scripts/run_analysis.py --approach 2 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 3 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 4
python scripts/run_analysis.py --approach all --data-dir /absolute/path/to/derived_data

python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

Approaches 2 and 3 consume an existing `clean_clonotypes_aaV.parquet`; Approach 4 reads
all upstream ranking tables. A fresh-start materialization command for that Parquet
interface has not yet been verified. Resource allocation is described in
[parallel execution](docs/parallel_execution.md).
