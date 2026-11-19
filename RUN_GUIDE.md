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

To fit exploratory V-segment plot cutoffs for each chain and panel:

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined --v-thresholds adaptive
```

This command refreshes the retention, dumbbell and localization figures for
TRA and TRB, including the CD4/CD8 composite for `all_combined`. It reuses
the previous run's depth and structural seed and saves the input count
snapshots on first use. Preparation sorts aaV before rarefaction so future
updates are reproducible. Legacy runs that sampled in hash-dependent row order
can have a different rarefaction draw at the same seed. Use `--v-baseline-results PATH` to select the directory
containing the previous run's `run_summary.json`. This update skips the
count tests, MDS, spike-in and seed search. To rebuild the snapshots after
changing input data, use `python scripts/run_v_updates.py prepare` with
`--baseline-results` and any required input path overrides.

The Otsu criterion maximizes between-class variance. The fit splits the
log10 initial share first, then splits retention for the segments above
that share cutoff. Each side of a split needs at least two segments.
The initial-share split leaves at least four segments for the retention
split. The cutoffs describe the observed distribution.

The run starts each fit from the full set of downsampled libraries and omits
one mouse from g1, g5 or g6, including all its organs and cell subsets.
Each eligible mouse in those three groups is tried once per panel. A fit is
skipped if the mouse has no libraries in that panel or its removal would
empty its group. The group cutoffs are arithmetic means of valid fits in
percentage units. The dashed lines are the arithmetic mean of the three
group means, giving all three groups equal weight. If a group has no
valid fit, the final cutoffs are blank in the CSV. The bubble plot shows
the measured points without cutoffs or candidate labels. The corresponding
dumbbell and localization panels mark the selection as unavailable.

Read the mouse, group and final rows in
`results/01_set_count_adaptive/all_combined/v_thresholds.csv` for the fitted
cutoffs and the number of attempted and successful fits. Selection and ranking
use the full-panel data with the aggregated cutoffs.
Figures and the archive go to `figures/01_set_count_adaptive/` and
`figures/01_set_count_adaptive.zip`. Standard outputs and spike-in settings
use the fixed mode.

### Exploratory calibration of Equation 7

After preparing the `all_combined` inputs:

```bash
python scripts/run_v_updates.py evidence
python scripts/run_v_updates.py weights --top 10 --ridge 0.1
```

The evidence stage uses a fixed rarefaction seed and caches single and
paired mouse exclusions from g1, g5 and g6. All libraries for the omitted
mice are removed before fitting the raw-UMI edgeR and Fisher tests or
calculating structural features. Each pair shares one evidence fit.
The structural term counts unique aaV left after subtracting g2, g4, g5
and g6.

The weights stage fits two Otsu cutoffs directly on the training profile
in each fold. It ranks eligible segments by the fraction of the V
segment's g1 UMI burden in `g1 only` relative to g5 and g6. Ties use the
number of distinct aaV with that V segment, then the segment name. A fold
with no valid cutoffs supplies zero eligible segments. The cutoffs from
the full-cohort plots are descriptive and independent of this procedure.

For each shortlisted V segment, the response is the fraction of distinct
aaV with that V segment on the omitted mouse that also occur in the
training `g1 only` set. Responses are positive for omitted g1 mice and
negative for omitted g5 and g6 mice. Rows with no observations for that
V segment on the omitted mouse are excluded from the fit. The fit uses
nonnegative ridge regression with standardized predictors and equal total
weight for each group, equal weight per usable mouse within a group and
equal weight per observed top V segment within a mouse.

For validation, each outer mouse is excluded from both the training
features and the responses used to fit its weights. The inner fold
excludes a second mouse, computes its predictors on the remaining mice
and uses that second mouse only to measure the response. The final
selection and scoring for the outer fold use its training profile and
the weights fitted from the inner folds. The outer mouse provides the
response after selection and scoring. The `*_outer_status.csv` reports
why an outer fold was skipped. For example, if only three g1 mice have
usable data, omitting one g1 mouse in each of the outer and inner folds
can leave too few g1 replicates for edgeR and no g1 calibration responses.

The weights stage reuses evidence when changing `--top` or `--ridge`.
Prespecify these settings to interpret the outer results as validation.
The reports in `results/01_v_update_inputs/all_combined/eq7_weights/`
include outer-fold coefficients and scores. The coefficients in
`*_coefficients.csv` use the single-exclusion rows from the full cohort.
The final descriptive shortlist refits the Otsu cutoffs on the full
cohort using the same fixed rarefaction and the saved full-cohort count
test verdicts. Use the outer reports for held-out results.

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

Use the manuscript dose grid to generate the four-panel sensitivity figure for both chains:

```bash
python scripts/run_analysis.py --mode spike-in --strata all_combined \
  --spike-run-id article_tra --spike-chain TRA --spike-v-gene TRAV8N-2 --spike-fractions article
python scripts/run_analysis.py --mode spike-in --strata all_combined \
  --spike-run-id article_trb --spike-chain TRB --spike-fractions article
```

An exact V segment requested with `--spike-v-gene` must satisfy the baseline eligibility rules. With no explicit gene, selection uses the baseline and seed. Slurm preparation accepts `--chain TRB` and `--fractions article`. Use a fresh run ID for each experiment. See [figures/README.md](figures/README.md) for the figure filenames and `surplus` layout.

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
