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

Standard runs use dynamic V-segment cutoffs and write them to
`results/01_set_count/<stratum>/v_thresholds.csv`. Use
`--v-thresholds fixed` to reproduce the historic 0.6% and 70% graphical cutoffs.

### Redraw article figures from a completed run

The adaptive figure directory initially has five per-chain figures after
`thresholds`. The ridge diagnostic figures are written during `weights`.
The two recovery histograms and per-segment boxplot are generated in
`stability` and stored in the chain's `surplus` directory.

```bash
python scripts/run_v_updates.py redraw --components structural stability
python scripts/run_v_updates.py redraw --components eq7
```

The redraws use saved count snapshots, structural and stability tables,
and the completed Equation 7 CSVs. The files are written as print-resolution PNG to `figures/01_set_count_mdpi/all_combined/`. The redraw stage leaves
analysis tables untouched. Omit `--components` for all three families.
Use `--figure-output PATH` for a different destination. The new text
sizes are calibrated to a manuscript figure width of approximately
180 mm.

If adaptive stability tables have never been computed, run:

```bash
python scripts/run_v_updates.py stability --figure-output figures/01_set_count_mdpi/all_combined
```

This computes structural draws from the saved snapshots and updates the
stability CSV files. It does not rerun differential-abundance tests or
regression. If Equation 7 is still running, redraw its figures after
`calibration_summary.csv` becomes available in the `eq7_weights` directory.

### Refresh thresholds and calibrate Equation 7

With an existing combined baseline, run:

```bash
python scripts/run_v_updates.py all
```

If the baseline is missing, first run
`python scripts/run_analysis.py --approach 1 --strata all_combined`.
Use `--baseline-results PATH` to point to a baseline directory with
`run_summary.json` and the TRA and TRB V-ranking CSVs. To rebuild the snapshots
when the underlying sample files change, use `python scripts/run_v_updates.py prepare`
with that baseline and any input path overrides.

The `all` command refreshes the V figures, caches the mouse-omission
evidence and fits Equation 7 for TRA and TRB. To run the steps
separately:

```bash
python scripts/run_v_updates.py thresholds
python scripts/run_v_updates.py stability
python scripts/run_v_updates.py evidence
python scripts/run_v_updates.py weights
```

`thresholds` refreshes the structural figures and the signed edgeR
summary from the saved counts and results. The equivalent
`run_analysis.py --approach 1 --strata all_combined --v-thresholds adaptive`
command also refreshes only these figures. They are written to
`figures/01_set_count_adaptive/`. The matching CSV in
`results/01_set_count_adaptive/all_combined/v_thresholds.csv` includes the
full-profile thresholds and the sensitivity of those thresholds to the
omission of each mouse. Both cutoffs are recomputed for each chain and
panel. Retention is the upper quartile among segments above the median
initial share.

`stability` uses the cached rarefied tables for the displayed shortlist and
the cached raw tables for the 100 library draws. It refits both adaptive
cutoffs on each draw and in 50 balanced-mouse draws, with no top-N cap.
It writes the CSVs under `results/01_set_count_adaptive/all_combined/` and
the PNGs under each chain’s `surplus` directory. It does not fit edgeR or
Equation 7. An existing baseline keeps its displayed seed in this stage.
A new full standard run chooses a seed using the adaptive TRAV shortlists.

`evidence` computes every triplet of omitted mice with one from each of
g1, g5 and g6. There are 75 triplets with the original mouse counts.
It also computes single omissions for internal holdout features and paired
omissions for the inner training data used during penalty selection. These
are separate split families. In a paired split, one omitted mouse is the
outer holdout and the other supplies an inner recovery response. Pairs
with two g1 omissions can be recorded as unavailable because one g1 mouse
remains and edgeR requires two per contrast arm. They supply no calibration
rows. The triplet calibration retains two g1 mice.
All libraries from each omitted mouse are removed together. The structural
counts use a fixed, per-library rarefaction. The step fits edgeR and Fisher
on raw training counts and caches the evidence for each split.

`weights` evaluates nonnegative ridge penalties of 0.001, 0.01, 0.1, 1
and 10. It uses the largest penalty within one estimated standard error
of the minimum group-balanced g5/g6 mouse loss. If fewer than two g5 or
g6 mice provide comparable losses, the recorded preset is 0.1. Use
`--ridge` to supply an explicit penalty. An optional `--top N` limits
calibration rows in each split by decreasing g1-only UMI share.
The final fit uses all usable triplets, and all eligible combined V
segments receive an E score. With three g1 mice, edgeR cannot fit nested
g1 holdouts that leave only one g1 replicate.

| Output | What it shows |
|---|---|
| `*_coefficients.csv` | The three learned slopes, their weighted standardization and the prediction intercept |
| `*_calibration_rows.csv` | The training predictors, their standardized values and the observed and fitted signed responses |
| `*_stratified_folds.csv` | The three omitted mice, the cutoffs and the number of usable observations per calibration split |
| `*_scored_top_v.csv` | Every combined candidate passing both displayed cutoffs, its E and the three component contributions |
| `*_ridge_grid.csv` and `*_ridge_mouse_losses.csv` | Penalty selection, losses for each mouse and the intercept-only baseline |
| `*_outer_scores.csv` and `*_outer_status.csv` | Internal predictions used during penalty selection and their coverage |
| `calibration_summary.csv` | Fit and scoring status for each chain |

These tables are under `results/01_v_update_inputs/all_combined/eq7_weights/`.
The archive `figures/01_set_count_adaptive.zip` includes four Equation 7
plots per chain. They show the standardized design matrix, calibration
scatter and internal predictions, penalty curve and V scores by contribution.
The fitted and predicted response adds an intercept to E.

### External marrow-transplant check in Zheng et al.

Build [TRUST4 v1.1.7](https://github.com/liulab-dfci/TRUST4/releases/tag/v1.1.7).
The default external study is [PRJNA650133](https://www.ebi.ac.uk/ena/browser/view/PRJNA650133).
The sample manifest selects three BALB/c spleens on day 9 after C57BL/6
marrow transplantation and pins both FASTQ checksums for each run.

```bash
python scripts/run_eq7_external.py run --study PRJNA650133 --trust4 ~/tools/TRUST4-1.1.7/run-trust4 \
  --expect-run 1249 1031 15000
```

The test uses the frozen TRA and TRB weights with internal splenic CD4+
and CD8+ libraries pooled by mouse. It uses the internal g1 versus g5
reference to compute cutoffs, candidates and predictions. From Zheng's
three spleens it computes observed signed recovery and mouse-level mean
absolute error for the g5-like arm, with the intercept-only baseline.
The archive has no corresponding intact BALB/c spleen group. This check
uses a single arm, so it exports no external group contrast, list Jaccard
index or E union loss. Read-fragment and UMI counts retain their own units.
See `results/eq7_external/PRJNA650133/validation/` for the summary and
chain-specific tables.

### Two-arm TRB validation in Tereshchenko et al.

[PRJNA954716](https://www.ebi.ac.uk/ena/browser/view/PRJNA954716)
includes three intact BALB/c mice and four
BALB/c recipients of C57BL/6 skin grafts, profiled in splenic CD4+ T cells.

```bash
python scripts/run_eq7_external.py run --study PRJNA954716 --trust4 ~/tools/TRUST4-1.1.7/run-trust4 \
  --expect-run 1249 1031 15000
```

The `prepare` stage downloads the seven paired FASTQ sets, checks their
ENA checksums and assembles them with the supplied TRUST4 mouse IMGT
reference. The `test` stage uses the fitted TRB coefficients and numeric
cutoffs from the internal g1 versus g5 reference restricted to splenic CD4+
libraries. It fits the count tests for that reference and the external group
contrast because their samples differ from the combined baseline. The
combined results are reused. The penalty and both cutoffs are frozen from
the internal data before computing the external results.

Under `results/eq7_external/PRJNA954716/validation/` the outputs include
the frozen coefficients, reference cutoffs, library depth, per-mouse recovery
errors and the comparison of both selected V lists. The union loss charges
one unit for a V segment passing in just one cohort. For shared segments
it reports the absolute E shift divided by the larger E, with zero when
both scores are zero. The Jaccard index and absolute E differences are
reported separately. The two PNGs show candidate coverage and absolute error by mouse, and paired
structural coordinates for the union of selected V segments. This study profiles
TRB. Its skin-graft condition tests transfer across interventions.

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
trb_rarefaction_v_profiles.csv
{tra,trb}_adaptive_only_g1_by_seed.csv
{tra,trb}_adaptive_only_g1_summary.csv
{tra,trb}_rarefaction_adaptive_stability.csv
{tra,trb}_balanced_mouse_adaptive_stability.csv
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

For the representative-seed strata, the stability outputs also include `TRA_g1_rarefaction_stability.png`, `TRB_g1_rarefaction_stability.png`, `TRA_g1_balanced_mouse_stability.png`, `TRB_g1_balanced_mouse_stability.png` and `{TRA,TRB}_g1_adaptive_only_g1_stability.png`.

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

### Repertoire QC and depth figures

The exploratory `presenting_repseq.ipynb` includes clonoset summaries and the
early F2/MDS depth checks. The paper's depth checks use
[`src/repertoire_depth.py`](src/repertoire_depth.py) from the standard analysis.
The redraw stage reads the saved depth and MDS tables. Library summaries and
exploratory CDR3 and diversity analyses do not need to run to redraw the paper.

For the manuscript anchor, verify the cached seeds and depth before drawing:

```bash
python scripts/run_v_updates.py redraw --expect-run 1249 1031 15000
```

The values refer to the structural snapshot, calibration rarefaction and
target depth, in that order. Redraw leaves the existing cache and coefficients
in place and writes print-resolution PNGs. Older PDFs are not overwritten.

To redraw an external evaluation from its saved validation tables:

```bash
python scripts/run_eq7_external.py plot --study PRJNA954716
python scripts/run_eq7_external.py plot --study PRJNA650133
```

Use `--work-dir` if the existing results are outside the default
`results/eq7_external/<BioProject>/` directory. The plot step checks the
coefficient checksum and does not run TRUST4 or refit the model.
