# Fisher direction sensitivity and external aaV-fraction ridge validation

## Scientific questions

1. Does defining a positive pooled Fisher effect by the 2x2 **odds ratio**
   improve the model over the manuscript's raw UMI ratio?
2. Does the aaV-fraction ridge transfer to an independent skin-graft TRB
   cohort that has both intact and transplant arms?

This code **does not rebuild the structural rarefaction or the 75 fixed
held-out triplets**. It uses the existing complete
`eq7_benchmark_1249/_inputs` (structural seed 1249, calibration seed 1031,
depth 15,000) and the original aaV-fraction result directory, including
its per-fold edgeR effects. New per-fold Fisher tests are necessary because
the saved features only contain V-level fractions, not each aaV's
2x2 counts and p-value. These new Fisher results are resumably cached
one existing split at a time.

For each aaV c the change is:

```text
a = pooled functional UMI for c in g1
b = pooled functional UMI for c in g5 and g6
N1 = all functional UMI in g1
N0 = all functional UMI in g5 and g6

OLD direction: sign(log2((a+1)/(b+1)))
NEW direction: sign(log2((a+0.5)*(N0-b+0.5)
                         / ((b+0.5)*(N1-a+0.5))))
```

The Fisher exact p-value, BH FDR, aaV-level denominator, edgeR feature,
log residual count, shortlist, joint recovery target, ridge penalties
(1 and 10) and fixed 75 held-out triplets are unchanged. Both old and
new full models are trained and evaluated on the *same* rows.

The program asserts that (a) raw-count Fisher fractions recalculated on
each retained fold equal the previous cached fractions exactly, and
(b) the old full-model predictions and frozen coefficients match the
previous run to numerical precision. A mismatch aborts the run.

## External evidence without raw FASTQ or TRUST4 reruns

Provide the latest Overleaf ZIP that contains:

- `analysis_outputs/eq7_continuous_1249/external_predictions.csv`
- `analysis_outputs/PRJNA954716/validation/train_reference_scores.csv`
- `analysis_outputs/PRJNA650133/validation/{TRA,TRB}/train_reference_scores.csv`

The program uses the **previously evaluated observed joint recovery**
for the external mice from this archive. It calculates only the new
aaV-level statistical features on the existing matched **internal**
spleen raw-UMI panels. Frozen full-model coefficients fitted on the
75 internal triplets then predict the prespecified external candidates.
No external reads are reprocessed, no external rarefaction is redone.

Only PRJNA954716 (TRB) has a two-arm joint-recovery contrast and supports
a fair external MAE comparison with the baseline. PRJNA650133 (TRA/TRB)
contains graft recipients only: predictions are saved for description,
but its absent intact arm means it cannot supply an observed two-arm
contrast or MAE.

## Run in a separate worktree on Aldan-3

```bash
conda activate mice-transplant-2025
cd ~/projects/mice_transplant_2025/mice_transplant_git
git fetch origin analysis/eq7-fisher-or-external
git worktree add --detach ../eq7_fisher_or_external FETCH_HEAD
cd ../eq7_fisher_or_external
python -m pytest -q tests/test_eq7_fisher_or.py
```

Upload the latest Overleaf ZIP via JupyterLab into this new folder
and rename it to `article_latest.zip`. It must contain the three
`analysis_outputs/` paths described above. Then:

```bash
cache="$HOME/projects/mice_transplant_2025/mice_transplant_git/results/eq7_benchmark_1249/_inputs"
previous="$HOME/projects/mice_transplant_2025/eq7_aav_fraction/results/eq7_aav_fraction_1249"
article="$PWD/article_latest.zip"
test -f "$cache/v_input_seeds.json" && test -f "$previous/internal_predictions.csv" && test -f "$article"

mkdir -p logs/slurm
sbatch \
  --cpus-per-task=24 --mem=64G --time=1-00:00:00 \
  --export=ALL,EQ7_CACHE="$cache",EQ7_PREVIOUS="$previous",EQ7_OUTPUT="$PWD/results/eq7_fisher_or_external_1249",EQ7_ARTICLE_ZIP="$article",EQ7_PYTHON="$(command -v python)" \
  --output="$PWD/logs/slurm/%x-%j.out" \
  scripts/run_eq7_fisher_or_external.sbatch
```

For internal-only sensitivity, omit `EQ7_ARTICLE_ZIP` from `--export`;
the other analyses will run and save all internal outcomes.

The job requests 24 CPUs, but the R edgeR fits for the **two matched
external panels** need not parallelize across all CPUs. Per-fold edgeR
is never refitted. Fisher can use up to 24 workers.

## Outputs

- `fisher_cache/TRA/<hash>` and `TRB/<hash>`: per-fold odds fractions
- `internal_predictions.csv`: observed and predicted joint recovery
- `internal_error_summary.csv`: MAE versus Constant for old and new fits
- `internal_ranking_summary.csv`: top-k recovery for TRA and TRB
- `frozen_coefficients.csv`: coefficients trained using the identical 75 folds
- `external_predictions.csv` (when archive supplied): frozen predictions
  using stored external outcomes
- `external_error_summary.csv`: skin-graft TRB external MAE only;
  marrow two-arm MAE is undefined and therefore omitted
- `external_topk_summary.csv`: external observed TRB top-k contrasts
- `run.json`: provenance, original seeds and audit flags

Do not overwrite the previous aaV-fraction run or manuscript.
These are sensitivity analyses and the 75 splits share only 13 biological
mice. A lower internal MAE is not evidence of generalization until the
external comparison is assessed.
