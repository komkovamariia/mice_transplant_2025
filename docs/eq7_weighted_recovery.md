# Equation 7 recovery contrasts

The original Equation 7 score predicts a signed fraction of unique aaV
recovered from the training reference. This analysis also compares a ridge
model for the g1-minus-control contrast in unique-aaV recovery with a ridge
model for the contrast in abundance-weighted recovery.

For each aaV found in retained g1 and absent from retained g5 and g6,
the second model weights its query count by its detection frequency
across retained g1 mice. Counts are divided by all functional molecules
in the query mouse and expressed per million. Group means give each animal
equal weight. The response is the g1 mean minus the average of the g5 and
g6 means. The original structural shortlists are used for all comparisons.

The script reuses the 1249/1031/15000 count snapshots and three-mouse
omission evidence. Each test triplet stays out of the predictors, reference
and fit used to score it. The inner triplets produce training-only
outcomes. The independent sample size is 13 animals, even though each
chain has 75 test triplets. The two contrast fits use fixed ridge penalties
of 1 and 10. Results include top-k rankings, MAE against a matched
Constant baseline, and one-predictor deletions for the weighted fit.

## Run with an existing cache

Activate the repository's analysis environment and run from the repository
root. The cache can be the manuscript-aligned benchmark's _inputs directory
or a regular analysis cache with the same seeds and completed evidence.

```bash
python scripts/run_eq7_weighted.py \
  --cache results/eq7_benchmark_1249/_inputs \
  --output results/eq7_weighted_1249
```

To add external evaluation, set article_outputs to the absolute path of
analysis_outputs from the manuscript archive and use the completed TRUST4
run directories.

```bash
python scripts/run_eq7_weighted.py \
  --cache results/eq7_benchmark_1249/_inputs \
  --article-inputs "$article_outputs" \
  --skin results/eq7_external/PRJNA954716 \
  --marrow results/eq7_external/PRJNA650133 \
  --output results/eq7_weighted_1249_external
```

Each external directory must contain trust4_reports. The script checks the
report hashes and reproduces the original unique-aaV recovery measures
before calculating abundance-weighted recovery. PRJNA954716 supports the
TRB contrast for g1 versus the graft group. PRJNA650133 provides
graft-arm recovery for TRA and TRB.

The default runs on CPU. Add --device cuda to use a CUDA-enabled PyTorch
installation for the unrestricted ridge solver. Data processing and
the nonnegative solver use CPU.

## Figures

Use the completed benchmark tables and the weighted recovery results.

```bash
python scripts/plot_eq7_comparison.py \
  --benchmark results/eq7_benchmark_1249 \
  --weighted results/eq7_weighted_1249_external \
  --output figures/eq7_weighted_1249
```

To include the direct V-usage edgeR rankers in the selection figure,
pass --v-usage with the directory containing internal_v_usage_edger_summary.csv.
