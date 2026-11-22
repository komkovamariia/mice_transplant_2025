# Benchmark: fraction of positively enriched aaV within each V segment

This experiment compares a new aaV-level evidence encoding with the existing
binary aaV indicators and the V-level continuous-effect model.

**No mouse partitions, rarefaction or structural cutoffs are recalculated.**
The script reads the existing 75 triplets per chain and 1249/1031/15000
snapshots under `results/eq7_benchmark_1249/_inputs`. Only the aaV-level
edgeR and Fisher tests that are necessary for the new continuous predictors
are computed. Every test excludes all libraries from the three held-out mice.
Per-triplet results are written atomically into an independent `effect_cache`,
so interrupted jobs can resume without repeating successful folds.

## Definition

For each V segment and each statistical test, calculate:

\[
x_j(v)=\frac{
 \#\{c:V(c)=v,\ q_{j,c}<0.05,\ \log_2FC_{j,c}>0\}
}{
 \#\{c:V(c)=v,\ c\text{ was tested by }j\}
}, \quad j\in\{\mathrm{edgeR},\mathrm{Fisher}\}.
\]

The denominator includes all tested aaV for this V, including those with
non-significant or negative effects. It differs between tests because edgeR
filters low-count aaV before testing, whereas Fisher tests all aaV with
positive pooled counts in the contrasted groups. If no aaV was tested,
the corresponding feature is set to zero. Fisher uses the **original
raw-pooled `Fisher_log2FC_g1`** definition to remain comparable to the
manuscript's binary Fisher indicator; unequal group sizes affect that
definition and require caution.

The model uses `(edgeR_positive_fraction, Fisher_positive_fraction,
log1p_n_g1_exclusive_aav)` and the **unchanged joint recovery contrast**
from the continuous-evidence experiment. A separate ridge model is
retrained after removing each feature. Penalties are fixed at 1 and 10
for the full model; ablations use 1. All preprocessing for the new
features excludes held-out mice.

## Run on Aldan-3 without the 6-hour preparation

Use the actual existing results paths. From a clean checkout or Git worktree
containing this branch:

```bash
conda activate mice-transplant-2025
bench="$HOME/projects/mice_transplant_2025/mice_transplant_git/results/eq7_benchmark_1249"
previous="$HOME/projects/mice_transplant_2025/mice_transplant_git/results/eq7_continuous_1249"
python -m pytest tests/test_eq7_aav_fraction.py -q
python scripts/run_eq7_aav_fraction.py \
  --cache "$bench/_inputs" \
  --reference-results "$previous" \
  --output results/eq7_aav_fraction_1249 \
  --workers 4
```

If `$previous/internal_predictions.csv` does not exist yet, omit the
`--reference-results` argument. The new model still runs, but the
combined comparison with the old binary and V-continuous models will
then be unavailable.

For a Slurm batch, request one node with sufficient RAM and CPUs, and run
the command above through `sbatch`. **Do not rerun** `prepare_eq7_1249.py`
or `build_weight_evidence`.

## Outputs

All outputs are under the new directory and leave the original benchmark
and manuscript unchanged:

* `effect_cache/{TRA,TRB}/<hash>/*.csv`: held-out, resumable aaV test fractions
* `internal_predictions.csv`: per-V observed and predicted joint recovery
* `internal_error_summary.csv`: MAE versus matched Constant baseline
* `internal_ranking_by_mouse_split.csv`: top-k mouse-level observations
* `internal_ranking_summary.csv`: recovered top-5, top-10 TRA and top-1,
  top-2 TRB joint contrast
* `frozen_coefficients.csv`: coefficients for descriptive full-data fits
* `run.json`: provenance and unchanged seed settings

Comparison with previous binary and V-continuous predictions requires
matching mouse splits, V candidates and held-out target values. The program
checks this and fails if observations differ.

A reduction in MAE relative to Constant and the existing models is
required before claiming additional predictive value. The 75 triplets
reuse 13 biological mice, so their mean errors are descriptive rather
than 75 statistically independent estimates. This experiment does not
establish TCR antigen specificity or GVHD association.
