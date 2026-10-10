# Computational methods

## Inputs and analysis units

Approach 1 reads the study metadata, the clonoset index and the MiXCR exports referenced by that index. Clonotypes are defined as exact `(CDR3 amino-acid sequence, V segment)` pairs with zero mismatches. The historical TRA sample exclusions and the original chain-plus-sample-number metadata join are kept. There is no clone-level minimum UMI threshold.

The analysis is run in nine strata. `all_combined` contains CD4+ and CD8+ samples from thymus and spleen. The remaining strata split or pool cell type and tissue as listed in the main README.

Structural and count-based analyses use the same functional clonotypes but handle depth differently. Structural set analysis downsamples each eligible library independently to 15,000 functional UMI. Libraries below 15,000 UMI are omitted from this branch. Exact-set calculations use the union of the retained sample-level sets.

Count-based inference uses the original functional integer UMI counts. In pooled strata, repeated samples are summed within treatment group and mouse before Fisher and edgeR. These pooled values are depth-weighted sums of the available libraries.

## Representative structural seed

For `all_combined`, `cd4_combined` and `cd8_combined`, structural rarefaction is repeated for 100 baseline seeds.

For each rerun, V segments are ranked by the preservation score `S(v)` defined below. A consensus top-10 TRAV list is built from how often a segment appears in the rerun top 10, with median `S(v)` used for tie-breaking. Each rerun is then represented by the g1-only regional-share values of those consensus TRAV segments.

The representative seed is the rerun with the smallest mean absolute distance from the across-rerun median profile. Ties are resolved by larger consensus top-10 overlap and then by lower seed. Only baseline data are used for this selection.

## V-segment preservation

For V segment `v`:

```text
f_full(v) = n_full(v) / N_full
r(v) = n_rem(v) / n_full(v)
S(v) = f_full(v) * r(v)
```

`n_full(v)` is the number of unique g1 aaV clonotypes using `v`, `n_rem(v)` is the number left after subtracting the union of g2, g4, g5 and g6, and `N_full` is the total number of unique g1 aaV clonotypes.

The graphical thresholds are:

```text
f_full(v) >= 0.006
r(v) >= 0.66
```

Eligible V segments are ordered by decreasing `S(v)`. The dumbbell plot uses this order. The bubble plot labels the seven highest-scoring displayed V segments. The localization heatmap uses up to ten threshold-eligible V segments and orders its rows by g1-only regional share, then by `S(v)`.

The heatmap uses four mutually exclusive g1-containing regions:

```text
g1 only
(g1 ∩ g5) without g6
(g1 ∩ g6) without g5
g1 ∩ g5 ∩ g6
```

For each V segment, regional values are fractions of its g1 UMI burden. Normalized regional entropy is:

```text
H(v) = -sum(p_region(v) * log2(p_region(v))) / log2(4)
```

Zero-probability terms contribute zero. This UMI-weighted entropy is separate from the unique-aaV preservation score `S(v)`.

## Differential abundance

Fisher testing compares pooled g1 with pooled g5+g6 counts. The displayed effect size is:

```text
log2((g1_count + 1) / (g5_count + g6_count + 1))
```

Benjamini-Hochberg correction is applied across the tested aaV clonotypes.

edgeR uses mouse-level analysis units with group levels `allogeneic` and `g1`. Positive logFC means enrichment in g1. The model sequence is `DGEList`, `calcNormFactors`, `filterByExpr`, `estimateDisp`, `glmQLFit` and `glmQLFTest`. At least two analysis units are required in each contrast group.

The final TRAV table keeps differential evidence separate from the structural score. Its evidence score is:

```text
4 * edgeR_g1_enriched
+ 2 * fisher_g1_enriched
+ log1p(n_g1_exclusive_aav)
```

Ties use significant-aaV counts, exclusive-aaV counts and V name. This ranking is not the same quantity as `S(v)`.

## Stability checks

The representative-seed strata include two structural checks.

The 100 rarefaction reruns report how many reference top-10 V segments are recovered in each rerun. TRA also has a boxplot of the g1-only regional share for each consensus top-10 TRAV segment across the 100 runs. The box shows the interquartile range and median; the mean and representative-seed value are overlaid.

A second analysis repeats structural subtraction after balancing the number of mice contributed by g1, g2, g4, g5 and g6 to the smallest group size. This is repeated for 50 seeds. It checks whether unequal mouse counts are driving the structural shortlist.

Neither stability analysis uses spike-in data.

## Validation

The test suite covers sample selection, metadata joins, mouse pooling, structural downsampling, `S(v)`, four-region entropy, output persistence and spike-in propagation. CI also runs the notebook edgeR cell on synthetic count data with R and rpy2.

Study-level values require the study data and MiXCR exports. The repository tests check the implementation and output contract, not the biological effect sizes.

The separate perturbation experiment is described in [spike_in.md](spike_in.md).
