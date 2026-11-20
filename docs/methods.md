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

Standard runs fit the two graphical thresholds separately for each chain
and panel using the displayed rarefied libraries:

```text
f_full(v) >= median(f_full)
r(v) >= Q_0.75(r | f_full >= median(f_full))
```

The quantiles use linear interpolation. Four segments with valid positive
initial shares are required. The retention cutoff includes ties and is
positive whenever a segment on the right was retained. If the upper
quartile is zero, the smallest positive retention on the right is used.
The thresholds describe a relative shortlist. They do not identify a
biological boundary or control a false discovery rate. Leave-one-mouse-out
cutoffs are reported to show sensitivity to the included mice.

The dumbbell plot orders V segments by decreasing `S(v)` after applying
both cutoffs. The combined bubble annotates up to 11 segments and the
separate CD8+ and CD4+ panels annotate up to 10 each. Each heatmap uses
its paired bubble's annotated segments and orders them by g1-only
regional share, then by `S(v)`. The subset panels restrict the combined
run's rarefied libraries to the indicated subset and recompute the set
intersections. The representative-seed selection and spike-in family
selection use their historical eligibility rules.

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

## Depth check and MDS

The depth check uses the same g1/g2 libraries at all five depths (10,000, 13,000, 15,000, 20,000 and 30,000 functional UMI). At each depth, 50 independent draws are compared with the full-depth profiles of those libraries. The figure shows the median and 5th–95th percentiles of Spearman correlation between pairwise F2 similarities and Procrustes disparity between the two-dimensional classical MDS configurations. The manuscript MDS scatter uses the main rarefaction depth and its eligible g1/g2 libraries. Classical MDS uses `-log10(F2)` as the dissimilarity. For zero overlap, F2 is replaced by half the lowest strictly positive off-diagonal similarity. The horizontal break is drawn when a large internal gap is present.

## Differential abundance

The combined analysis fits Fisher and edgeR separately for TRA and TRB.

Fisher testing compares pooled g1 with pooled g5+g6 counts. The displayed effect size is:

```text
log2((g1_count + 1) / (g5_count + g6_count + 1))
```

Benjamini-Hochberg correction is applied across the tested aaV clonotypes.

edgeR uses mouse-level analysis units with group levels `allogeneic` and `g1`. Positive logFC indicates higher abundance in g1. Negative logFC indicates higher abundance in the pooled g5+g6 reference. The model sequence is `DGEList`, `calcNormFactors`, `filterByExpr`, `estimateDisp`, `glmQLFit` and `glmQLFTest`. At least two analysis units are required in each contrast group. The figure shows up to eight V segments per direction of the mean significant-aaV effect. The support panel counts significant aaV in both directions.

The baseline V-ranking tables report the historic Equation 7 score:

```text
4 * edgeR_g1_enriched
+ 2 * fisher_g1_enriched
+ log1p(n_g1_exclusive_aav)
```

Ties use significant-aaV counts, exclusive-aaV counts and V name.

### Calibrated Equation 7 in the combined stratum

The calibration fits three nonnegative coefficients. Its predictors are
the edgeR and Fisher indicators for a positive mean effect among
FDR-significant aaV and the natural log of one plus the number of unique
residual aaV. A row represents a V segment and an omitted mouse. Its
response is the fraction of that mouse's unique aaV with the segment
that occur in the training `g1 only` set. The fraction is signed
positive for g1 and negative for g5 and g6. Rows without observations
for the segment in the omitted mouse are excluded.

The objective is weighted squared error plus a ridge penalty on the
slopes in standardized predictor space. Each group has equal total
weight, each contributing mouse has equal weight within its group and
each usable row has equal weight within its mouse. Constant predictors
receive zero coefficients. Back-transformed slopes define the additive
ranking score `E`. The fitted response adds a separate, unpenalized
intercept. The response and metric describe overlaps in this cohort.

In outer validation, both the outer and inner mouse are excluded
from the structural and count calculations for every calibration row.
The cutoffs and predictor standardization use the respective training
data. The outer scores then use predictors computed on all mice except
the outer test mouse. With three g1 mice, holding out two leaves one
g1 replicate and prevents calibration of outer g1 folds. The
full-cohort fit uses one-mouse omissions. Its learned coefficients are
applied to every candidate passing the cutoffs in the same rarefaction
profile used for the displayed combined figures. See
[RUN_GUIDE.md](../RUN_GUIDE.md) for the three diagnostic plots.

## Stability checks

The representative-seed strata include two structural checks.

The 100 rarefaction reruns report how many reference top-10 V segments are recovered in each rerun. TRA also has a boxplot of the g1-only regional share for each consensus top-10 TRAV segment across the 100 runs. The box shows the interquartile range and median; the mean and representative-seed value are overlaid.

A second analysis repeats structural subtraction after balancing the number of mice contributed by g1, g2, g4, g5 and g6 to the smallest group size. This is repeated for 50 seeds. It checks whether unequal mouse counts are driving the structural shortlist.

Neither stability analysis uses spike-in data.

## Validation

The test suite covers sample selection, metadata joins, mouse pooling, structural downsampling, `S(v)`, four-region entropy, output persistence and spike-in propagation. CI also runs the notebook edgeR cell on synthetic count data with R and rpy2.

Study-level values require the study data and MiXCR exports. The repository tests check the implementation and output contract, not the biological effect sizes.

The separate perturbation experiment is described in [spike_in.md](spike_in.md).
