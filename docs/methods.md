# Computational methods

## Inputs and analysis units

Approach 1 reads the study metadata, the clonoset index and the MiXCR exports referenced by that index. Clonotypes are defined as exact `(CDR3 amino-acid sequence, V segment)` pairs with zero mismatches. The historical TRA sample exclusions and the original chain-plus-sample-number metadata join are kept. There is no clone-level minimum UMI threshold.

The analysis is run in nine strata. `all_combined` contains CD4+ and CD8+ samples from thymus and spleen. The remaining strata split or pool cell type and tissue as listed in the main README.

Structural and count-based analyses use the same functional clonotypes but handle depth differently. Structural set analysis downsamples each eligible library independently to 15,000 functional UMI. Libraries below 15,000 UMI are omitted from this branch. Exact-set calculations use the union of the retained sample-level sets.

Count-based inference uses the original functional integer UMI counts. In pooled strata, repeated samples are summed within treatment group and mouse before Fisher and edgeR. These pooled values are depth-weighted sums of the available libraries.

## Representative structural seed

For `all_combined`, `cd4_combined` and `cd8_combined`, structural rarefaction is repeated for 100 baseline seeds.

Each draw uses its own adaptive cutoffs and retains every passing TRAV segment.
The representative draw minimizes the mean Jaccard distance between its
shortlist and the other valid draws. The smaller seed breaks exact ties.
The same seed is used for TRB figures. Only baseline data enter this selection.

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

The three predictors are edgeR and Fisher indicators for a positive mean
effect among FDR-significant aaV and the natural log of one plus the number
of unique residual aaV. The calibration omits every triplet consisting of
one g1, one g5 and one g6 mouse. With 3, 5 and 5 mice, this gives 75 training
splits with 2, 4 and 4 animals in those groups. All libraries of each omitted
mouse stay together. The other groups retain their original role in set
subtraction. The structural and count-based inputs and both dynamic cutoffs
are computed from the remaining training samples. Calibrating at a fixed
per-library rarefaction seed avoids representative-seed selection on these
folds.

A calibration row represents a split, an omitted mouse and a V segment
passing both training cutoffs. The response divides the number of that
mouse's aaV using the segment and present in the training g1-only set by the
number of its observed aaV using the segment. The training g1-only set is
defined relative to g5 and g6. The fraction is signed positive for g1 and
negative for g5 and g6. An unobserved V segment has a missing response.

The objective is weighted squared error plus a ridge penalty on the
slopes in standardized predictor space. Weight is apportioned equally to
groups, contributing mice within a group, eligible splits of a mouse and
usable V segments within that split. Constant predictors receive zero
slopes. Back-transformed slopes define the additive score `E` in the
original units. A separate unpenalized intercept is added for response
prediction. The final weights score every combined segment passing the
cutoffs of the displayed rarefaction profile.

The penalty is selected on 10 potential control-mouse validation folds.
For each held g5 or g6 mouse, the inner training folds omit one further
mouse at a time for calibration. Both mice are excluded before constructing
the structural and statistical features for that inner fit. When predicting
the held control mouse, predictors and candidate cutoffs use only its
complement. The five ridge penalties are evaluated by mean squared error,
averaged first by V within a mouse, then by mouse within group and equally
across the two control groups. The largest candidate within one estimated
standard error of the minimum is selected. This is a regularization heuristic
for dependent folds. The code reports an intercept-only reference. If a
comparison across the grid lacks two evaluable mice per control group,
the preset penalty is 0.1. A g1 tuning holdout with an additional g1 omission
leaves one mouse and violates edgeR's replication requirement. The g1
responses participate in the final triplet calibration.

### External marrow-transplant check in Zheng et al.

The primary marrow-transfer check uses the three day-9 BALB/c recipient
spleen libraries from Zheng et al. (2020), BioProject
[PRJNA650133](https://www.ebi.ac.uk/ena/browser/view/PRJNA650133). The donor
strain is C57BL/6. The published transplant included mature donor T cells
with the marrow. The libraries contain unsorted splenocytes. Both chains
are extracted with TRUST4 v1.1.7, retaining productive calls with a
chain-consistent J and an unambiguous V at gene resolution. External
structural depth is measured in read fragments, with the same nominal
depth used for the internal libraries measured in UMI.

The internal reference pools splenic CD4+ and CD8+ samples by mouse
for the g1 versus g5 count contrast. For the structural branch, the
reference is the union of eligible rarefied splenic libraries in each
group. Candidate cutoffs and predictions use this internal reference,
and the combined-run TRA and TRB coefficients stay frozen. For each
eligible external mouse, the response is the negative fraction of its
unique aaV for a candidate V segment observed in the internal g1-only
set. The output reports coverage and mean absolute error averaged by
mouse within the g5-like arm, with an intercept-only reference. An
intact BALB/c spleen arm is required to reproduce the two-group count
contrast and compare independently selected candidate lists. The
Zheng cohort supports the single-arm recovery endpoint.

### Two-arm external TRB transfer

The external panel uses Tereshchenko et al. (2023), BioProject
[PRJNA954716](https://www.ebi.ac.uk/ena/browser/view/PRJNA954716).
Three intact BALB/c and four BALB/c mice given C57BL/6 skin grafts
provide the two conditions. Both are bulk libraries from splenic CD4+
cells. Paired reads are assembled independently with TRUST4 v1.1.7 and
the bundled mouse IMGT reference in bulk non-UMI mode. Productive TRB aaV
are collapsed at V-gene resolution and use integer read-fragment counts.

The transfer reference restricts the internal data to splenic CD4+
libraries of g1 and g5. Each cohort uses an intact versus allograft count
contrast and subtracts its corresponding graft group only. Structural
rarefaction uses the original internal depth in UMI and the same nominal
depth in external read fragments. Both cohorts use the fitted combined TRB
coefficients and the cutoffs obtained from the internal reference.

The mouse-level endpoint compares frozen internal predictions with signed
aaV recovery from the internal intact-only set in each external mouse.
The mean absolute error is averaged across segments within each mouse,
then within each group and equally across groups, with an intercept-only
reference. The cohort-level endpoint compares the selected lists and
scores after re-estimating edgeR, Fisher and residual counts on the external
cohort. A segment passing exactly one list contributes one unit to loss.
A shared segment contributes the absolute score difference divided by its
larger score, or zero if both scores are zero. The loss is averaged across
the union. The design also reports Jaccard overlap, absolute shifts for
shared segments and group-level coverage. See [RUN_GUIDE.md](../RUN_GUIDE.md)
for commands and outputs.

## Stability checks

The representative-seed strata include two structural checks.

All V segments passing both adaptive cutoffs in the displayed profile define
the reference for that chain. Each of 100 library rarefaction reruns refits
both cutoffs. We report shortlist sizes, intersection counts, reference
recovery and Jaccard similarity. The per-V boxplots report the UMI share in
the g1-only overlap region relative to g5 and g6 for every reference segment
and distinguish a missing V segment from a measured zero share.

A second analysis balances g1, g2, g4, g5 and g6 to the smallest group size.
It retains all eligible libraries of each sampled mouse and refits both
cutoffs in each of 50 draws. The CSVs identify unavailable fits and
comparisons. The perturbations in both checks use baseline data.

## Validation

The test suite covers sample selection, metadata joins, mouse pooling, structural downsampling, `S(v)`, four-region entropy, output persistence and spike-in propagation. CI also runs the notebook edgeR cell on synthetic count data with R and rpy2.

Study-level values require the study data and MiXCR exports. The repository tests check the implementation and output contract, not the biological effect sizes.

The separate perturbation experiment is described in [spike_in.md](spike_in.md).
