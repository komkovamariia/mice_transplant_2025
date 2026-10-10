# Computational methods

## Analysis units and clonotype definition

Approach 1 uses the study metadata, a clonoset-index CSV and referenced MiXCR exports.
Exact clonotypes are defined by `(CDR3 amino-acid sequence, V segment)` with zero
mismatches. No minimum threshold of three UMI is applied. Historical sample exclusions
and metadata-key handling are documented in the [equivalence audit](approach1_historical_equivalence.md).

The first stratum, `all_combined`, selects the union of CD4 and CD8 samples from thymus
and spleen, including source names mapped to thymus by the existing metadata rules.
It restores the original combined sample scope. The other eight strata are retained.
For the structural branch, MiXCR clonotype exports are first restricted to functional
rearrangements and interpreted using UMI counts. Each eligible library is independently
downsampled to 15,000 functional UMI before the exact aaV sets are constructed. Libraries
below this depth are excluded from the structural calculation. For all_combined, cd4_combined and cd8_combined, 100 independent baseline rarefactions
are generated. The ten TRAV segments most consistently present among the top-ten S(v)
candidates define a ten-dimensional g1-only regional-share profile. For each rerun, the
mean absolute deviation of this profile from the across-rerun median profile is
calculated. The rerun with the smallest deviation is selected as the representative
seed; ties are resolved by higher overlap with the consensus top-ten TRAV list and then
by lower seed. The main structural figures are generated from this representative rerun.
Spike-in data are not used in seed selection. Exact-set calculations then use the union of the retained
sample-level sets. No clone-level minimum UMI threshold is imposed.

For count-based inference, the same functional clonotypes are retained at their original
integer UMI counts rather than rarefied. Pooled strata sum these counts within treatment
group and mouse before edgeR and Fisher analysis. They represent depth-weighted pools of
the available libraries, not equally weighted tissue averages. Strata overlap and should
not be treated as independent biological replicates of one another.

## V-segment preservation and figures

For V segment `v`, let `n_full(v)` be the number of unique aaV in g1 and `n_rem(v)`
the number remaining after subtraction of the union of g2, g4, g5 and g6. Let
`N_full` be the total number of unique g1 aaV.

```text
f_full(v) = n_full(v) / N_full
r(v) = n_rem(v) / n_full(v)
S(v) = f_full(v) * r(v)
```

Graphical candidates require `f_full(v) >= 0.006` and `r(v) >= 0.66`. Up to ten are
selected in descending preservation score `S(v)`, with deterministic tie handling.
The dumbbell plot is ordered explicitly by decreasing `S(v)`. The localization heatmap uses the complete threshold-eligible top-ten shortlist, while
the bubble plot labels the seven highest S(v) values among all displayed V segments to
avoid label crowding. Heatmap rows are ordered by decreasing
g1-only regional share and then by `S(v)`. The bubble plot's lower Y limit adapts to the
displayed values.

Only the four mutually exclusive g1-containing regions contribute to the heatmap:
`g1 only`, `(g1 intersect g5) without g6`, `(g1 intersect g6) without g5`, and the
triple intersection. Within each V segment, regional probabilities are the fractions
of its g1 UMI burden in these regions. Normalized entropy is

```text
H(v) = -sum(p_region(v) * log2(p_region(v))) / log2(4)
```

Zero-probability contributions are zero. UMI-weighted regional entropy is distinct
from the unique-aaV preservation score. There is one TRA heatmap and one TRB heatmap,
with no cell grid and annotation colors chosen for contrast. Figure titles identify
chain, treatment context, cell subset and organ. Auxiliary pairwise and g2-centered
figures remain outside the article output; relevant non-g2 numerical tables are retained.

## Differential abundance

Fisher tests compare pooled g1 versus pooled g5+g6 UMI counts. The historical displayed
effect is `log2((g1_count + 1) / (g5_count + g6_count + 1))`. Benjamini-Hochberg correction
is applied across tested aaV. Fisher results are a sensitivity analysis of pooled counts.

The edgeR model uses independent analysis units with group levels `allogeneic` and
`g1`; positive logFC denotes enrichment in g1. The retained sequence is `DGEList`,
`calcNormFactors`, `filterByExpr`, `estimateDisp`, `glmQLFit` and `glmQLFTest`.
At least two analysis units per contrast group are required. The workflow requires a
successful fitted model and records its R session information.

Final TRAV tables distinguish edgeR/Fisher concordance, single-method support and
exact-set descriptive evidence. The existing final evidence score is
`4 * edgeR_g1_enriched + 2 * fisher_g1_enriched + log1p(n_g1_exclusive_aav)`;
ties use significant-aaV counts, exclusive-aaV counts and V name. This evidence score
has a different purpose from `S(v)`, which selects bubble-plot and heatmap genes.

## Robustness checks

The all_combined, cd4_combined and cd8_combined strata report structural robustness
checks in addition to the ordinary analysis. First, the 15,000-UMI rarefaction is
repeated for 100 independent seeds. Top-ten overlap is measured against a consensus
top-ten list derived from those baseline reruns and displayed as the frequency
distribution of the number of consensus segments recovered. A separate TRA boxplot
summarizes the g1-only regional share of the consensus top-ten TRAV segments across all
100 reruns; boxes show the interquartile range and median, while the mean and the
representative-seed value are overlaid explicitly. Second, structural subtraction is
repeated 50 times after balancing the number of mice
contributed by g1, g2, g4, g5 and g6 to the smallest group size. This tests whether the
shortlist is driven by unequal opportunities to observe rare clonotypes in groups with
different numbers of animals. Both stability analyses use baseline data only and are not
computed from spike-in perturbations.

## Computational validation

Regression tests check sample selection, mouse pooling, metadata joins, article score,
four-region entropy, output persistence, functional 15,000-UMI structural processing and
spike-in propagation into the original notebook calculations. A separate CI job executes
the actual notebook edgeR cell on synthetic count data with R and rpy2. The source
notebooks are committed without outputs.

Study-level numerical results must be reproduced from the supplied data. Passing CI
establishes the tested software behavior; it does not establish biological effect sizes
or numerical equivalence with the historical unpooled model. See the
[spike-in protocol](spike_in.md) for the separate perturbation experiment.
