# Approach 1 historical-equivalence audit

The active [`venn_original.ipynb`](../venn_original.ipynb) was reconstructed from the
last verified notebook-first implementation at commit
`6678766a0e8913ec2feacb3efb2daebd6a27eda4`, blob
`043ec6bf8227b72ea3e3bfc6683aaac1622fe341`.

## Preserved operations

| Stage | Preserved notebook behavior |
|---|---|
| Inputs | Read `metadata_mice_transplant.csv`, `clonosets_mice_transplant_2025_df.csv`, and its referenced MiXCR exports. |
| Historical key | Construct `<metadata chain>-<sample_no>` before chain normalization, preserving keys such as `alpha-100` and `beta-100`. |
| Exclusions | Retain the three prespecified TRA sample exclusions from the historical notebook. |
| Clonotype | Use exact `(CDR3aa, V)` features through `overlap_type="aaV"` and `mismatches=0`. |
| TRA/TRB sets | Retain g1 to g6 count tables, the g1-centered five-circle set diagrams, exact g1 subtraction, V-segment retention tables, and TRA/TRB bubble plots. g2-centered duplicates are excluded from this study. |
| Repertoire similarity | Retain non-g2 mouse-level TRA/TRB Jensen-Shannon, rank-correlation, and UMI-correlation tables. Auxiliary pairwise figures are excluded from the article output. |
| Fisher | Compare pooled g1 and g5+g6 feature counts with the historical pseudocount-adjusted log2 ratio and Benjamini-Hochberg correction. |
| edgeR | Retain TMM normalization, `filterByExpr`, quasi-likelihood fitting, and the g1 coefficient relative to pooled g5+g6 annotation. |
| Thresholds | Introduce no minimum UMI or read-count threshold. |

The active differential section contains edgeR and the Fisher sensitivity analysis
only. No additional differential-expression model is fitted.

The runner treats edgeR fitting as a required output contract. It reports a failed
stratum when the quasi-likelihood model does not complete, which prevents an earlier
plotting or dependency error from being interpreted as a complete differential result.

## Eight-stratum extension

The historical notebook analyzed all available compartments together. The active
notebook applies the same operations independently to eight prespecified strata. Sample
selection occurs before the first `repseq.intersections.count_table` call. In the four
combined strata, multiple selected samples from the same mouse are summed before the
edgeR model; exact-set operations still use the union of selected samples.

Each execution receives one `MICE_TCR_STRATUM` value and creates one executed notebook,
log, figure set, table set, and conclusion. The source notebook contains no hard-coded
TRAV conclusion from the unstratified historical run.

## Output-completeness extension

The historical notebook saved some figures explicitly and displayed others only
inline. The restored notebook registers a Matplotlib figure hook before analytical
cells execute. Explicit and displayed figures are saved as PNG and PDF under
`figures/01_set_count/<stratum>/`. A per-stratum manifest is checked by the runner.

For each chain, the bubble plot labels the ten V segments with the highest preservation
score from the article, `S(v) = f_full(v) × r(v)`. The associated heatmap is restricted
to that same ordered set after the graphical thresholds of 0.6% initial share and 66%
retention, and to the four regions that contain g1. Entropy is normalized
by `log2(4)`. Heatmap annotations switch between light and dark text using rendered cell
luminance. No pairwise similarity, UMI-correlation, cross-stratum, or mouse-level Venn
heatmaps are generated.

After all eight executions, the final run consolidates the stratum-level TRAV tables
and creates one `figures/01_set_count.zip` archive. Descriptive exact-set evidence is
kept separate from FDR-supported inference.

## Verification boundary

Repository tests verify the sample-key merge, eight run names, notebook structure,
English-only active content, plot-capture contract, and output audit. Final numerical
verification requires the HPC study files and the edgeR-enabled environment because the
MiXCR exports are not stored in this repository.
