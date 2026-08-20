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
| TRA/TRB sets | Retain g1 to g6 count tables, five-circle set diagrams, exact subtraction, V-segment retention tables, bubble plots, heatmaps, and intragroup Venn diagrams. |
| Repertoire similarity | Retain mouse-level TRA/TRB Jensen-Shannon, rank-correlation, and UMI-correlation analyses. |
| Fisher | Compare pooled g1 and g5+g6 feature counts with the historical pseudocount-adjusted log2 ratio and Benjamini-Hochberg correction. |
| edgeR | Retain TMM normalization, `filterByExpr`, quasi-likelihood fitting, and the g1 coefficient relative to pooled g5+g6 annotation. |
| Thresholds | Introduce no minimum UMI or read-count threshold. |

The active differential section contains edgeR and the Fisher sensitivity analysis
only. No additional differential-expression model is fitted.

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

After all eight executions, the final run consolidates the stratum-level TRAV tables
and builds the cross-stratum evidence heatmap. Descriptive exact-set evidence is kept
separate from FDR-supported inference.

## Verification boundary

Repository tests verify the sample-key merge, eight run names, notebook structure,
English-only active content, plot-capture contract, and output audit. Final numerical
verification requires the HPC study files and the edgeR-enabled environment because the
MiXCR exports are not stored in this repository.
