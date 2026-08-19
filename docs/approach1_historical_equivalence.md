# Approach 1 historical-equivalence audit

This audit records which operations in the active count-based pipeline reproduce the
last executable `venn_original.ipynb` workflow and which operations extend it to the
eight prespecified biological strata.

## Operations preserved from the historical notebook

| Stage | Historical operation retained in the active pipeline |
|---|---|
| Input files | Read `metadata_mice_transplant.csv`, `clonosets_mice_transplant_2025_df.csv`, and the MiXCR exports referenced by the clonoset index. |
| Legacy sample key | Construct `<metadata chain>-<sample_no>` before chain standardization. The study files therefore match keys such as `alpha-100` and `beta-100`. |
| Metadata merge | Replace the clonoset-index legacy identifier with the article-level `sample_id` and retain the chain reported by the clonoset index. |
| Chain and groups | Select TRA repertoires from g1, g2, g5, and g6 for count-table construction. g2 remains available in the count table but is outside the primary g1 versus g5 + g6 contrast. |
| Exclusions | Exclude `g1_m3_thymus_cd8_80_alpha`, `g5_m1_spleen_cd4_93_alpha`, and `g5_m4_thymus_cd8_100_alpha`. |
| Feature definition | Call `repseq.intersections.count_table` with `overlap_type="aaV"` and `mismatches=0`. |
| Abundance filtering | Apply no minimum UMI threshold. Invalid amino-acid sequences and non-TRAV features are excluded before modeling. |
| Fisher test | Pool g5 and g6 counts, use the historical two-sided 2 x 2 contingency table, calculate `log2((g1 + 1)/(g5 + g6 + 1))`, and control FDR by Benjamini-Hochberg. |
| edgeR | Pool g5 and g6 in the group annotation, retain separate count columns, run TMM normalization, `filterByExpr`, quasi-likelihood fitting, and coefficient 2 for g1 relative to the allogeneic reference. |
| Active comparison | Report g1 versus pooled g5 + g6. DESeq2 is excluded from the active workflow. |

## Required extensions

The historical notebook evaluated the combined TRA count table. The active workflow
applies the same count-based analysis independently to eight biological strata. The
four single-compartment strata retain sample-level count columns. The four combined
strata sum counts within the same mouse before inference so that paired tissues or cell
subsets do not become independent replicates.

The active implementation also adds deterministic filenames, complete PNG and PDF
figure persistence, per-stratum conclusions, a consolidated TRAV table, durable cell
logs, checkpoint recovery, and explicit input validation. These changes affect
orchestration and reporting; they do not introduce an abundance threshold or change
the primary biological contrast.

## Verification boundary

Unit and synthetic tests verify the merge contract, feature conversion, historical
Fisher effect definition, eight-stratum selection, output persistence, and repository
structure. Final numerical equivalence for the original unstratified contrast must be
checked on the HPC study files because those files are not stored in this repository.
