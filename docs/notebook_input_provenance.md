# Notebook input provenance

This file records which data each stage reads and where the derived Parquet interface enters the workflow.

## Input map

| Notebook or stage | Direct inputs | Parquet use |
|---|---|---|
| `presenting_repseq.ipynb` | Raw paired FASTQ files, `metadata_mice_transplant.csv`, MiXCR outputs, `table_report.csv` | None in the primary preparation block |
| Historical `venn_original.ipynb` | `metadata_mice_transplant.csv`, `clonosets_mice_transplant_2025_df.csv`, referenced MiXCR exports | None |
| Historical `mirpy_analysis.ipynb` | Flattened clonotype tables including `clean_clonotypes_aaV.parquet` and aaVJ count tables | Derived analysis input |
| Historical `results_summary.ipynb` | CSV tables from upstream analyses | None |
| Historical Study 2 notebook | `study2_data/clonosets_study2.csv`, repseq-compatible tables and derived CSV/Parquet layers | Derived analysis input |
| Active Approach 1 | Metadata, clonoset index and referenced MiXCR exports | None |
| Active Approaches 2 and 3 | `clean_clonotypes_aaV.parquet` through `src.strata.load_repertoire()` | Required |
| Active Approach 4 | Standardized outputs from Approaches 1 to 3 | None |

The historical audit used `presenting_repseq.ipynb`, the edgeR/Fisher version of `venn_original.ipynb` from commit `6678766a0e8913ec2feacb3efb2daebd6a27eda4`, `mirpy_analysis.ipynb`, `results_summary.ipynb`, the Study 2 alloreactivity notebook and the active notebooks in this repository.

## Approach 1 inputs

The default paths are:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
/projects/mice_transplant_2025/test_run/mixcr/
```

`clonosets_mice_transplant_2025_df.csv` maps samples to exported clonotype files. The `filename` entries must resolve to existing MiXCR tables. `repseq.intersections.count_table()` reads those tables and builds the aaV count matrices in memory.

Raw FASTQ files are only needed when the MiXCR outputs have to be regenerated. A normal Approach 1 rerun starts from the existing exports and clonoset index.

The active notebook keeps the original chain-plus-sample-number metadata merge, the historical TRA exclusions and exact aaV matching with zero mismatches. It writes one executed notebook per biological stratum under `audit_runs/`.

## Derived-data boundary

Approaches 2 and 3 require `clean_clonotypes_aaV.parquet`. There is currently no checked fresh-start command in this repository that rebuilds that file from the clonoset index. For that reason, the Parquet file is treated as an existing derived input for those approaches.

A future materialization step would need to preserve sample ID, mouse ID, group, tissue, CD4/CD8 subtype, chain, CDR3 amino-acid sequence, V segment and abundance.
