# Notebook input provenance

This document records the data interface used by each notebook in the original HPC workflow and distinguishes source data from derived analysis artifacts.

## Audit basis

The input audit used the following sources:

- `presenting_repseq.ipynb`, which prepares MiXCR outputs and the clonoset index;
- the edgeR-and-Fisher `venn_original.ipynb` produced in commit `6678766a0e8913ec2feacb3efb2daebd6a27eda4`;
- `mirpy_analysis.ipynb` and `results_summary.ipynb` from the original three-method study;
- `study2_alloreactivity/alloreactivity_study.ipynb` from the original clone-level study;
- the restored root notebook and the three later source notebooks under `notebooks/`.

## Input map

| Notebook or stage | Direct inputs | Role of Parquet |
|---|---|---|
| `presenting_repseq.ipynb` | Raw paired FASTQ files, `metadata_mice_transplant.csv`, MiXCR outputs, `table_report.csv` | No Parquet input is used in the primary preparation block. |
| Historical `venn_original.ipynb` | `metadata_mice_transplant.csv`, `clonosets_mice_transplant_2025_df.csv`, and the MiXCR export files referenced by the clonoset index | No Parquet file is read. |
| Historical `mirpy_analysis.ipynb` | Derived flattened clonotype tables including `clean_clonotypes_aaV.parquet` and aaVJ count tables | Parquet is a derived analysis format used by this notebook. It is not an input to the first pass. |
| Historical `results_summary.ipynb` | CSV result tables produced by the upstream analyses | No primary repertoire input is read. |
| Historical Study 2 notebook | `study2_data/clonosets_study2.csv`, repseq-compatible clonotype tables, and previously generated CSV/Parquet result layers | Several Parquet files are derived results consumed by later sections. |
| Active Approach 1 | `metadata_mice_transplant.csv`, `clonosets_mice_transplant_2025_df.csv`, and its referenced MiXCR exports, read directly in `venn_original.ipynb` | No Parquet input is used. |
| Active Approaches 2 and 3 | `clean_clonotypes_aaV.parquet` through `src.strata.load_repertoire()` | Parquet remains the later derived interface for sequence-space and clone-evidence analyses. |
| Active Approach 4 | Standardized tables produced by active Approaches 1 to 3 | No repertoire file is read directly. |

## Verified first-pass data contract

The original first pass uses the following paths inside `venn_original.ipynb`:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
/projects/mice_transplant_2025/test_run/mixcr/
```

`clonosets_mice_transplant_2025_df.csv` is an index of samples and exported clonotype files. Its `filename` values must resolve to existing MiXCR clonotype tables. The notebook passes this index to `repseq.intersections.count_table()`, which reads the referenced clonotype tables and constructs aaV count matrices in memory.

The raw FASTQ directory is used by `presenting_repseq.ipynb` when MiXCR outputs need to be regenerated. It is not read during a normal rerun of `venn_original.ipynb` when the MiXCR exports and clonoset index already exist.

## Active first-pass implementation

`python scripts/run_analysis.py --approach 1` executes the restored notebook once for
each biological stratum. The notebook reproduces the original
chain-plus-sample-number metadata merge, applies the historical TRA exclusions, and
invokes `repseq.intersections.count_table()` in aaV mode with zero mismatches. No
minimum abundance threshold is applied. The runner produces eight separately executed
notebooks under `audit_runs/`.

## Remaining derived-data boundary

Approaches 2 and 3 still require `clean_clonotypes_aaV.parquet`. The repository does not currently contain a verified materialization command that recreates this later derived table from the clonoset index. The Parquet requirement therefore applies only to those later methods. A future materialization step must preserve sample identifiers, mouse identifiers, group, tissue, CD4/CD8 subtype, chain, CDR3 amino-acid sequence, V segment, and abundance before the complete four-method workflow can be treated as a fresh-start pipeline.
