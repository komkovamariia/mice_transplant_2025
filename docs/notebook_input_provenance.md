# Notebook input provenance

This document records the data interface used by each notebook in the original HPC workflow and distinguishes source data from derived analysis artifacts.

## Audit basis

The input audit used the following sources:

- `presenting_repseq.ipynb`, which prepares MiXCR outputs and the clonoset index;
- `venn_original.ipynb` immediately before removal of DESeq2 from the notebook;
- the DESeq2-free `venn_original.ipynb` produced in commit `6678766a0e8913ec2feacb3efb2daebd6a27eda4`;
- `mirpy_analysis.ipynb` and `results_summary.ipynb` from the original three-method study;
- `study2_alloreactivity/alloreactivity_study.ipynb` from the original clone-level study;
- the four modular notebooks under `approaches/`.

## Input map

| Notebook or stage | Direct inputs | Role of Parquet |
|---|---|---|
| `presenting_repseq.ipynb` | Raw paired FASTQ files, `metadata_mice_transplant.csv`, MiXCR outputs, `table_report.csv` | No Parquet input is used in the primary preparation block. |
| Historical `venn_original.ipynb` | `metadata_mice_transplant.csv`, `clonosets_mice_transplant_2025_df.csv`, and the MiXCR export files referenced by the clonoset index | No Parquet file is read. |
| Historical `mirpy_analysis.ipynb` | Derived flattened clonotype tables including `clean_clonotypes_aaV.parquet` and aaVJ count tables | Parquet is a derived analysis format used by this notebook. It is not an input to the first pass. |
| Historical `results_summary.ipynb` | CSV result tables produced by the upstream analyses | No primary repertoire input is read. |
| Historical Study 2 notebook | `study2_data/clonosets_study2.csv`, repseq-compatible clonotype tables, and previously generated CSV/Parquet result layers | Several Parquet files are derived results consumed by later sections. |
| Modular Approaches 1 to 3 | `clean_clonotypes_aaV.parquet` through `src.strata.load_repertoire()` | This interface was introduced by the modular refactor. |
| Modular Approach 4 | Standardized tables produced by modular Approaches 1 to 3 | No repertoire file is read directly. |

## Verified first-pass data contract

The original first pass uses the following paths inside `venn_original.ipynb`:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
/projects/mice_transplant_2025/test_run/mixcr/
```

`clonosets_mice_transplant_2025_df.csv` is an index of samples and exported clonotype files. Its `filename` values must resolve to existing MiXCR clonotype tables. The notebook passes this index to `repseq.intersections.count_table()`, which reads the referenced clonotype tables and constructs aaV count matrices in memory.

The raw FASTQ directory is used by `presenting_repseq.ipynb` when MiXCR outputs need to be regenerated. It is not read during a normal rerun of `venn_original.ipynb` when the MiXCR exports and clonoset index already exist.

## Parquet limitation in the modular branch

The modular notebooks currently call `load_repertoire()`, which expects `clean_clonotypes_aaV.parquet`. The repository does not currently contain a verified materialization command that creates this file from `clonosets_mice_transplant_2025_df.csv` and its MiXCR exports. Therefore:

1. the Parquet file must not be presented as a prerequisite for the historical first pass;
2. `python scripts/run_analysis.py --approach 1` is not equivalent to the verified `venn_original.ipynb` audit run unless that derived file has already been prepared;
3. the original HPC command remains the verified route for rerunning the first pass from the existing project data;
4. a future MiXCR-to-Parquet materialization step must preserve sample identifiers, mouse identifiers, group, tissue, CD4/CD8 subtype, chain, CDR3 amino-acid sequence, V segment, and UMI count before the modular runner is treated as a fresh-start workflow.

This distinction prevents a derived cache format from being confused with the study's primary data source.
