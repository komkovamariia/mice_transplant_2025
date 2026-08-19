# Mouse TCR repertoire analysis after transplantation

This repository contains a reproducible analysis of mouse T-cell receptor alpha-chain
(TRA) repertoires. The primary contrast is g1 versus the pooled g5+g6 allogeneic
reference. Exact clonotypes are defined by CDR3 amino-acid sequence and TRAV segment
(`aaV`).

## Analytical workflow

| Order | Source notebook | Scientific purpose |
|---|---|---|
| 1 | [`01_set_count_analysis.ipynb`](notebooks/01_set_count_analysis.ipynb) | Exact aaV set subtraction, edgeR differential counts, Fisher cross-check, and TRAV prioritization |
| 2 | [`02_sequence_embedding_analysis.ipynb`](notebooks/02_sequence_embedding_analysis.ipynb) | Local sequence-space enrichment and whole-repertoire geometry |
| 3 | [`03_clone_alloreactivity_analysis.ipynb`](notebooks/03_clone_alloreactivity_analysis.ipynb) | Mouse-normalized clone-level alloreactivity evidence |
| 4 | [`04_cross_approach_comparison.ipynb`](notebooks/04_cross_approach_comparison.ipynb) | Direction-aware comparison of the first three methods |

The word *approach* refers to an independent analytical method. Numbering defines the
recommended execution order and does not represent repeated runs of the same model.

## Eight biological strata

Each stratum is analyzed separately and receives its own tables, figures, and conclusion.
Pooling is performed within mouse so paired tissues or cell subsets are not treated as
independent biological replicates.

| Stratum | Independent analysis unit |
|---|---|
| `cd4_thymus` | One CD4 thymus sample |
| `cd8_thymus` | One CD8 thymus sample |
| `cd4_spleen` | One CD4 spleen sample |
| `cd8_spleen` | One CD8 spleen sample |
| `cd4_combined` | CD4 thymus and spleen pooled within mouse |
| `cd8_combined` | CD8 thymus and spleen pooled within mouse |
| `thymus_combined` | Thymus CD4 and CD8 counts pooled within mouse |
| `spleen_combined` | Spleen CD4 and CD8 counts pooled within mouse |

## Repository layout

```text
mice_transplant_2025/
├── notebooks/                 # Four numbered source notebooks
├── src/                       # Tested analytical implementation
├── scripts/                   # Runner and repository validation
├── figures/                   # The only figure root
│   └── <numbered_method>/<stratum>/
├── results/                   # Tables, conclusions, logs, and executed notebooks
│   ├── <numbered_method>/
│   ├── executed_notebooks/
│   ├── logs/
│   └── cache/
├── docs/                      # Provenance and execution notes
└── tests/                     # Regression tests
```

Generated figures never appear under `notebooks/`, `results/`, or method-specific source
directories. Every saved plot is indexed in `results/figure_inventory.csv`.

## Input provenance

Approach 1 directly reproduces the verified first-pass data route:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
MiXCR exports referenced by the clonoset-index CSV
```

It constructs the aaV count table in memory with the pinned `repseq` implementation. It
does not require a Parquet file. No minimum UMI or read-count threshold is applied;
non-functional amino-acid sequences and non-TRAV records are excluded during input
validation.

Approaches 2 and 3 use the later derived sample-resolved interface
`clean_clonotypes_aaV.parquet`. Approach 4 reads the standardized result tables created by
Approaches 1 to 3. Full details are provided in
[`docs/notebook_input_provenance.md`](docs/notebook_input_provenance.md). The exact
correspondence between the historical first pass and active Approach 1 is documented in
[`docs/approach1_historical_equivalence.md`](docs/approach1_historical_equivalence.md).

## Quick start on HPC

Update the active branch and environment:

```bash
git switch reproducible-article-environment
git pull --ff-only origin reproducible-article-environment

conda env update \
  -n mice-transplant-2025 \
  -f environment.yml \
  --prune
conda activate mice-transplant-2025
```

Run Approach 1 across all eight strata:

```bash
python scripts/run_analysis.py --approach 1
```

Run one or several selected strata:

```bash
python scripts/run_analysis.py \
  --approach 1 \
  --strata cd4_thymus,cd8_thymus
```

Run the later methods when the derived Parquet table is available:

```bash
python scripts/run_analysis.py --approach 2 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 3 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 4
```

Run the complete workflow:

```bash
python scripts/run_analysis.py \
  --approach all \
  --data-dir /absolute/path/to/derived_data
```

See [`RUN_GUIDE.md`](RUN_GUIDE.md) for explicit path overrides, checkpoint recovery, log
monitoring, and validation commands.

Long runs can be resumed with `--resume`. The complete progress stream is written to
`results/logs/pipeline.log`.

## Approach 1 outputs

Every biological stratum produces:

- a count matrix;
- complete edgeR and Fisher clonotype statistics;
- a TRAV ranking with evidence class and signed effect direction;
- a dedicated machine-readable JSON conclusion;
- a concise publication-oriented Markdown conclusion;
- ten matched PNG/PDF figure pairs.

The final notebook cell creates:

```text
results/01_set_count/eight_stratum_summary.csv
results/01_set_count/distinctive_trav_summary.csv
results/01_set_count/cross_stratum_conclusion.md
figures/01_set_count/distinctive_trav_across_eight_strata.png
figures/01_set_count/distinctive_trav_across_eight_strata.pdf
```

## Reproducibility safeguards

- Source notebooks contain no committed outputs and no inline-only plotting code.
- All Markdown, code comments, labels, titles, filenames, and repository documentation are
  written in English.
- edgeR is the prespecified sample-resolved model; Fisher testing is an independent pooled
  count check. DESeq2 is not part of the active workflow.
- Cell-level logs record start time, elapsed time, completion, failure, and remaining cell
  count. Successful cells are checkpointed for `--resume`.
- Static validation rejects Cyrillic text, figures outside `figures/`, multiple figure
  roots, missing strata, source notebook outputs, and plotting calls without persistent
  figure output.
