# Mouse TCR repertoire analysis after transplantation

This repository contains a reproducible analysis of mouse T-cell receptor repertoires
after transplantation. The primary exact-set and count-based analysis is implemented
directly in [`venn_original.ipynb`](venn_original.ipynb).

The notebook was restored from the verified pre-refactor implementation at commit
`6678766a0e8913ec2feacb3efb2daebd6a27eda4`. It retains the original TRA/TRB exact
`aaV` intersections, V-segment retention profiles, mouse-level similarity analyses,
UMI correlations, edgeR model, and Fisher sensitivity analysis. It does not use a
derived Parquet table.

## Analytical approaches

| Order | Source notebook | Purpose |
|---|---|---|
| 1 | [`venn_original.ipynb`](venn_original.ipynb) | Exact aaV sets, TRA/TRB repertoire structure, edgeR, Fisher, and TRAV prioritization |
| 2 | [`02_sequence_embedding_analysis.ipynb`](notebooks/02_sequence_embedding_analysis.ipynb) | Sequence-space enrichment and repertoire geometry |
| 3 | [`03_clone_alloreactivity_analysis.ipynb`](notebooks/03_clone_alloreactivity_analysis.ipynb) | Mouse-normalized clone-level alloreactivity evidence |
| 4 | [`04_cross_approach_comparison.ipynb`](notebooks/04_cross_approach_comparison.ipynb) | Direction-aware comparison of the completed approaches |

Approach 1 is notebook-first by design. One parameterized source notebook is executed
independently for each biological stratum, producing eight self-contained audit
notebooks rather than one notebook containing eight opaque pipeline calls.

## Eight biological strata

| Order | Stratum key | Biological subset | Executed notebook |
|---|---|---|---|
| 1 | `cd4_thymus` | CD4, thymus | `audit_runs/01_cd4_thymus.executed.ipynb` |
| 2 | `cd8_thymus` | CD8, thymus | `audit_runs/02_cd8_thymus.executed.ipynb` |
| 3 | `cd4_spleen` | CD4, spleen | `audit_runs/03_cd4_spleen.executed.ipynb` |
| 4 | `cd8_spleen` | CD8, spleen | `audit_runs/04_cd8_spleen.executed.ipynb` |
| 5 | `cd4_combined` | CD4, thymus and spleen pooled within mouse | `audit_runs/05_cd4_thymus_spleen.executed.ipynb` |
| 6 | `cd8_combined` | CD8, thymus and spleen pooled within mouse | `audit_runs/06_cd8_thymus_spleen.executed.ipynb` |
| 7 | `thymus_combined` | Thymus, CD4 and CD8 pooled within mouse | `audit_runs/07_thymus_cd4_cd8.executed.ipynb` |
| 8 | `spleen_combined` | Spleen, CD4 and CD8 pooled within mouse | `audit_runs/08_spleen_cd4_cd8.executed.ipynb` |

Stratification is applied to the sample index before any count table is constructed.
For the four combined strata, counts are pooled within mouse before the edgeR model so
paired samples are not treated as independent replicates.

## Input provenance for Approach 1

The default HPC inputs are the same CSV and MiXCR route used by the historical
notebook:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
MiXCR clonotype exports referenced by the clonoset-index CSV
```

The historical `alpha-N` and `beta-N` identifiers are matched to metadata before chain
normalization. Exact clonotypes are defined as `(CDR3 amino-acid sequence, V segment)`
with `overlap_type="aaV"` and `mismatches=0`. The three prespecified historical sample
exclusions are retained. No minimum UMI threshold is introduced.

Approaches 2 and 3 use the later derived `clean_clonotypes_aaV.parquet` interface. That
file is not required for Approach 1.

## Repository layout

```text
mice_transplant_2025/
├── venn_original.ipynb        # Complete source notebook for Approach 1
├── notebooks/                 # Approaches 2 to 4
├── src/                       # Shared code for later approaches and runtime control
├── scripts/                   # Execution and validation commands
├── audit_runs/                # Eight generated executed notebooks
├── logs/                      # Pipeline and per-notebook cell logs
├── figures/
│   └── 01_set_count/
│       ├── <stratum>/         # Every figure displayed for that stratum
│       └── cross_stratum/     # Eight-stratum TRAV comparison
└── results/
    └── 01_set_count/
        ├── <stratum>/         # Tables, manifest, and conclusion
        └── eight_stratum_distinctive_trav_summary.csv
```

`figures/` is the only figure root. A notebook-level figure hook saves every explicit
and displayed Matplotlib figure as PNG and PDF. Each stratum writes a manifest, and the
runner fails if a listed file is missing.

## Quick start on HPC

Update the branch and environment:

```bash
cd ~/mice_transplant_env_test
git status --short
git switch reproducible-article-environment
git pull --ff-only origin reproducible-article-environment

conda env update \
  -n mice-transplant-2025 \
  -f environment.yml \
  --prune
conda activate mice-transplant-2025
```

Run the complete first approach across all eight strata:

```bash
python scripts/run_analysis.py --approach 1
```

Run only one stratum:

```bash
python scripts/run_analysis.py \
  --approach 1 \
  --strata cd4_thymus
```

Monitor the first stratum:

```bash
tail -f logs/01_cd4_thymus.log
```

Each log entry reports the current cell, total cell count, elapsed time, and remaining
cells. See [`RUN_GUIDE.md`](RUN_GUIDE.md) for safe handling of local notebook changes,
all commands, path overrides, and direct `nbconvert` execution.

Replay a checkpointed run after interruption:

```bash
python scripts/run_analysis.py --approach 1 --resume
```

## Output contract

Each stratum produces:

- one executed notebook under `audit_runs/`;
- one cell-level log under `logs/`;
- all displayed figures under `figures/01_set_count/<stratum>/` in PNG and PDF;
- exact-set, similarity, UMI-correlation, Fisher, and edgeR tables under
  `results/01_set_count/<stratum>/`;
- a complete TRAV ranking and a concise stratum-specific conclusion;
- a figure manifest verified by the runner.

After all eight runs, the final notebook creates
`results/01_set_count/eight_stratum_distinctive_trav_summary.csv` and the cross-stratum
TRAV evidence heatmap. A null inferential result is reported as such; descriptive
exact-set rankings remain visible and are never relabeled as statistically supported.

Run the later approaches only when their derived input is available:

```bash
python scripts/run_analysis.py --approach 2 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 3 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 4
python scripts/run_analysis.py --approach all --data-dir /absolute/path/to/derived_data
```

## Validation

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

Repository validation checks the restored notebook structure, all eight run names,
English-only active content, absence of committed notebook outputs, centralized figure
persistence, and output-verification logic. Full numerical execution requires the study
data and the edgeR-enabled Conda environment on HPC.
