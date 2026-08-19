# Execution guide

## 1. Update the repository safely

From the HPC checkout:

```bash
cd ~/mice_transplant_env_test
git status --short
git switch reproducible-article-environment
git pull --ff-only origin reproducible-article-environment
git log -1 --oneline
```

If `git status --short` lists modified tracked files, save them before pulling:

```bash
git stash push -u -m "HPC notebooks and outputs before repository update"
git pull --ff-only origin reproducible-article-environment
git stash list
```

Do not restore an old full-tree stash automatically after the update. The legacy notebook
layout differs from the active layout and may produce path conflicts.

## 2. Update the environment

```bash
conda env update \
  -n mice-transplant-2025 \
  -f environment.yml \
  --prune
conda activate mice-transplant-2025
```

The environment contains the pinned `repseq` revision, R, rpy2, and edgeR required by
Approach 1.

## 3. Verify the direct inputs for Approach 1

```bash
test -s /projects/mice_transplant_2025/metadata_mice_transplant.csv \
  && echo "Metadata CSV: OK"

test -s /projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv \
  && echo "Clonoset index: OK"
```

The clonoset-index CSV contains the MiXCR export paths. Approach 1 validates every
referenced file before creating the aaV count table. Parquet is not used by this method.

## 4. Run Approach 1 across all eight strata

```bash
python scripts/run_analysis.py --approach 1
```

The command uses the documented HPC paths by default. To provide different locations:

```bash
python scripts/run_analysis.py \
  --approach 1 \
  --metadata-csv /absolute/path/to/metadata_mice_transplant.csv \
  --clonoset-index /absolute/path/to/clonosets_mice_transplant_2025_df.csv
```

The eight analyses are executed in this order:

1. `cd4_thymus`
2. `cd8_thymus`
3. `cd4_spleen`
4. `cd8_spleen`
5. `cd4_combined`
6. `cd8_combined`
7. `thymus_combined`
8. `spleen_combined`

The four pooled strata use mouse-level analysis units.

## 5. Run selected Approach 1 strata

One stratum:

```bash
python scripts/run_analysis.py \
  --approach 1 \
  --strata cd4_thymus
```

Several strata:

```bash
python scripts/run_analysis.py \
  --approach 1 \
  --strata cd4_thymus,cd8_thymus,thymus_combined
```

The consolidated table contains the strata completed by that invocation. Run all eight
strata for the final manuscript-level comparison.

## 6. Monitor progress

Follow the complete runner log:

```bash
tail -f results/logs/pipeline.log
```

Follow Approach 1 at cell level:

```bash
tail -f results/logs/01_set_count.log
```

Each entry reports the cell position, stable cell ID, analysis step, elapsed time, and
remaining cell count:

```text
[CELL 003/010] id=a1-stratum-1 | step=Analyze CD4 T cells: thymus | START
[CELL 003/010] id=a1-stratum-1 | step=Analyze CD4 T cells: thymus | RUNNING | elapsed=300s
[CELL 003/010] id=a1-stratum-1 | step=Analyze CD4 T cells: thymus | DONE in 412.7s
```

## 7. Resume after interruption

```bash
python scripts/run_analysis.py --approach 1 --resume
```

The runner reloads the executed-notebook checkpoint, rebuilds the bootstrap state, and
skips completed analytical cells with matching IDs.

## 8. Locate Approach 1 results

| Output | Location |
|---|---|
| Per-stratum tables and conclusions | `results/01_set_count/` |
| Eight-stratum overview | `results/01_set_count/eight_stratum_summary.csv` |
| Consolidated distinctive TRAV table | `results/01_set_count/distinctive_trav_summary.csv` |
| Executed notebook | `results/executed_notebooks/01_set_count.executed.ipynb` |
| Cell-level log | `results/logs/01_set_count.log` |
| Figures | `figures/01_set_count/<stratum>/` |
| Figure inventory | `results/figure_inventory.csv` |

Each stratum produces ten plot types in both PNG and PDF format. A statistically null
result still produces a labeled figure panel, which keeps the eight output sets complete
and directly comparable.

## 9. Run Approaches 2 to 4

Approaches 2 and 3 use the later derived sample-resolved table:

```text
/absolute/path/to/derived_data/clean_clonotypes_aaV.parquet
```

Run Approach 2:

```bash
python scripts/run_analysis.py \
  --approach 2 \
  --data-dir /absolute/path/to/derived_data
```

Run Approach 3:

```bash
python scripts/run_analysis.py \
  --approach 3 \
  --data-dir /absolute/path/to/derived_data
```

Compare completed outputs from Approaches 1 to 3:

```bash
python scripts/run_analysis.py --approach 4
```

Run all methods in publication order:

```bash
python scripts/run_analysis.py \
  --approach all \
  --data-dir /absolute/path/to/derived_data
```

Run a selected method set:

```bash
python scripts/run_analysis.py \
  --approach 1,3 \
  --data-dir /absolute/path/to/derived_data
```

## 10. Understand the active notebooks

The four source notebooks are stored directly under `notebooks/` and are numbered by
execution order. The former nested `approaches/` source directory is no longer used.

An older executed `venn_original.ipynb` may remain in a local `audit_runs/` directory or
Git stash. Its Russian Markdown and inline figures are historical outputs and are not
modified by pulling the active branch. The current first-pass notebook is:

```text
notebooks/01_set_count_analysis.ipynb
```

## 11. Validate the checkout

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

These checks validate structure, English-only active content, all eight strata, notebook
cell IDs, centralized figure persistence, Python syntax, and regression behavior. Full
numerical validation still requires the study data and the edgeR-enabled Conda environment.

## 12. Control CPU use

```bash
ARTICLE_N_JOBS=8 python scripts/run_analysis.py --approach 1
```

Do not request more workers than the scheduler allocation. Additional scheduler-specific
guidance is provided in [`docs/parallel_execution.md`](docs/parallel_execution.md).
