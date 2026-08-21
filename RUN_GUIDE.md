# Execution guide

## 1. Update the branch without losing local work

From the HPC checkout:

```bash
cd ~/mice_transplant_env_test
git status --short
```

If the status is empty, update directly:

```bash
git switch reproducible-article-environment
git pull --ff-only origin reproducible-article-environment
git log -1 --oneline
```

If tracked notebooks or generated directories are listed, save them before pulling:

```bash
git stash push -u -m "HPC notebooks and outputs before restored notebook update"
git switch reproducible-article-environment
git pull --ff-only origin reproducible-article-environment
git stash list
```

Do not immediately restore an old full-tree stash. First inspect it:

```bash
git stash show --stat stash@{0}
```

Old modified copies of `venn_original.ipynb` can overwrite the restored source notebook.
Recover only specific personal files if they are still needed.

## 2. Update and activate the environment

```bash
conda env update \
  -n mice-transplant-2025 \
  -f environment.yml \
  --prune
conda activate mice-transplant-2025
```

The environment includes the pinned `repseq` revision, R, rpy2, and edgeR required by
the restored first approach.

## 3. Verify the direct inputs

```bash
test -s /projects/mice_transplant_2025/metadata_mice_transplant.csv \
  && echo "Metadata CSV: OK"

test -s /projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv \
  && echo "Clonoset index: OK"
```

Approach 1 reads these CSV files and the MiXCR exports referenced by the second file.
It does not read Parquet data.

## 4. Run all eight independent notebooks

```bash
python scripts/run_analysis.py --approach 1
```

The command executes `venn_original.ipynb` eight times in this order:

1. `cd4_thymus`
2. `cd8_thymus`
3. `cd4_spleen`
4. `cd8_spleen`
5. `cd4_combined`
6. `cd8_combined`
7. `thymus_combined`
8. `spleen_combined`

The generated audit notebooks are:

```text
audit_runs/01_cd4_thymus.executed.ipynb
audit_runs/02_cd8_thymus.executed.ipynb
audit_runs/03_cd4_spleen.executed.ipynb
audit_runs/04_cd8_spleen.executed.ipynb
audit_runs/05_cd4_thymus_spleen.executed.ipynb
audit_runs/06_cd8_thymus_spleen.executed.ipynb
audit_runs/07_thymus_cd4_cd8.executed.ipynb
audit_runs/08_spleen_cd4_cd8.executed.ipynb
```

The command ends with an output audit. It fails if an executed notebook, result table,
manifest, PNG, or PDF is missing.

## 5. Run only one stratum

```bash
python scripts/run_analysis.py \
  --approach 1 \
  --strata cd4_thymus
```

Run several selected strata:

```bash
python scripts/run_analysis.py \
  --approach 1 \
  --strata cd4_thymus,cd8_thymus,thymus_combined
```

The full cross-stratum summary is finalized only when all eight strata are requested in
one invocation.

## 6. Override the default HPC paths

```bash
python scripts/run_analysis.py \
  --approach 1 \
  --metadata-csv /absolute/path/to/metadata_mice_transplant.csv \
  --clonoset-index /absolute/path/to/clonosets_mice_transplant_2025_df.csv \
  --working-dir /absolute/path/to/test_run
```

Relative MiXCR filenames are resolved against the clonoset-index directory and its
`mixcr/` subdirectory.

## 7. Monitor cell-level progress

Follow the complete pipeline:

```bash
tail -f logs/pipeline.log
```

Follow one executed notebook:

```bash
tail -f logs/01_cd4_thymus.log
```

The log format includes both progress and remaining cells:

```text
[CELL 012/032] id=venn-historical-code-17 | step=Historical notebook analysis block 17 | START | 20 code cell(s) remain after this cell
[CELL 012/032] id=venn-historical-code-17 | step=Historical notebook analysis block 17 RUNNING | elapsed=300s
[CELL 012/032] id=venn-historical-code-17 | step=Historical notebook analysis block 17 | DONE in 418.4s | 20 code cell(s) not yet checkpointed
```

## 8. Resume after interruption

```bash
python scripts/run_analysis.py --approach 1 --resume
```

A Jupyter checkpoint contains cell outputs but not the live Python kernel. Therefore a
resume safely replays completed cells to reconstruct notebook state, then continues
with durable checkpointing. This is slower than serializing opaque Python objects but
preserves notebook semantics.

## 9. Locate figures and tables

| Output | Location |
|---|---|
| Executed notebooks | `audit_runs/` |
| Cell-level logs | `logs/` |
| Per-stratum figures | `figures/01_set_count/<stratum>/` |
| Complete figure archive | `figures/01_set_count.zip` |
| Per-stratum tables and conclusion | `results/01_set_count/<stratum>/` |
| Figure manifest | `results/01_set_count/<stratum>/figure_manifest.csv` |
| Consolidated TRAV table | `results/01_set_count/eight_stratum_distinctive_trav_summary.csv` |

Every figure displayed by Matplotlib is saved in PNG and PDF. Explicit legacy paths are
redirected to the same central figure root. Non-figure tables remain under `results/`.
The TRA and TRB V-segment heatmaps contain only the ten segments with the highest
article score `S(v) = f_full(v) × r(v)` and only the four g1-containing regions. The
ranking candidates satisfy at least 0.6% initial share and at least 66% retention. The
normalized entropy therefore uses `log2(4)`. These are the only two heatmaps produced
per stratum. Auxiliary pairwise and mouse-level Venn figures are omitted. The final
output audit fails explicitly if edgeR does not complete or any additional heatmap is
present.

The runner clears the generated figure and result directories for each requested
stratum before execution, then creates `figures/01_set_count.zip` after validation.

## 10. Direct `nbconvert` execution

The runner is recommended because it provides cell-by-cell logs, consistent names, and
output verification. A direct single-stratum audit can still be run in the historical
style:

```bash
mkdir -p audit_runs logs

MICE_TCR_STRATUM=cd4_thymus \
MICE_TCR_FINALIZE_SUMMARY=0 \
jupyter nbconvert \
  --to notebook \
  --execute venn_original.ipynb \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=-1 \
  --output-dir audit_runs \
  --output 01_cd4_thymus.executed.ipynb \
  2>&1 | tee logs/01_cd4_thymus.nbconvert.log
```

This command uses the same CSV and MiXCR inputs. The notebook itself still saves all
figures and tables, but `nbconvert` does not provide the runner's per-cell progress log
or final output audit.

## 11. Run later approaches

Approaches 2 and 3 use the later sample-resolved Parquet interface:

```bash
python scripts/run_analysis.py \
  --approach 2 \
  --data-dir /absolute/path/to/derived_data

python scripts/run_analysis.py \
  --approach 3 \
  --data-dir /absolute/path/to/derived_data

python scripts/run_analysis.py --approach 4
```

Run every approach in publication order only when the derived input for Approaches 2
and 3 is available:

```bash
python scripts/run_analysis.py \
  --approach all \
  --data-dir /absolute/path/to/derived_data
```

## 12. Validate the checkout

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

The lightweight validation does not substitute for execution on the study data. It
checks notebook structure, English-only source text, eight-stratum naming, centralized
figure persistence, Python syntax, and runner output verification.

## 13. Control CPU use

```bash
ARTICLE_N_JOBS=8 python scripts/run_analysis.py --approach 1
```

Do not request more workers than the scheduler allocation.
