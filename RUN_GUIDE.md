# Execution guide

## 1. Prepare the environment

From the repository root:

```bash
conda env create -f environment.yml
conda activate mice-transplant-2025
python scripts/validate_repository.py
python -m pytest -q
```

For an existing environment:

```bash
conda env update -n mice-transplant-2025 -f environment.yml --prune
conda activate mice-transplant-2025
```

## 2. Verify the data used by the original first pass

`venn_original.ipynb` uses the existing HPC project layout:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
/projects/mice_transplant_2025/test_run/mixcr/
```

Check the two direct tabular inputs before execution:

```bash
test -s /projects/mice_transplant_2025/metadata_mice_transplant.csv
test -s /projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
```

The second CSV is a sample and file index. Its `filename` column points to the exported MiXCR clonotype tables. The notebook passes this index to repseq, which builds the exact aaV count tables during execution.

No Parquet file is required for this first pass. The raw FASTQ directory is needed only when `presenting_repseq.ipynb` is used to regenerate MiXCR outputs.

## 3. Run the verified first pass on HPC

From the audit checkout:

```bash
cd ~/mice_transplant_env_test
conda activate mice-transplant-2025
mkdir -p audit_runs logs

jupyter nbconvert \
  --to notebook \
  --execute venn_original.ipynb \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=-1 \
  --output-dir audit_runs \
  --output 01_venn_original.executed.ipynb \
  2>&1 | tee logs/01_venn_original.log
```

The executed notebook is written to:

```text
audit_runs/01_venn_original.executed.ipynb
```

The complete console record is written to:

```text
logs/01_venn_original.log
```

If `venn_original.ipynb` is absent because the modular branch removed legacy root notebooks, restore the DESeq2-free version from repository history:

```bash
git show 6678766a0e8913ec2feacb3efb2daebd6a27eda4:venn_original.ipynb \
  > venn_original.ipynb
```

The notebook attached during the audit represents the immediately preceding DESeq2-containing state. It confirms the same CSV and MiXCR input route. The command above restores the subsequent version in which the DESeq2 section was removed.

## 4. Understand the modular Parquet interface

The concise notebooks under `approaches/` were added later. Modular Approaches 1 to 3 call `src.strata.load_repertoire()` and expect:

```text
data/clean_clonotypes_aaV.parquet
```

This is a derived, sample-resolved analysis table with the following fields:

```text
cdr3, v_gene, umi, group, sample_id, mouse_id, source, subtype
```

The repository currently has no verified command that creates this table from the MiXCR clonoset index. Therefore the modular runner is valid only when the derived table has already been materialized. It is not the required route for rerunning `venn_original.ipynb`.

When the derived table exists, its directory can be supplied with:

```bash
export MICE_TCR_DATA_DIR=/absolute/path/to/derived_data
```

or:

```bash
python scripts/run_analysis.py \
  --approach all \
  --data-dir /absolute/path/to/derived_data
```

The complete provenance audit is available in [docs/notebook_input_provenance.md](docs/notebook_input_provenance.md).

## 5. Understand the six modular strata

| Stratum | Independent analysis unit |
|---|---|
| `cd4_thymus` | One CD4 thymus sample |
| `cd4_spleen` | One CD4 spleen sample |
| `cd8_thymus` | One CD8 thymus sample |
| `cd8_spleen` | One CD8 spleen sample |
| `cd4_combined` | CD4 thymus and spleen pooled within each mouse |
| `cd8_combined` | CD8 thymus and spleen pooled within each mouse |

The combined analyses preserve mouse-level independence. CD4 and CD8 cells are always analyzed separately.

## 6. Run the modular analyses

Modular Approach 1: exact aaV sets, edgeR, and Fisher cross-check. Use this command only after the derived Parquet table has been prepared:

```bash
python scripts/run_analysis.py --approach 1
```

Approach 2: sequence-space density, RFF-MMD, PERMANOVA, and witness scoring:

```bash
python scripts/run_analysis.py --approach 2
```

Approach 3: mouse-level clone evidence and V-segment prioritization:

```bash
python scripts/run_analysis.py --approach 3
```

Cross-approach comparison after Approaches 1–3:

```bash
python scripts/run_analysis.py --approach 4
```

Complete publication pipeline:

```bash
python scripts/run_analysis.py --approach all
```

Selected approaches:

```bash
python scripts/run_analysis.py --approach 1,2
```

Selected strata:

```bash
python scripts/run_analysis.py \
  --approach 3 \
  --strata cd4_spleen,cd8_spleen
```

## 7. Monitor modular progress

The master log records approach-level progress:

```bash
tail -f outputs/logs/pipeline.log
```

Each approach also has a cell-level log:

```bash
tail -f outputs/logs/02_sequence_embedding.log
```

A cell entry contains:

```text
[CELL 003/009] id=... | step=... | START | 6 position(s) remain
[CELL 003/009] id=... | step=... | RUNNING | elapsed=300s
[CELL 003/009] id=... | step=... | DONE in 412.7s | 6 code cell(s) not yet checkpointed
```

The executed notebook is checkpointed after every successful code cell.

## 8. Resume a modular run after failure

Use the same approach and strata with `--resume`:

```bash
python scripts/run_analysis.py \
  --approach 2 \
  --strata cd4_thymus,cd8_thymus \
  --resume
```

The runner restores outputs from completed cells, re-executes bootstrap cells to rebuild Python state, and skips completed analytical cells. Approach 2 reuses its cached coordinate basis only when the stored clonotype fingerprint matches the current input.

## 9. Locate results

| Result | Location |
|---|---|
| Statistical tables | `outputs/tables/<approach>/` |
| Per-stratum conclusions | `outputs/tables/<approach>/<stratum>_conclusion.md` |
| Executed notebooks | `outputs/notebooks/` |
| Logs | `outputs/logs/` |
| Cached heavy objects | `outputs/cache/` |
| Figures | `figures/<approach>/<stratum>/` |

Every matplotlib figure is saved as PNG and PDF. The validation script rejects active plotting code that has no persistent figure output.

## 10. Interpret the cross-approach output

Approach-specific V-segment ranks are compared within the same biological stratum. The comparison reports:

1. pairwise rank correlation;
2. top-10 overlap;
3. signed effect correlation on a common g1-versus-allogeneic axis;
4. directional agreement;
5. direction-aware consensus rank.

A segment can therefore receive repeated ranking support while retaining a mixed effect direction. The pipeline reports that disagreement explicitly.

## 11. Control CPU usage

The runtime respects Linux CPU affinity and Slurm allocation variables. To reduce CPU use within an allocation:

```bash
ARTICLE_N_JOBS=8 python scripts/run_analysis.py --approach 2
```

Do not set `ARTICLE_N_JOBS` above the scheduler allocation. See [docs/parallel_execution.md](docs/parallel_execution.md) for monitoring commands and memory policy.

## 12. Final validation

After code or notebook changes:

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

A full scientific validation requires the correct notebook-specific input. For the historical first pass, this means the metadata CSV, the clonoset-index CSV, and its referenced MiXCR exports. For the modular pipeline, this means a traceably derived sample-resolved Parquet table. Regenerated tables, figures, and conclusions should be reviewed together before manuscript use.
