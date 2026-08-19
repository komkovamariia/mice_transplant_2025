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

## 2. Provide the canonical input

Default location:

```text
data/clean_clonotypes_aaV.parquet
```

Required columns:

```text
cdr3, v_gene, umi, group, sample_id, mouse_id, source, subtype
```

Use an external directory:

```bash
export MICE_TCR_DATA_DIR=/absolute/path/to/data
```

or specify it during execution:

```bash
python scripts/run_analysis.py \
  --approach all \
  --data-dir /absolute/path/to/data
```

The loader validates UMI values, functional CDR3 amino-acid sequences, TRA identity when a chain column is present, group assignments, mouse identifiers, cell subsets, tissues, and sample-level metadata consistency before analysis begins.

## 3. Understand the six strata

| Stratum | Independent analysis unit |
|---|---|
| `cd4_thymus` | One CD4 thymus sample |
| `cd4_spleen` | One CD4 spleen sample |
| `cd8_thymus` | One CD8 thymus sample |
| `cd8_spleen` | One CD8 spleen sample |
| `cd4_combined` | CD4 thymus and spleen pooled within each mouse |
| `cd8_combined` | CD8 thymus and spleen pooled within each mouse |

The combined analyses preserve mouse-level independence. CD4 and CD8 cells are always analyzed separately.

## 4. Run the analyses

Approach 1: exact aaV sets, edgeR, and Fisher cross-check:

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

## 5. Monitor progress

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

## 6. Resume after failure

Use the same approach and strata with `--resume`:

```bash
python scripts/run_analysis.py \
  --approach 2 \
  --strata cd4_thymus,cd8_thymus \
  --resume
```

The runner restores outputs from completed cells, re-executes bootstrap cells to rebuild Python state, and skips completed analytical cells. Approach 2 reuses its cached coordinate basis only when the stored clonotype fingerprint matches the current input.

## 7. Locate results

| Result | Location |
|---|---|
| Statistical tables | `outputs/tables/<approach>/` |
| Per-stratum conclusions | `outputs/tables/<approach>/<stratum>_conclusion.md` |
| Executed notebooks | `outputs/notebooks/` |
| Logs | `outputs/logs/` |
| Cached heavy objects | `outputs/cache/` |
| Figures | `figures/<approach>/<stratum>/` |

Every matplotlib figure is saved as PNG and PDF. The validation script rejects active plotting code that has no persistent figure output.

## 8. Interpret the cross-approach output

Approach-specific V-segment ranks are compared within the same biological stratum. The comparison reports:

1. pairwise rank correlation;
2. top-10 overlap;
3. signed effect correlation on a common g1-versus-allogeneic axis;
4. directional agreement;
5. direction-aware consensus rank.

A segment can therefore receive repeated ranking support while retaining a mixed effect direction. The pipeline reports that disagreement explicitly.

## 9. Control CPU usage

The runtime respects Linux CPU affinity and Slurm allocation variables. To reduce CPU use within an allocation:

```bash
ARTICLE_N_JOBS=8 python scripts/run_analysis.py --approach 2
```

Do not set `ARTICLE_N_JOBS` above the scheduler allocation. See [docs/parallel_execution.md](docs/parallel_execution.md) for monitoring commands and memory policy.

## 10. Final validation

After code or notebook changes:

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

A full scientific validation additionally requires execution with the canonical repertoire table. Regenerated tables, figures, and conclusions should be reviewed together before manuscript use.
