# Mouse TCR repertoire analysis after allogeneic transplantation

This repository contains a reproducible, multi-approach analysis of mouse T-cell receptor repertoires in the BALB/c → C57BL/6 transplantation model. The computational design deliberately separates three analytical questions and reconciles them only after each method has produced an independent result.

## Analytical design

| Stage | Directory | Question | Primary output |
|---|---|---|---|
| 1 | `approaches/01_set_count/` | Which exact aaV clonotypes and V segments differ by set membership and abundance? | edgeR/Fisher clonotype statistics and V-segment ranks |
| 2 | `approaches/02_sequence_embedding/` | Which sequence-space neighborhoods and whole repertoires shift between g1 and allogeneic groups? | density enrichment, RFF-MMD, PERMANOVA, witness ranks |
| 3 | `approaches/03_clone_alloreactivity/` | Which individual clonotypes are repeatedly supported by independent allogeneic prevalence, abundance, and convergence signals? | clone-level candidate table and V-segment ranks |
| 4 | `approaches/04_cross_approach/` | Which V-segment signals are stable across methods and biological compartments? | per-stratum concordance and consensus ranks |

The analyses are performed independently for six prespecified biological strata:

`cd4_thymus`, `cd4_spleen`, `cd8_thymus`, `cd8_spleen`, `cd4_combined`, and `cd8_combined`.

The two combined strata contain thymus and spleen samples of the indicated T-cell subset. They do not pool CD4 and CD8 cells.

## Repository layout

```text
.
├── approaches/
│   ├── 01_set_count/
│   │   ├── README.md
│   │   └── set_count_analysis.ipynb
│   ├── 02_sequence_embedding/
│   │   ├── README.md
│   │   └── sequence_embedding_analysis.ipynb
│   ├── 03_clone_alloreactivity/
│   │   ├── README.md
│   │   └── clone_alloreactivity_analysis.ipynb
│   └── 04_cross_approach/
│       ├── README.md
│       └── cross_approach_comparison.ipynb
├── src/                         # Shared analysis code
├── scripts/
│   ├── run_analysis.py          # Cell-aware execution driver
│   └── validate_repository.py   # Static repository checks
├── figures/
│   ├── 01_set_count/
│   ├── 02_sequence_embedding/
│   ├── 03_clone_alloreactivity/
│   └── 04_cross_approach/
├── outputs/
│   ├── tables/
│   ├── notebooks/
│   ├── logs/
│   └── cache/
├── docs/
│   └── parallel_execution.md
├── supplementary/
├── environment.yml
├── requirements.txt
└── RUN_GUIDE_RU.md
```

All generated figures are written to the single `figures/` hierarchy. Every plotting routine writes both a 300-dpi PNG and a vector PDF. Notebook-only figures are not used as the sole copy of a result.

## Input data

The canonical input is a sample-resolved aaV clonotype table:

```text
data/clean_clonotypes_aaV.parquet
```

The minimum required columns are:

```text
cdr3, v_gene, umi, group, sample_id, mouse_id, source, subtype
```

`v_germ` and `treatment` are used when available. The code derives CD4/CD8 identity from `subtype` with `sample_id` as a fallback, and derives thymus/spleen identity from `source` with `sample_id` as a fallback. Samples that cannot be assigned unambiguously are not silently inserted into a stratum.

The input location can be changed without editing a notebook:

```bash
export MICE_TCR_DATA_DIR=/path/to/data
```

or:

```bash
export MICE_TCR_CLEAN_PARQUET=/path/to/clean_clonotypes_aaV.parquet
```

Primary sequencing data are not embedded in this repository. A clean execution therefore requires the canonical sample-resolved input table.

## Environment

Create the complete Python/R environment:

```bash
conda env create -f environment.yml
conda activate mice-transplant-2025
```

Update an existing environment:

```bash
conda env update -n mice-transplant-2025 -f environment.yml --prune
conda activate mice-transplant-2025
```

The count-based model uses edgeR through rpy2. DESeq2 is not part of the active pipeline. `repseq` is installed from a fixed Git commit, and the runtime uses the CPU allocation visible to the current process rather than the physical-node core count.

Validate the repository before a long run:

```bash
python scripts/validate_repository.py
```

## Running the analyses

The recommended interface is `scripts/run_analysis.py`. It executes notebooks cell by cell, writes a durable executed-notebook checkpoint after every successful code cell, and appends progress to `outputs/logs/`.

Run Approach 1 only:

```bash
python scripts/run_analysis.py --approach 1
```

Run Approach 2 only:

```bash
python scripts/run_analysis.py --approach 2
```

Run Approach 3 only:

```bash
python scripts/run_analysis.py --approach 3
```

Run the cross-approach comparison after all three upstream approaches:

```bash
python scripts/run_analysis.py --approach 4
```

Run the complete pipeline in publication order:

```bash
python scripts/run_analysis.py --approach all
```

Run more than one selected approach:

```bash
python scripts/run_analysis.py --approach 1,2
```

Run only selected biological strata:

```bash
python scripts/run_analysis.py \
  --approach 2 \
  --strata cd4_thymus,cd8_thymus
```

Specify the data directory on the command line:

```bash
python scripts/run_analysis.py \
  --approach all \
  --data-dir /path/to/data
```

### Resume after a failed long run

The execution driver checkpoints completed code-cell IDs and all upstream tables/figures. To resume:

```bash
python scripts/run_analysis.py --approach 2 --resume
```

Bootstrap cells are always re-executed to recreate Python state. Completed analytical cells are skipped only when they are present in the checkpoint. Approach 2 also caches the shared embedding basis under `outputs/cache/02_sequence_embedding/`, so a restart does not require refitting the global sequence-space basis when the canonical clonotype set is unchanged.

### Live progress

A log line is written when every code cell starts and completes, including the number of cells remaining. Long-running cells emit a heartbeat every 30 seconds. Follow the current analysis with:

```bash
tail -f outputs/logs/02_sequence_embedding.log
```

A typical sequence is:

```text
[CELL 003/009] START | 6 remaining after this cell
[CELL 003/009] RUNNING | elapsed=300s
[CELL 003/009] DONE in 412.7s | 6 code cell(s) not yet checkpointed
```

For scheduler-level CPU and memory monitoring, see `docs/parallel_execution.md`.

## Outputs

Each analytical approach writes stratum-specific tables to:

```text
outputs/tables/<approach>/
```

Each stratum produces a standardized V-segment ranking containing at least:

```text
v_gene, rank, score, stratum
```

Approach-specific evidence is retained in additional columns. This shared schema is the only input used by the cross-approach ranking layer.

Figures are written to:

```text
figures/<approach>/<stratum>/
```

and the global cross-stratum figure is written to:

```text
figures/04_cross_approach/
```

Executed notebooks and progress logs are written to:

```text
outputs/notebooks/
outputs/logs/
```

## Statistical interpretation

Approach 1 uses exact aaV clonotypes and an edgeR quasi-likelihood model for g1 versus the allogeneic g5+g6 reference; Fisher's exact test is retained as an independent count-based cross-check.

Approach 2 uses one shared TCREmp/PCA coordinate system, exact k-d-tree neighborhood enrichment, RFF-MMD repertoire distances, deterministic 9,999-permutation PERMANOVA, and witness scoring. The permutation schedule is generated before parallel execution, so the p-value does not depend on worker completion order.

Approach 3 is a prioritization analysis rather than an antigen-specificity classifier. Its integrated score combines allogeneic mouse prevalence, allogeneic-to-g1 abundance shift, and same-V local CDR3 convergence. Sequence physicochemical descriptors are reported for interpretation but are not treated as direct evidence of alloreactivity.

Approach 4 does not force agreement. It reports rank correlation, top-10 overlap, and a normalized consensus rank separately for each biological stratum. A robust biological signal is one that is reproduced across analytical lines and, ideally, across more than one compartment.

Historical PERMANOVA values from the earlier notebook implementation must not be reused: the original helper ignored permuted labels. The active implementation in `src/runtime.py` performs the intended deterministic permutation test, and only regenerated p-values should be reported.

## Reproducibility notes

The exact mirpy density backend remains `kdtree`; no approximate-nearest-neighbor substitution is used. CPU parallelism is scheduler-aware and does not alter clonotype definitions, biological contrasts, FDR thresholds, or permutation schedules.

The active repository contains only the current English computational workflow. Earlier notebook and manuscript variants remain available through Git history but are not part of the executable pipeline.
