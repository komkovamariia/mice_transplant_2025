# Reproducible analysis of mouse T-cell receptor repertoires after transplantation

[![Repository validation](https://github.com/komkovamariia/mice_transplant_2025/actions/workflows/validation.yml/badge.svg)](https://github.com/komkovamariia/mice_transplant_2025/actions/workflows/validation.yml)

**Authors:** Komkova M., Andreev V., Chernov P., Kofiadi I.

This repository contains the reproducible analysis workflow for TRA and TRB repertoires
in the mouse transplantation study. The primary analysis combines exact clonotype-set
comparisons, V-segment retention, clonotype abundance and differential-abundance
inference. A count-matrix spike-in experiment provides an explicit sensitivity analysis
for recovery of a known TRA family through the same workflow.

The default analysis begins with the combined **CD4+ and CD8+ T-cell repertoires from
thymus and spleen** and then runs eight compartment-specific analyses. Each run records
its executed notebook, numerical results, publication-ready figures and software
provenance.

[Run guide](RUN_GUIDE.md) · [Methods](docs/methods.md) ·
[Spike-in sensitivity analysis](docs/spike_in.md) ·
[Aldan-3 parallel execution](docs/aldan3_parallel.md) ·
[Input provenance](docs/notebook_input_provenance.md) ·
[Repository writing and naming style](docs/repository_style.md) ·
[Change history](CHANGELOG.md)

## Installation and first analysis

Use the study Linux/HPC environment and existing MiXCR exports:

```bash
git clone https://github.com/komkovamariia/mice_transplant_2025.git
cd mice_transplant_2025
conda env create -f environment.yml
conda activate mice-transplant-2025
python scripts/validate_repository.py
python scripts/run_analysis.py --approach 1
```

`python scripts/run_analysis.py` also selects the primary exact-set and count-based
analysis. The nine biological strata run in sequence, beginning with `all_combined`.
To run only the combined analysis:

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined
```

The primary analysis reads these existing study files by default:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
```

The clonoset index must reference accessible MiXCR exports. Override paths with
`--metadata-csv`, `--clonoset-index` and `--working-dir`, or with the environment
variables in [config/example.env](config/example.env). Raw FASTQ files and derived
Parquet tables are not required for the primary analysis.

## Biological strata

The combined stratum pools the sampled CD4+ and CD8+ subsets. It does not represent a
CD4/CD8 double-positive cell population. Combined count-based analyses sum UMI counts
within each mouse so that tissues from the same mouse remain one biological replicate.

| Run order | Key | Cell subsets and organs | Executed notebook |
|---|---|---|---|
| 1 | `all_combined` | CD4+ + CD8+, thymus + spleen | `audit_runs/00_all_combined.executed.ipynb` |
| 2 | `cd4_thymus` | CD4+, thymus | `audit_runs/01_cd4_thymus.executed.ipynb` |
| 3 | `cd8_thymus` | CD8+, thymus | `audit_runs/02_cd8_thymus.executed.ipynb` |
| 4 | `cd4_spleen` | CD4+, spleen | `audit_runs/03_cd4_spleen.executed.ipynb` |
| 5 | `cd8_spleen` | CD8+, spleen | `audit_runs/04_cd8_spleen.executed.ipynb` |
| 6 | `cd4_combined` | CD4+, thymus + spleen | `audit_runs/05_cd4_thymus_spleen.executed.ipynb` |
| 7 | `cd8_combined` | CD8+, thymus + spleen | `audit_runs/06_cd8_thymus_spleen.executed.ipynb` |
| 8 | `thymus_combined` | CD4+ + CD8+, thymus | `audit_runs/07_thymus_cd4_cd8.executed.ipynb` |
| 9 | `spleen_combined` | CD4+ + CD8+, spleen | `audit_runs/08_spleen_cd4_cd8.executed.ipynb` |

`venn_original.ipynb` remains the stable executable source for every run of the primary
analysis. The historical filename is retained because it is part of the reproducibility
contract. In documentation and article-oriented descriptions it is referred to as the
**primary exact-set and count-based analysis notebook**.

## Spike-in sensitivity analysis

A local sequential run can be started with:

```bash
python scripts/run_analysis.py --mode spike-in --spike-run-id sensitivity_01
```

The experiment first creates a fresh combined baseline, selects a noncandidate TRAV
family containing ten eligible observed aaV clonotypes and freezes that biological
identity before any perturbation. Each dose is generated independently from the frozen
baseline count matrices. Only g1 counts are modified; subtraction controls remain
unchanged.

The current default dose grid is dense and approximately logarithmic, from **0.00001% to
1% added family UMI relative to the original g1 sample library**. Custom fractions can
be supplied explicitly, including larger doses when required for ranking-sensitivity
experiments. Target selection is reproducible with a fixed seed.

The protocol evaluates five linked outcomes: recovery of the introduced aaV clonotypes
in the g1 remainder, increased TRAV remainder and retention, increased TRAV share in g1,
positive edgeR log2 fold change with FDR < 0.05 for introduced aaV, and improvement of
the selected TRAV in the final integrated ranking. See
[Spike-in sensitivity analysis](docs/spike_in.md) for eligibility rules, rounding,
parallel execution and interpretation limits.

For bounded Aldan-3 execution, use `scripts/slurm_spike.py`. A prepared run follows the
explicit dependency graph **baseline -> dose array -> final summary**, with the dose
array throttled to the configured CPU budget.

## Results and figures

| Scientific output | Location |
|---|---|
| Primary analysis source | `venn_original.ipynb` |
| Executed notebooks and cell logs | `audit_runs/`, `logs/` |
| Standard publication figures | `figures/01_set_count/<stratum>/` |
| Complete standard figure archive | `figures/01_set_count.zip` |
| Standard numerical results | `results/01_set_count/<stratum>/` |
| Cross-stratum TRAV summary | `results/01_set_count/all_strata_distinctive_trav_summary.csv` |
| Spike-in executed notebooks | `audit_runs/spike_in/<run_id>/<case>.executed.ipynb` |
| Spike-in numerical results and target selection | `results/01_spike_in/<run_id>/` |
| Spike-in sensitivity figures and archive | `figures/01_spike_in/<run_id>/`, `<run_id>.zip` |

Every displayed Matplotlib figure is persisted in PNG and PDF. Standard runs include
TRA and TRB V-segment heatmaps built from the ordered candidates used by the article
retention score. The formulas and biological interpretation are defined in
[Methods](docs/methods.md).

A standard rerun replaces generated figures and tables for the selected strata.
Spike-in runs require a fresh identifier and preserve ordinary analysis results. The
runner verifies required numerical outputs, figure manifests and successful edgeR
fitting before accepting a run as complete. A failed notebook execution retains its
checkpoint and cell log for diagnosis.

```bash
tail -f logs/00_all_combined.log
# original first compartment
tail -f logs/01_cd4_thymus.log
# replay a standard run after interruption
python scripts/run_analysis.py --approach 1 --resume
```

## Additional analyses

The optional analyses retain their stable command-line identifiers while using
article-oriented names in documentation:

| Command identifier | Scientific description |
|---|---|
| `--approach 2` | Sequence-space enrichment analysis |
| `--approach 3` | Clone-level alloreactivity analysis |
| `--approach 4` | Cross-method evidence comparison |

The sequence-space and clone-level analyses require an existing, sample-resolved
`clean_clonotypes_aaV.parquet`. Its required columns and preparation boundary are
described in [Input provenance](docs/notebook_input_provenance.md).

```bash
python scripts/run_analysis.py --approach 2 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 3 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 4
python scripts/run_analysis.py --approach all --data-dir /absolute/path/to/derived_data
```

All analyses expose the same nine biological strata. The additional analyses remain
optional when reproducing the primary notebook results.

## Reproducibility and validation

```bash
python scripts/normalize_comments.py
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

Continuous integration checks source notebooks, mouse-level pooling, spike-in count
propagation, target selection, dose rounding and output contracts. A separate R-enabled
job fits the notebook's edgeR cell on synthetic counts. These checks verify software
behavior; numerical reproduction of the biological results still requires the study
data and review of the executed notebooks.

The repository provides analysis code and environment specifications. Study MiXCR
exports and the derived Parquet table are supplied separately. Runtime results record
source hashes, package versions and R session information. Cite the exact repository
revision used for an analysis; [CITATION.cff](CITATION.cff) contains the software
citation metadata and the study author list.

Use [GitHub issues](https://github.com/komkovamariia/mice_transplant_2025/issues) for
reproducible problems and [CONTRIBUTING.md](CONTRIBUTING.md) for change validation.
