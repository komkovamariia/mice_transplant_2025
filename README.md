# Mouse T-cell receptor repertoire analysis after transplantation

[![Tests](https://github.com/komkovamariia/mice_transplant_2025/actions/workflows/validation.yml/badge.svg)](https://github.com/komkovamariia/mice_transplant_2025/actions/workflows/validation.yml)

**Authors:** Komkova M., Andreev V., Chernov P., Kofiadi I.

This repository contains the analysis code for TRA and TRB repertoires in the mouse transplantation study. The main workflow combines exact aaV clonotype-set comparisons, V-segment retention, UMI abundance and differential-abundance analysis with edgeR. A separate spike-in experiment is used to measure how much added TRA signal is required for recovery by the same pipeline.

[Run guide](RUN_GUIDE.md) · [Methods](docs/methods.md) · [Spike-in analysis](docs/spike_in.md) · [Aldan-3 runs](docs/aldan3_parallel.md) · [Input provenance](docs/notebook_input_provenance.md) · [Change history](CHANGELOG.md)

## Quick start

```bash
git clone https://github.com/komkovamariia/mice_transplant_2025.git
cd mice_transplant_2025
conda env create -f environment.yml
conda activate mice-transplant-2025
python scripts/validate_repository.py
python scripts/run_analysis.py --approach 1
```

Approach 1 is the main analysis. It runs nine biological strata, starting with the combined CD4+ and CD8+ repertoires from thymus and spleen. To run only the combined analysis:

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined
```

The default input paths are:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
```

The clonoset index points to the MiXCR exports used by the notebook. Paths can be changed with `--metadata-csv`, `--clonoset-index` and `--working-dir`, or through [config/example.env](config/example.env). Approach 1 does not require the derived Parquet table.

## Biological strata

| Run order | Key | Samples | Executed notebook |
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

For pooled strata, UMI counts from repeated samples of the same mouse are summed before count-based inference. Each treatment-group/mouse pair remains one analysis unit.

The main source notebook is `venn_original.ipynb`. The filename is kept because it is used by the runner and by archived analyses.

## Spike-in analysis

A local spike-in run can be started with:

```bash
python scripts/run_analysis.py --mode spike-in --spike-run-id sensitivity_01
```

The experiment starts from a fresh `all_combined` baseline, fixes an eligible TRA family and then adds counts to g1 independently at each dose. The control groups are unchanged. The default dose grid spans 0.00001% to 1% added family UMI relative to the original g1 library size. Target selection uses a fixed seed.

The analysis records clonotype recovery, TRAV retention, TRAV share in g1, edgeR log2 fold change and FDR, and the final TRAV rank. Full selection and dose definitions are in [docs/spike_in.md](docs/spike_in.md).

For Aldan-3, use `scripts/slurm_spike.py`; the current cluster setup is described in [docs/aldan3_parallel.md](docs/aldan3_parallel.md).

## Results and figures

| Output | Location |
|---|---|
| Main analysis source | `venn_original.ipynb` |
| Executed notebooks | `audit_runs/` |
| Logs | `logs/` |
| Figures | `figures/01_set_count/<stratum>/` |
| Figure archive | `figures/01_set_count.zip` |
| Numerical results | `results/01_set_count/<stratum>/` |
| Cross-stratum TRAV summary | `results/01_set_count/all_strata_distinctive_trav_summary.csv` |
| Spike-in results | `results/01_spike_in/<run_id>/` |
| Spike-in figures | `figures/01_spike_in/<run_id>/` |

Figures from the primary analysis are saved as PNG and PDF. The retained TRA and TRB heatmaps use the same candidate ordering as the article analysis. Formulas and score definitions are documented in [docs/methods.md](docs/methods.md).

Progress can be followed directly from the logs:

```bash
tail -f logs/00_all_combined.log
tail -f logs/01_cd4_thymus.log
```

To rerun a standard analysis after interruption:

```bash
python scripts/run_analysis.py --approach 1 --resume
```

## Additional analyses

| Command | Analysis |
|---|---|
| `--approach 2` | Sequence-space enrichment |
| `--approach 3` | Clone-level alloreactivity |
| `--approach 4` | Cross-method comparison |

Approaches 2 and 3 use the sample-resolved `clean_clonotypes_aaV.parquet` input described in [docs/notebook_input_provenance.md](docs/notebook_input_provenance.md).

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

The GitHub Actions workflow runs the same repository checks and test suite, plus an R-enabled edgeR integration test. Reproducing the biological results still requires the study data and the corresponding MiXCR exports.

[CITATION.cff](CITATION.cff) contains the software citation metadata and the study author list.
