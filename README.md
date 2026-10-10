# Mouse T-cell receptor repertoire analysis after transplantation

[![Tests](https://github.com/komkovamariia/mice_transplant_2025/actions/workflows/validation.yml/badge.svg)](https://github.com/komkovamariia/mice_transplant_2025/actions/workflows/validation.yml)

Authors: Komkova M., Andreev V., Chernov P., Kofiadi I.

This repository contains the TRA and TRB repertoire analysis used in the mouse transplantation study. Approach 1 is the main workflow. It combines exact aaV clonotype sets, V-segment retention, UMI abundance, Fisher testing and edgeR. A separate spike-in experiment measures how an added TRA family propagates through the same pipeline.

Documentation: [run guide](RUN_GUIDE.md), [methods](docs/methods.md), [spike-in protocol](docs/spike_in.md), [Aldan-3 guide](docs/aldan3_parallel.md), [input provenance](docs/notebook_input_provenance.md).

## Quick start

```bash
git clone https://github.com/komkovamariia/mice_transplant_2025.git
cd mice_transplant_2025
conda env create -f environment.yml
conda activate mice-transplant-2025
python scripts/validate_repository.py
python scripts/run_analysis.py --approach 1
```

Approach 1 reads the study metadata, the clonoset index and the MiXCR exports referenced by that index. Its default paths are:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
```

Override them with `--metadata-csv`, `--clonoset-index` and `--working-dir`, or set the matching variables from [config/example.env](config/example.env).

To run only the combined repertoire:

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined
```

To run the combined, pooled CD4 and pooled CD8 strata:

```bash
python scripts/run_analysis.py --approach 1 \
  --strata all_combined,cd4_combined,cd8_combined
```

## Approach 1 strata

| Order | Key | Samples |
|---:|---|---|
| 1 | `all_combined` | CD4+ + CD8+, thymus + spleen |
| 2 | `cd4_thymus` | CD4+, thymus |
| 3 | `cd8_thymus` | CD8+, thymus |
| 4 | `cd4_spleen` | CD4+, spleen |
| 5 | `cd8_spleen` | CD8+, spleen |
| 6 | `cd4_combined` | CD4+, thymus + spleen |
| 7 | `cd8_combined` | CD8+, thymus + spleen |
| 8 | `thymus_combined` | CD4+ + CD8+, thymus |
| 9 | `spleen_combined` | CD4+ + CD8+, spleen |

The source notebook is `venn_original.ipynb`. The name is kept for compatibility with the runner and archived analyses.

Structural set analysis uses functional clonotypes and downsamples each eligible library to 15,000 UMI. `all_combined`, `cd4_combined` and `cd8_combined` use 100 baseline rarefactions to choose the representative structural seed. Fisher and edgeR use the corresponding functional raw UMI counts. Pooled strata sum repeated samples within treatment group and mouse before count-based inference. See [docs/methods.md](docs/methods.md) for the exact definitions.

## Spike-in experiment

A local run starts with:

```bash
python scripts/run_analysis.py --mode spike-in --spike-run-id sensitivity_01
```

The target TRA family is selected from the unmodified `all_combined` baseline and then reused at every dose. Controls are unchanged. The experiment records aaV recovery, TRAV retention and abundance, edgeR logFC/FDR and final TRAV rank. Dose definitions and interpretation are in [docs/spike_in.md](docs/spike_in.md).

For Slurm runs on Aldan-3, use `scripts/slurm_spike.py` with [docs/aldan3_parallel.md](docs/aldan3_parallel.md).

## Outputs

| Output | Location |
|---|---|
| Executed notebooks | `audit_runs/` |
| Logs | `logs/` |
| Approach 1 figures | `figures/01_set_count/<stratum>/<TRA|TRB>/` |
| Approach 1 tables | `results/01_set_count/<stratum>/` |
| Approach 1 figure archive | `figures/01_set_count.zip` |
| Spike-in results | `results/01_spike_in/<run_id>/` |
| Spike-in figures | `figures/01_spike_in/<run_id>/` |

Generated results, figures, logs and executed notebooks are not tracked.

## Other approaches

```bash
python scripts/run_analysis.py --approach 2 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 3 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 4
python scripts/run_analysis.py --approach all --data-dir /absolute/path/to/derived_data
```

Approaches 2 and 3 use the derived `clean_clonotypes_aaV.parquet` interface described in [docs/notebook_input_provenance.md](docs/notebook_input_provenance.md). Approach 4 reads the standardized outputs of the preceding approaches.

## Checks

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

GitHub Actions runs the repository checks, the Python test suite and an R-enabled edgeR integration test. Reproducing the study-level numbers still requires the study data and MiXCR exports.

Software citation metadata are in [CITATION.cff](CITATION.cff).
