# Mouse TCR repertoires after transplantation

[![Repository validation](https://github.com/komkovamariia/mice_transplant_2025/actions/workflows/validation.yml/badge.svg)](https://github.com/komkovamariia/mice_transplant_2025/actions/workflows/validation.yml)

Analysis notebooks for TRA/TRB repertoires in the mouse transplantation study. The
primary workflow measures exact clonotype overlap, V-segment preservation and
clonotype abundance. It compares g1 with g5 and g6 using Fisher tests and edgeR.

The default analysis begins with **CD4 + CD8 T cells from thymus + spleen** and then
runs the eight compartment-specific analyses. Each run produces an executed notebook,
figures, numerical tables and a record of its software environment. An optional
**spike-in experiment** tests recovery of a fixed family of observed TRA clonotypes
across increasing UMI doses.

[Run guide](RUN_GUIDE.md) · [Methods](docs/methods.md) ·
[Spike-in protocol](docs/spike_in.md) · [Aldan-3 parallel runs](docs/aldan3_parallel.md) · [Changes](CHANGELOG.md) ·
[Input provenance](docs/notebook_input_provenance.md)

## Installation and first run

Use the study's Linux/HPC environment and existing MiXCR exports:

```bash
git clone https://github.com/komkovamariia/mice_transplant_2025.git
cd mice_transplant_2025
conda env create -f environment.yml
conda activate mice-transplant-2025
python scripts/validate_repository.py
python scripts/run_analysis.py --approach 1
```

`python scripts/run_analysis.py` also selects Approach 1. The nine analyses run in
sequence, beginning with `all_combined`. To run only the combined analysis:

```bash
python scripts/run_analysis.py --approach 1 --strata all_combined
```

Approach 1 reads these existing study files by default:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
```

The clonoset index must reference accessible MiXCR exports. Override paths with
`--metadata-csv`, `--clonoset-index` and `--working-dir`, or the environment variables
in [config/example.env](config/example.env). Raw FASTQ files and derived Parquet tables
are not required for Approach 1.

## Biological strata

The combined stratum pools the sampled CD4 and CD8 subsets. It does not identify a
CD4/CD8 double-positive population. Combined analyses sum UMI counts within each
mouse for inference, so tissues from the same mouse remain one biological replicate.

| Run order | Key | Cell subsets and organs | Executed notebook |
|---|---|---|---|
| 1 | `all_combined` | CD4 + CD8, thymus + spleen | `audit_runs/00_all_combined.executed.ipynb` |
| 2 | `cd4_thymus` | CD4, thymus | `audit_runs/01_cd4_thymus.executed.ipynb` |
| 3 | `cd8_thymus` | CD8, thymus | `audit_runs/02_cd8_thymus.executed.ipynb` |
| 4 | `cd4_spleen` | CD4, spleen | `audit_runs/03_cd4_spleen.executed.ipynb` |
| 5 | `cd8_spleen` | CD8, spleen | `audit_runs/04_cd8_spleen.executed.ipynb` |
| 6 | `cd4_combined` | CD4, thymus + spleen | `audit_runs/05_cd4_thymus_spleen.executed.ipynb` |
| 7 | `cd8_combined` | CD8, thymus + spleen | `audit_runs/06_cd8_thymus_spleen.executed.ipynb` |
| 8 | `thymus_combined` | CD4 + CD8, thymus | `audit_runs/07_thymus_cd4_cd8.executed.ipynb` |
| 9 | `spleen_combined` | CD4 + CD8, spleen | `audit_runs/08_spleen_cd4_cd8.executed.ipynb` |

The original eight filenames remain stable. `venn_original.ipynb` is the executable
source for every Approach 1 run. Sample selection precedes count-matrix construction.

## Spike-in sensitivity experiment

```bash
python scripts/run_analysis.py --mode spike-in --spike-run-id sensitivity_01
```

This runs a fresh combined baseline, chooses a noncandidate TRAV with ten eligible
observed aaV clonotypes, and adds the same family to g1 at **0.01%, 0.1% and 1%** of
each sample's original UMI library. Each dose starts from frozen baseline matrices.
Control counts remain unchanged. A fixed seed makes target selection reproducible.

The protocol checks recovery in the g1 remainder, V-segment preservation and share,
positive edgeR log2 fold change with FDR < 0.05, and movement in the final TRAV ranking.
It records hypotheses that fail as well as those that pass. UMI share and unique-aaV
share are reported separately. See [the protocol](docs/spike_in.md) for target
eligibility, rounding, custom doses and interpretation of sensitivity bounds.

## Results and figures

| Artifact | Location |
|---|---|
| Primary analysis source | `venn_original.ipynb` |
| Executed notebooks and logs | `audit_runs/`, `logs/` |
| Standard figures, PNG and PDF | `figures/01_set_count/<stratum>/` |
| Complete standard figure archive | `figures/01_set_count.zip` |
| Standard numerical results | `results/01_set_count/<stratum>/` |
| Nine-stratum TRAV summary | `results/01_set_count/all_strata_distinctive_trav_summary.csv` |
| Spike-in notebooks | `audit_runs/spike_in/<run_id>/<case>.executed.ipynb` |
| Spike-in results and selection | `results/01_spike_in/<run_id>/` |
| Spike-in figures and archive | `figures/01_spike_in/<run_id>/`, `<run_id>.zip` |

Every displayed Matplotlib figure is saved in PNG and PDF. There are two V-segment
heatmaps per standard stratum or spike-in case: TRA and TRB. They use the ordered
bubble-plot candidates selected by the article preservation score, with four
regions containing g1. The [methods](docs/methods.md) specify the formulas and thresholds.

A standard rerun replaces generated figures and tables for the selected strata.
Spike-in runs require a fresh identifier and preserve ordinary results. The runner
checks figure manifests and successful edgeR fitting before accepting completion.
A failed run retains its executed notebook and cell log for diagnosis.

```bash
tail -f logs/00_all_combined.log
# For the original first compartment:
tail -f logs/01_cd4_thymus.log
# Replay a standard run after interruption:
python scripts/run_analysis.py --approach 1 --resume
```

## Additional approaches

Approaches 2 and 3 require an existing, sample-resolved
`clean_clonotypes_aaV.parquet`. Its required columns and the remaining preparation
boundary are described in [input provenance](docs/notebook_input_provenance.md).

```bash
python scripts/run_analysis.py --approach 2 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 3 --data-dir /absolute/path/to/derived_data
python scripts/run_analysis.py --approach 4
python scripts/run_analysis.py --approach all --data-dir /absolute/path/to/derived_data
```

Approach 2 examines sequence-space enrichment; Approach 3 integrates clone-level
evidence; Approach 4 compares their TRAV rankings with Approach 1. All expose the same
nine strata. They remain optional when reproducing the primary notebook analysis.

## Reproducibility and availability

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

CI checks the source notebooks, pooling, injected-count propagation, target selection,
rounding and output contracts. A separate R-enabled job fits the notebook's actual
edgeR cell on synthetic counts. These checks do not establish the study's biological
results: numerical reproduction requires the study data and review of executed notebooks.

The repository provides analysis code and environment specifications. Study MiXCR
exports and the derived Parquet table are supplied separately; no public data accession
or manuscript DOI is recorded here. Runtime results include source hashes, package
versions and R session information. Cite the repository revision used for an analysis;
[CITATION.cff](CITATION.cff) provides software citation metadata.

Use [GitHub issues](https://github.com/komkovamariia/mice_transplant_2025/issues) for
reproducible problems and [CONTRIBUTING.md](CONTRIBUTING.md) for change validation.
