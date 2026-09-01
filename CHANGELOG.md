# Changelog

## 2026-09-01: Bounded Aldan-3 parallel sensitivity runs

- Added preparation, read-only cluster diagnostics and explicit Slurm submission for the combined spike-in experiment.
- Added a 24-CPU experiment budget: baseline up to 24 workers, then at most three independent 8-CPU dose notebooks.
- Added sixteen scheduler-default family doses spanning 0.00001% to 1%.
- Separated baseline, dose and finalization stages with success dependencies, isolated logs and exclusive case claims.
- Added source/configuration/selection checks and retained immutable baseline matrix checksums.
- Prevented explicit worker overrides from exceeding Slurm or affinity limits; defaulted unscheduled execution to one worker.
- Parallelized historical two-sided Fisher tests in ordered batches, retaining one global BH correction.
- Preserved the per-notebook edgeR model, mouse-level pooling, g1/control identities and all biological definitions.
- Added scheduling and numerical-equivalence regression tests; real-data timing and peak memory require HPC measurement.

## 2026-09-01: Combined analysis and spike-in sensitivity

### Analysis

- Added `all_combined` for CD4 + CD8, thymus + spleen and made it the first default stratum.
- Retained all eight existing strata and their executed-notebook filenames; the new file is `00_all_combined.executed.ipynb`.
- Made Approach 1 the CLI default. The full sequence remains available through `--approach all`.
- Reused group count matrices for both exact-set analysis and inference so injected counts cannot be lost in a second input read.
- Preserved mouse-level pooling, historical exclusions, exact aaV identity, Fisher, edgeR and the existing final TRAV ranking.
- Corrected the heatmap legend to identify its UMI-weighted regional proportions.
- Fixed deferred lambda expressions in the Approach 2 and 3 notebooks so each selected stratum actually executes, and added the combined stratum to all optional approaches.

### Spike-in mode

- Added a separate baseline-plus-dose experiment through the full primary notebook.
- Added deterministic selection of an observed aaV family sharing a baseline noncandidate TRAV.
- Required absence from subtraction controls and a declared maximum baseline frequency in g1.
- Added family doses of 0.01%, 0.1% and 1% of each original g1 sample library, configurable through CLI options.
- Added integer allocation, realized-dose reporting, immutable matrix snapshots and checksum validation.
- Added clonotype recovery, retention, unique-aaV and UMI shares, edgeR filtering/FDR results and final-rank changes.
- Added observed sensitivity brackets and explicit reporting when a limit is outside the tested range or not established.
- Isolated experimental figures, ZIP archives, executed notebooks, logs and tables from ordinary article outputs.

### Reproducibility and presentation

- Added source hashes, package versions, analysis-unit tables and R session information.
- Updated the README, run guide, method descriptions, data-availability boundaries and software citation metadata.
- Added regression coverage and an R-enabled CI test of the original edgeR notebook cell.
- Renamed the combined result table to `all_strata_distinctive_trav_summary.csv`.

## Earlier notebook restoration

- Restored the primary notebook workflow from historical commit `6678766` and removed DESeq2 from active analysis.
- Corrected historical metadata keys before chain normalization.
- Added independent executed notebooks, cell logs and PNG/PDF persistence for the original eight strata.
- Aligned bubble labels with the article preservation score and restricted their companion heatmaps to four g1-containing regions.
- Removed g2-centered duplicate figures and auxiliary pairwise/mouse-overlap figure panels from article outputs.
- Added output validation and a complete `figures/01_set_count.zip` archive.
