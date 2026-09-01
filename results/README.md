# Numerical results

Standard Approach 1 writes tables, conclusions, analysis-unit metadata, software
provenance and figure manifests to `01_set_count/<stratum>/`. A complete run also
writes `01_set_count/all_strata_distinctive_trav_summary.csv`.

Spike-in experiments write to `01_spike_in/<run_id>/`. This includes the fixed target
selection, hashed count snapshots, the full outputs of baseline and dose notebooks,
per-clonotype recovery and a sensitivity summary. Ordinary results remain separate.

The optional approaches use `02_sequence_embedding/`, `03_clone_alloreactivity/` and
`04_cross_approach/`. Shared runtime caches may appear under `cache/`.

Executed notebooks are under `audit_runs/`, logs under `logs/`, and all figures under
`figures/` at the repository root. Generated results are excluded from Git.
