# Results

Generated numerical outputs are excluded from Git.

Approach 1 writes to `01_set_count/<stratum>/`. A full nine-stratum run also writes `01_set_count/all_strata_distinctive_trav_summary.csv`.

Spike-in runs write to `01_spike_in/<run_id>/`, including the target selection, matrix snapshots, per-dose metrics and sensitivity summary.

The optional approaches use `02_sequence_embedding/`, `03_clone_alloreactivity/` and `04_cross_approach/`. Shared caches may appear under `cache/`.

Executed notebooks are stored in `audit_runs/`, logs in `logs/`, and figures in `figures/`.
