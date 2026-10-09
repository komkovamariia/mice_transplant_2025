# Figures

All analysis figures are stored beneath this directory. The primary notebook saves
every retained figure as PNG; each case has a verified figure manifest.

| Analysis | Figure directory | Archive |
|---|---|---|
| Standard Approach 1 | `01_set_count/<stratum>/<TRA|TRB>/` | `01_set_count.zip` |
| Spike-in experiment | `01_spike_in/<run_id>/<case>/<TRA|TRB>/` | `01_spike_in/<run_id>.zip` |
| Sequence embedding | `02_sequence_embedding/<stratum>/TRA/` | Not generated automatically |
| Clone alloreactivity | `03_clone_alloreactivity/<stratum>/TRA/` | Not generated automatically |
| Method comparison | `04_cross_approach/<stratum>/TRA/` | Not generated automatically |

TRA and TRB publication figures are kept in separate chain directories. Each Approach 1 stratum has one TRA and one TRB V-segment heatmap using the genes selected by its bubble plot. Spike-in summary curves are stored as separate TRA PNG panels rather than a composite figure. Generated outputs are excluded from Git.
