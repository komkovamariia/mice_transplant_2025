# Figures

All analysis figures are stored beneath this directory. The primary notebook saves
every retained figure as PNG and PDF; each case has a verified figure manifest.

| Analysis | Figure directory | Archive |
|---|---|---|
| Standard Approach 1 | `01_set_count/<stratum>/` | `01_set_count.zip` |
| Spike-in experiment | `01_spike_in/<run_id>/<case>/` | `01_spike_in/<run_id>.zip` |
| Sequence embedding | `02_sequence_embedding/<stratum>/` | Not generated automatically |
| Clone alloreactivity | `03_clone_alloreactivity/<stratum>/` | Not generated automatically |
| Method comparison | `04_cross_approach/<stratum>/` | Not generated automatically |

Each Approach 1 stratum or experimental case has one TRA and one TRB V-segment
heatmap, using the genes selected by its bubble plot. Experimental summary curves
are stored at the experiment's figure root. Generated outputs are excluded from Git.
