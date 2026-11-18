# Figures

Generated figures are written below this directory and are excluded from Git.

| Analysis | Directory | Archive |
|---|---|---|
| Approach 1 | `01_set_count/<stratum>/<TRA|TRB>/` | `01_set_count.zip` |
| Spike-in | `01_spike_in/<run_id>/<case>/<TRA|TRB>/` | `01_spike_in/<run_id>.zip` |
| Sequence embedding | `02_sequence_embedding/<stratum>/TRA/` | none |
| Clone alloreactivity | `03_clone_alloreactivity/<stratum>/TRA/` | none |
| Cross-approach comparison | `04_cross_approach/<stratum>/TRA/` | none |

Approach 1 keeps TRA and TRB in separate directories. File names use the chain, group and analysis name, for example `TRA_g1_v_segment_representation.png` and `TRA_g1_v_segment_localization.png`.


The standard `all_combined` run writes nine manuscript PNGs to each chain directory:

| Reference | Filename after `TRA_g1_` or `TRB_g1_` |
|---|---|
| Figure 2 | `depth_check.png` |
| Figure 3 | `mds.png` |
| Figure 4 | `clonotype_overlap.png` |
| Figure 5 | `v_segment_representation.png` |
| Figure 6 | `v_segment_retention.png` |
| Figure 7 | `v_segment_localization.png` |
| Figure 8 | `edger_summary.png` |
| Figure 9 | `fisher_edger_concordance.png` |
| Figure 11 | `cd4_cd8_structure.png` |

Figure 10 is generated from the spike-in dose series in `01_spike_in/<run_id>/<TRA|TRB>/` as `<chain>_g1_spike_in_sensitivity.png`. Auxiliary plots, individual panel exports and per-dose plots are saved in the corresponding chain's `surplus/` directory.
