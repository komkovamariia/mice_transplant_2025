# Approach 2: sequence-space enrichment and repertoire geometry

This approach tests the g1 versus g5+g6 contrast at local sequence-space and whole-repertoire scales.

One TCREmp/PCA basis is fitted as a shared technical reference. Each biological stratum is then analyzed independently with exact k-d-tree density enrichment, RFF-MMD, deterministic two-group PERMANOVA, witness scoring, and motif summaries. Combined strata pool tissues within mouse before repertoire embedding.

Run:

```bash
python scripts/run_analysis.py --approach 2
```

Tables and conclusions are written to `outputs/tables/02_sequence_embedding/`. Figures are written to `figures/02_sequence_embedding/<stratum>/`. The validated embedding cache is stored under `outputs/cache/02_sequence_embedding/`.

