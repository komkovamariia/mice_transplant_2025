# Approach 2: sequence embedding and repertoire geometry

This approach asks whether biological groups occupy different regions of TCR sequence space and whether the shift is detectable at both local-clonotype and whole-repertoire scales.

A single TCREmp/PCA coordinate basis is fitted across the canonical repertoire. Each biological stratum is then tested independently with exact k-d-tree local-density enrichment, RFF-MMD repertoire distances, deterministic 9,999-permutation PERMANOVA, witness scoring, and a compact motif-convergence summary. The common coordinate basis is a technical reference only; all biological tests remain stratum-specific.

Run from the repository root:

```bash
python scripts/run_analysis.py --approach 2
```

Tables are written to `outputs/tables/02_sequence_embedding/`; figures are written to `figures/02_sequence_embedding/<stratum>/`. The global embedding basis is cached under `outputs/cache/02_sequence_embedding/`.
