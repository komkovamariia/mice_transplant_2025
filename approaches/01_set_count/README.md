# Approach 1: exact sets and differential counts

This approach compares g1 with the allogeneic g5+g6 reference using exact aaV clonotypes.

For every biological stratum it reports exact g1/g5/g6 set overlap, edgeR quasi-likelihood differential abundance, Fisher's exact-test support, and a V-segment ranking. In the combined CD4 and CD8 strata, thymus and spleen counts are pooled within each mouse before model fitting.

Run:

```bash
python scripts/run_analysis.py --approach 1
```

Tables and result-specific conclusions are written to `outputs/tables/01_set_count/`. Figures are written to `figures/01_set_count/<stratum>/`.

