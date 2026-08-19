# Approach 1: exact sets and differential counts

This approach asks whether exact aaV clonotypes and their V-segment representation differ between g1 and the allogeneic g5+g6 reference.

For each of the six biological strata, the notebook performs three linked but distinct calculations: exact g1/g5/g6 set overlap, edgeR quasi-likelihood differential abundance, and Fisher's exact test as a count-based cross-check. V segments are ranked primarily by the number of edgeR-significant aaV clonotypes, with effect size and g1-exclusive support retained as secondary evidence.

Run from the repository root:

```bash
python scripts/run_analysis.py --approach 1
```

Tables are written to `outputs/tables/01_set_count/`; figures are written to `figures/01_set_count/<stratum>/`.
