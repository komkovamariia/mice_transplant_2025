# Approach 3: clone-level alloreactivity prioritization

This approach asks which individual aaV clonotypes receive repeated, sample-resolved support from signals that are independent of the embedding-density test.

For each biological stratum, clonotypes are scored using allogeneic mouse prevalence, allogeneic-to-g1 abundance shift, and same-V CDR3 convergence. CDR3 length, hydrophobicity, and approximate net charge are retained for interpretation. The integrated score is a prioritization statistic and must not be interpreted as a calibrated probability of antigen specificity or alloreactivity.

Run from the repository root:

```bash
python scripts/run_analysis.py --approach 3
```

Tables are written to `outputs/tables/03_clone_alloreactivity/`; figures are written to `figures/03_clone_alloreactivity/<stratum>/`.
