# Approach 3: clone-level evidence prioritization

This approach ranks individual aaV clonotypes using three mouse-resolved signals: allogeneic prevalence, within-mouse relative-abundance shift, and same-V CDR3 convergence.

Each mouse receives equal abundance weight. CDR3 length, hydrophobicity, and approximate charge are reported for interpretation only. The integrated score is an evidence ranking and is not a calibrated probability of alloreactivity or antigen specificity.

Run:

```bash
python scripts/run_analysis.py --approach 3
```

Tables and conclusions are written to `outputs/tables/03_clone_alloreactivity/`. Figures are written to `figures/03_clone_alloreactivity/<stratum>/`.

