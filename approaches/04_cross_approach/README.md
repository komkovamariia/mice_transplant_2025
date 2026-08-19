# Cross-approach comparison

This stage compares the three upstream V-segment rankings without recomputing their underlying statistics.

For each biological stratum, it reports pairwise Spearman rank correlation, top-10 overlap, and a consensus ranking based on normalized within-approach ranks. The final cross-stratum table records how often a V segment appears among the top consensus results.

Run after Approaches 1–3:

```bash
python scripts/run_analysis.py --approach 4
```

or run the complete pipeline:

```bash
python scripts/run_analysis.py --approach all
```

Tables are written to `outputs/tables/04_cross_approach/`; figures are written to `figures/04_cross_approach/`.
