# Cross-approach comparison

This stage compares the three upstream V-segment rankings within each biological stratum.

It reports pairwise rank correlation, top-10 overlap, signed-effect correlation, directional agreement, and a direction-aware consensus rank. All effects use the same g1-versus-allogeneic orientation. Rank support and effect direction therefore remain distinct.

Run after Approaches 1–3:

```bash
python scripts/run_analysis.py --approach 4
```

Run the complete pipeline:

```bash
python scripts/run_analysis.py --approach all
```

Tables and conclusions are written to `outputs/tables/04_cross_approach/`. Figures are written to `figures/04_cross_approach/`.

