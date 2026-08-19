# Approach 1: exact sets and differential counts

This approach compares g1 with the allogeneic g5+g6 reference using exact aaV clonotypes.

For every biological stratum it reports exact g1/g5/g6 set overlap, edgeR quasi-likelihood differential abundance, Fisher's exact-test support, and a V-segment ranking. In the combined CD4 and CD8 strata, thymus and spleen counts are pooled within each mouse before model fitting.

## Verified historical HPC execution

The original first pass is `venn_original.ipynb`. It reads:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
```

The clonoset index points to the MiXCR export files. This execution path does not use Parquet.

```bash
mkdir -p audit_runs logs
jupyter nbconvert \
  --to notebook \
  --execute venn_original.ipynb \
  --ExecutePreprocessor.kernel_name=python3 \
  --ExecutePreprocessor.timeout=-1 \
  --output-dir audit_runs \
  --output 01_venn_original.executed.ipynb \
  2>&1 | tee logs/01_venn_original.log
```

## Modular execution

The concise notebook in this directory uses the later sample-resolved interface and requires a previously materialized `clean_clonotypes_aaV.parquet`:

```bash
python scripts/run_analysis.py --approach 1
```

The repository does not yet contain the MiXCR-to-Parquet materialization step, so this command must not be substituted for the verified historical first pass when the derived file is absent.

Tables and result-specific conclusions are written to `outputs/tables/01_set_count/`. Figures are written to `figures/01_set_count/<stratum>/`.
