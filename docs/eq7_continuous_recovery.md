# Joint recovery: binary versus continuous evidence

Both ridge fits predict the same g1 minus mean control contrast. The
mouse-level index is the geometric mean of two quantities:

* Fraction of distinct aaV in the observed V segment recovered from
  the retained g1 reference after g5/g6 subtraction.
* Recovered UMI per million functional UMI. Each aaV is weighted by
  its prevalence in the retained g1 reference mice.

The two fits keep the same adaptive shortlist and log residual count.
The binary model uses the original aaV-level edgeR and Fisher flags.
The continuous model derives V-usage log2 fold change from mouse-level
edgeR and a log2 odds ratio from pooled-count Fisher with a 0.5 correction.
Each signed effect is multiplied by one minus its BH FDR. Clonotype
counts are summed by V before testing. This comparison changes both
the aggregation level and the encoding of the quantitative predictors.

The 75 held-out triplets per chain reuse 13 independent animals.
All libraries from each held-out mouse are excluded from features,
cutoffs and training outcomes. Ridge uses penalties 1 and 10.
The previously run model and its output directories remain available.

## Run on Aldan-3

```bash
conda activate mice-transplant-2025
cd ~/mice_transplant_env_test/mice_transplant_2025
git fetch origin eq7-continuous-recovery
git checkout --detach FETCH_HEAD
bench="$HOME/projects/mice_transplant_2025/mice_transplant_git/results/eq7_benchmark_1249"
python scripts/run_eq7_continuous.py --cache "$bench/_inputs" \
  --output results/eq7_continuous_1249
python scripts/plot_eq7_continuous.py --results results/eq7_continuous_1249 \
  --benchmark "$bench" --output figures/eq7_continuous_1249
```

The first execution fits V-usage edgeR and Fisher on raw counts and
stores their fold-level features in effect_cache. It loads the existing
seed 1249/1031 and depth 15000 inputs. Repeat full runs reuse the cached
V-usage effects. Use `--device cuda` with CUDA-enabled PyTorch for the
ridge solve. Count processing and tests run on CPU.

For external scoring add the following to the model command, with
article_outputs pointing to the extracted Overleaf analysis_outputs:

```bash
--article-inputs "$article_outputs" \
--skin results/eq7_external/PRJNA954716 \
--marrow results/eq7_external/PRJNA650133
```

The external directories must have the original TRUST4 reports. The
script validates report checksums against the archived QC. The matched
spleen shortlists and model parameters are fixed from internal data. V
counts are normalized by the sequencing totals of each mouse. External
counts are productive read fragments. PRJNA650133 has only the graft
arm, so it yields descriptive recovery without a contrast error.
