# Splenic joint-recovery models

Each of the three ridge models is fitted for TRA and TRB in splenic CD4+, CD8+
and pooled CD4+ and CD8+ libraries. The pools are sums within the same mouse.
For a given chain, the same mouse set enters all three configurations. Each
mouse contributes only if both cell subsets have a library with at least 15,000
functional UMIs. Mice are held out in triplets with one mouse from each of g1,
g5 and g6. TRA uses 45 triplets, and TRB uses 75. Libraries, the structural
shortlist and the abundance evidence are processed within the corresponding
training set. Raw eligible counts enter edgeR and Fisher. Structural and
recovery analyses use 15,000 functional UMIs per library with seed 1031.

## Run

Use the environment in `environment.yml` and the existing analysis cache with
`mouse_metadata.csv` and `raw/count_snapshots`. For external evidence, pass the
two directories with TRUST4 reports in `trust4_reports`. Optional archived QC
tables under the manuscript's `analysis_outputs` directory verify report hashes.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
for chain in TRA TRB; do
  for subset in cd4 cd8 pooled; do
    python scripts/run_eq7_splenic.py \
      --cache "$input_cache" \
      --chain "$chain" --subset "$subset" \
      --skin "$skin_study" --marrow "$marrow_study" \
      --article-inputs "$archived_outputs" \
      --output "$run_dir/${chain}_${subset}" || exit
  done
done

python scripts/plot_eq7_splenic.py \
  --runs "$run_dir" --structural "$archived_trb_structure" --output "$figure_dir"
```

The six jobs can run independently. The pooled TRB job has the largest memory
footprint. The legacy all-tissue structural coordinates supplied to the plotting
command are included as `archived_trb_structure.csv` in the results archive.

The skin-graft comparison uses the splenic CD4+ TRB fit, a candidate list frozen
on the training cohort and a two-arm contrast. The training-cohort target
includes two transplant arms. In the public marrow cohort, recovery is measured
against pooled internal references for both chains. The cohort contains only
the recipient arm, so the contrast and its MAE are unavailable.

The outputs include the eligible libraries, per-split and per-V predictions,
rankings, error summaries at both penalties and frozen coefficients. The
manuscript reports the primary penalty of 1 in the figures and both 1 and 10 in
the error table. The original whole-tissue structural and selection benchmark
remains documented separately.

## Tissue panels and manuscript figures

Pass `--tissue thymus` or `--tissue combined` to run the same configuration on
thymic libraries or on spleen and thymus together. The combined panel requires
an eligible CD4+ and CD8+ library in each tissue for a mouse to enter that
panel. All eligible native and donor thymus libraries are included in the
thymus stratum as annotated in the input metadata. Counts from one mouse
are pooled before the abundance test. Mouse IDs are held out across all
libraries from that animal.

When a stratum has only two eligible intact mice, full-cohort abundance
predictors are fitted to recovery contrasts obtained from triplet omissions.
This yields calibration coefficients. There is no outer nested MAE for these
strata, because an outer fold would leave a single intact mouse and a further
intact-mouse omission would exhaust its reference. Their `run.json` records
`calibration_only` and the inputs to the final fit are written to
`calibration_only.csv`. Only splenic configurations are eligible for external
comparisons.

Recompute the spleen-only descriptive analyses and draw the article panels
using the new figure scripts. The same `input_cache` as above supplies the
raw MiXCR-derived count snapshots and metadata. Render with the Palatino
(P052) font available to Matplotlib.

```bash
for chain in TRA TRB; do
  python scripts/run_splenic_figures.py \
    --cache "$input_cache" --chain "$chain" --output "$figure_data/$chain"
done
python scripts/plot_splenic_manuscript.py \
  --data "$figure_data" --spleen-models "$run_dir" --output "$figure_dir"
```
