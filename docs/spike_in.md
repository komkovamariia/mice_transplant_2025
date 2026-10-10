# Spike-in sensitivity protocol

The spike-in experiment measures how an added TRA family moves through the Approach 1 pipeline in the `all_combined` stratum. It is a computational sensitivity test on observed count matrices. It is not an assay-wide limit of detection and does not assign antigen specificity or biological alloreactivity.

## Baseline and target family

The runner first executes the unmodified combined baseline. TRA and TRB group matrices are saved with checksums and sample metadata.

A target V segment must be outside the baseline graphical candidate set defined by:

```text
frequency_in_full_g1 >= 0.006
retained_cdr3_share >= 0.66
```

Segments already supported by positive baseline edgeR or Fisher evidence in the final ranking are excluded.

Eligible family members must:

1. be observed exact `(CDR3 amino-acid sequence, V segment)` keys;
2. share the selected TRAV annotation;
3. be absent from g2, g4, g5 and g6;
4. have baseline UMI frequency at most `1e-5` in every g1 sample.

The default family size is 10. Selection prefers clonotypes absent from g1, so observed g3 clonotypes can supply g1-absent, control-absent keys. The default seed is 1031. The selected TRAV and aaV keys are fixed before any dose is run.

If no family meets the requested criteria, the experiment stops after the baseline. It does not create synthetic sequences or swap in a different enriched TRAV.

## Dose definition

For g1 sample `s`, let `L_s` be the original TRA UMI library size, `f` the requested family fraction and `K` the family size.

```text
B_s = floor(f * L_s + 0.5)
```

`B_s` is the total number of added family UMI for that sample. Counts are split as evenly as possible across the ordered family; the integer remainder is assigned deterministically.

The recorded dose includes both:

```text
B_s / L_s
B_s / (L_s + B_s)
```

The sequential CLI defaults to `0.0001,0.001,0.01`, which correspond to 0.01%, 0.1% and 1% of the original g1 TRA library. These are family-level doses, not per-clonotype doses.

Each dose starts from the same baseline. Original counts remain in place. At very small fractions, rounding can produce zero added UMI or fewer added UMI than selected clonotypes.

All g1 libraries receive the same requested fraction. Control groups and TRB matrices are unchanged.

## Structural and differential branches

The perturbed functional UMI matrices feed both branches of the same source notebook.

Structural analysis keeps the baseline library-eligibility set and rarefies each eligible library to 15,000 UMI before exact-set construction. A spike-in therefore cannot make a previously low-depth library enter the structural analysis.

Fisher and edgeR use the perturbed raw integer UMI counts. For the g1 versus g5+g6 differential analysis, samples are pooled within mouse after the perturbation. g2 remains outside that contrast.

Every dose reruns the same retention, Fisher, edgeR, ranking and figure code as the baseline.

## Recorded measurements

| Question | Output |
|---|---|
| Are the selected aaV recovered in the g1 remainder? | Presence, remainder membership and recovery fraction per aaV |
| Does the selected TRAV retain more unique aaV? | Remainder count and retention change |
| Does its g1 share increase? | Unique-aaV share and UMI share, with baseline differences |
| Does edgeR recover the added family? | Tested status, logFC, FDR and number with logFC > 0 and FDR < 0.05 |
| Does the TRAV move in the final ranking? | Rank, rank improvement and entry from an unranked baseline |

The structural article score is recomputed as:

```text
S(v) = f_full(v) * r(v)
```

The final evidence ranking is the existing Approach 1 ranking and is not tuned for the perturbation.

A selected clonotype can be absent from the edgeR test after expression filtering. That case is recorded separately from a tested clonotype with nonsignificant FDR.

## Sensitivity summary

`spike_in_summary.csv` contains the baseline and one row per dose.

`sensitivity_bounds.json` records the lowest tested passing fraction, the largest tested failure below the first pass and all observed transitions for:

- recovery of all selected aaV in the g1 remainder;
- positive FDR-significant edgeR recovery of any selected aaV;
- positive FDR-significant edgeR recovery of all selected aaV.

These values describe the tested grid. If no dose passes, no passing threshold is reported. If all doses pass, the lower transition lies below the tested range. A higher-dose failure after a lower-dose pass is recorded as a nonmonotone response.

Unique-aaV metrics can plateau once all selected aaV are present, even when their UMI counts continue to increase. edgeR and the final rank also depend on filtering, dispersion, replicate structure, competing features and TMM normalization.

The experiment currently uses one fixed family and uniform within-sample dosing. Other families or heterogeneous prevalence patterns need separate runs.

## Outputs

A run writes:

```text
results/01_spike_in/<run_id>/
figures/01_spike_in/<run_id>/
audit_runs/spike_in/<run_id>/
logs/spike_in/<run_id>/
```

The baseline directory contains the matrix snapshots and selection files. Each dose directory contains `spike_dose_by_sample.csv`, `spike_clonotype_recovery.csv` and `spike_metrics.json`. The run root contains `spike_in_summary.csv`, `sensitivity_bounds.json` and the final conclusion file.

For Slurm execution, see [aldan3_parallel.md](aldan3_parallel.md). The scheduler default uses 16 fractions from `1e-7` to `0.01`, corresponding to 0.00001% through 1%. Scheduling changes how cases are launched, not how they are analyzed.
