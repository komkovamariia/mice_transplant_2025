# Spike-in sensitivity protocol

The experiment asks whether the primary analysis recovers an introduced TRA signal
in the combined CD4 + CD8, thymus + spleen stratum. It is a computational sensitivity
experiment on observed count matrices. It does not establish biological alloreactivity
or an assay-wide limit of detection.

## Baseline and target selection

The runner first executes the complete `venn_original.ipynb` on the unmodified combined
stratum. It freezes the six group matrices for each chain, checksums their files and
records sample metadata identity. It selects the family only after baseline retention,
Fisher, edgeR and TRAV rankings are available.

A target V segment must fail the baseline graphical candidate rule of both
`frequency_in_full_g1 >= 0.006` and `retained_cdr3_share >= 0.66`. Segments with positive
baseline edgeR or Fisher enrichment in the existing final ranking are also excluded.
This definition distinguishes a candidate from a descriptive entry in the ranking.

Eligible family members:

1. Are exact `(CDR3 amino-acid sequence, V segment)` keys already observed in the input matrices.
2. Share the selected exact TRAV annotation. This does not imply common antigen specificity.
3. Are absent from the union of the subtraction groups g2, g4, g5 and g6.
4. Have baseline UMI frequency at most `1e-5` in every g1 sample, including absence.

The default selects ten family members, preferentially those absent from g1. Observed
g3 clonotypes can supply g1-absent, control-absent keys. The seed is 1031. Both the
chosen TRAV and exact family are fixed before any perturbed analysis and reused for
all doses. No targets are selected by their eventual recovery performance.

If no valid family is available, the experiment stops after preserving the baseline
and explains which explicit options can be revised. It never fabricates new sequences
or silently substitutes a currently enriched TRAV.

## Dose definition and allocation

Let `L_s` be the original TRA UMI library size of g1 sample `s`, `f` the requested
family fraction, and `K` the number of selected aaV clonotypes. The added family budget is

```text
B_s = floor(f * L_s + 0.5)
```

`B_s` is split as evenly as integer counts permit across the fixed ordered family.
Any remainder is assigned deterministically. A sample's final library is `L_s + B_s`.
The recorded dose includes both `B_s / L_s` and `B_s / (L_s + B_s)`.

The three default fractions are 0.0001, 0.001 and 0.01: 0.01%, 0.1% and 1%. These are
**total family doses**, not doses per clonotype. Every dose starts from the same
baseline. Original UMI counts remain in place; added counts increase library size.
No additional UMI filter, resampling or library-depth matching is applied.

At low depth, a requested dose may round to zero or supply fewer UMI than family
members. Such aaV remain absent until they actually receive a positive count. The
experiment does not force a minimum of one UMI per selected clonotype.

All g1 sample libraries receive the same requested fraction, before their counts are
pooled within mouse. The control matrices and all TRB matrices are unchanged. Standard
TMM normalization, expression filtering and edgeR inference are rerun for each dose.

## Reuse of the primary analysis

A baseline plus each dose executes the same full source notebook, including its
TRA/TRB figures, retention tables, Fisher analysis, edgeR and final TRAV aggregation.
Exact-set analysis and differential inference consume the same perturbed matrices.
The latter uses the union of g1, g2, g5 and g6 matrices and then pools samples within
mouse. g2 is excluded from the g1 versus g5+g6 edgeR contrast.

Snapshots are loaded for both chains in perturbed runs. Matrix checksums and metadata
identity are verified before reuse. `selection.json`, `experiment.json`, source hashes,
Python package versions and `R_sessionInfo.txt` make the experiment auditable.

## Hypotheses and measurements

| Question | Recorded measurements |
|---|---|
| Are selected aaV retained in the g1 remainder? | Presence, remainder membership and recovery fraction for every selected key |
| Does the TRAV remainder increase? | Unique aaV before/after subtraction, absolute remainder change and retention change |
| Does its share in g1 increase? | Unique-aaV share and UMI share, each with its own baseline difference |
| Does edgeR recover the signal? | Tested/filtered status, logFC, FDR and count with logFC > 0 and FDR < 0.05 |
| Does the TRAV rise in the final ranking? | Rank, rank improvement, entry from an unranked baseline and unchanged article score |

The article score is recomputed as `S(v) = f_full(v) * r(v)`. The final evidence
ranking retains the existing notebook aggregation; neither score is tuned for the
spike-in. A baseline-unranked TRAV has a missing rank, not an invented numeric rank.
Entry into the ranking is reported separately from a numeric rank change.

A failed hypothesis is a result, so it does not stop the experiment. Technical failure
such as a missing figure or unfitted edgeR model does stop it. A selected aaV that
fails expression filtering is distinguished from one tested with nonsignificant FDR.

## Sensitivity bounds and interpretation

`spike_in_summary.csv` contains one baseline row and one row per dose.
`sensitivity_bounds.json` reports the lowest passing tested dose, the largest failing
dose below it, and every observed pass/fail transition for:

- recovery of all selected aaV in the g1 remainder;
- positive, FDR-significant edgeR recovery of any selected aaV;
- positive, FDR-significant edgeR recovery of all selected aaV.

If every tested dose passes, the lower transition is below the tested range. If none
passes, a passing threshold is not established. The highest dose tested is not an
upper sensitivity limit. An upper failure is recorded only when a higher dose loses
recovery after a lower passing dose. A nonmonotone response is retained in the output.

Increasing UMI for an aaV already present in g1 does not increase its unique-aaV count.
Exact-set retention and the article score may therefore plateau even while UMI share
and edgeR evidence increase. A high baseline retention can also leave no room for an
increase. FDR and final rank depend on dispersion, replicate counts, competing genes
and TMM normalization; the hypotheses are not guaranteed by construction.

This experiment uses one family and uniform within-sample dosing. Results describe
that configuration. Other families, heterogeneous mouse prevalence and repeated
simulation designs require explicit additional experiments.

## Outputs and reruns

Independent doses can run as a bounded Slurm array after baseline completion; see
[Aldan-3 parallel execution](aldan3_parallel.md). The scheduler default expands the
grid to sixteen fractions from `1e-7` to `0.01` (0.00001% to 1%). The sequential CLI
retains its three-dose default. Scheduling does not change any selection or test rule.

Use the commands and file map in [RUN_GUIDE.md](../RUN_GUIDE.md). Each run uses a fresh
identifier under `results/01_spike_in/`, `figures/01_spike_in/`, `audit_runs/spike_in/`
and `logs/spike_in/`. Standard article results under `01_set_count` are preserved.
The source matrices remain in the baseline directory and are never overwritten by a dose.
