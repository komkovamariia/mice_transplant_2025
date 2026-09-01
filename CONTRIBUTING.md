# Contributing

Use a branch and a focused pull request describing the biological or computational
problem, resulting behavior, validation and any limits of the evidence.

Before submitting changes, run:

```bash
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

Keep source notebooks free of outputs and execution counts. Preserve stable cell IDs
and existing output names where practical. Write comments, axes and documentation in
English. Add regression coverage when changing biological selection, statistical
calculations or data propagation.

Changes to exact aaV identity, sample exclusions, mouse pooling, contrasts, graphical
candidate thresholds or score definitions need a documented scientific rationale.
Do not infer missing metadata, silently relax thresholds or substitute simulated
results for study results. Synthetic tests must remain clearly identified as tests.

Report failures with the repository revision, relevant command, environment details
and the failing notebook cell or log excerpt. Keep input data and generated analysis
outputs outside source commits unless their publication is explicitly intended.
