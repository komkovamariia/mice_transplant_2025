# Contributing

Use a branch and a focused pull request describing the biological or computational
problem, resulting behavior, validation and any limits of the evidence.

Before submitting changes, run:

```bash
python scripts/normalize_comments.py
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

Keep source notebooks free of outputs and execution counts. Preserve stable cell IDs
and existing output names where practical. Write comments, axes and documentation in
English. Add regression coverage when changing biological selection, statistical
calculations or data propagation.

Code comments should be concise, begin with a lower-case ordinary word when possible
and omit a full stop at the end of the final sentence. Scientific abbreviations and
proper names retain their normal capitalization. Move extended biological,
statistical and implementation rationale to Markdown under `docs/`. The complete
writing and naming convention is documented in [Repository writing and naming style](docs/repository_style.md).

Changes to exact aaV identity, sample exclusions, mouse pooling, contrasts, graphical
candidate thresholds or score definitions need a documented scientific rationale.
Do not infer missing metadata, silently relax thresholds or substitute simulated
results for study results. Synthetic tests must remain clearly identified as tests.

Report failures with the repository revision, relevant command, environment details
and the failing notebook cell or log excerpt. Keep input data and generated analysis
outputs outside source commits unless their publication is explicitly intended.
