# Repository writing and naming style

This repository is maintained as research software accompanying a scientific study. User-facing text should read like material prepared for a manuscript or a reproducible methods supplement: concise, specific and biologically interpretable.

## Authors

Use the following author list and order wherever repository authorship is stated:

**Komkova M., Andreev V., Chernov P., Kofiadi I.**

Do not replace this list with generic contributor labels.

## Code comments

Keep executable files focused on implementation. Short comments should explain why a step is required, especially when the reason is biological, statistical or related to reproducibility.

For ordinary source-code comments:

- begin the first ordinary word with a lower-case letter when grammatically possible
- preserve scientific abbreviations and proper names such as UMI, TRA, TRB, edgeR, MiXCR and SLURM
- do not end the final sentence with a full stop
- avoid restating the code directly
- keep comments short enough to remain useful during code review
- move extended rationale, protocol descriptions and interpretation to a Markdown document under `docs/`

Technical directives such as shebangs, `#SBATCH`, `# noqa`, type-checker directives and formatter directives are exempt because their syntax is machine-readable.

The maintenance command below normalizes ordinary comments without altering executable statements:

```bash
python scripts/normalize_comments.py --write
```

Run it before repository validation when comments have been added or edited.

## Documentation

Long explanations belong in Markdown documents rather than inside scripts. Documentation should state the biological question first, then the computational method, parameters, outputs and interpretation limits.

Prefer titles such as:

- `Spike-in sensitivity analysis`
- `Primary exact-set and count-based analysis`
- `Sequence-space enrichment analysis`
- `Clone-level alloreactivity analysis`
- `Cross-method evidence comparison`

These titles are intended for readers of a scientific article. Internal executable filenames may remain stable when renaming them would break imports, scripts, notebooks, archived results or reproducibility links.

## Filenames and generated outputs

Use short descriptive filenames that reveal the scientific content without requiring knowledge of the implementation. Prefer names such as `spike_in_sensitivity.png`, `candidate_clones.csv` and `trav_ranking.csv` over temporary names, numbered scratch files or tool-generated labels.

Stable filenames that are already part of the reproducibility contract should not be renamed only for aesthetics. In those cases, improve the human-readable title, caption or README description while preserving the path.

## Validation

A wording-only cleanup must not change analytical behavior. Before merging style changes, run:

```bash
python scripts/normalize_comments.py
python scripts/validate_repository.py
python -m compileall -q src scripts
python -m pytest -q
```

Any change that affects biological selection, count propagation, statistical testing, filenames used by downstream code or output contracts requires separate scientific review and regression testing.
