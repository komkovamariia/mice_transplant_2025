# Figures

This is the only figure root in the repository. Every Approach 1 plot displayed by the
restored notebook is saved as a 300-dpi PNG and an editable PDF.

Generated files follow this layout:

```text
figures/
├── 01_set_count/
│   ├── <stratum>/
│   └── cross_stratum/
├── 02_sequence_embedding/<stratum>/
├── 03_clone_alloreactivity/<stratum>/
└── 04_cross_approach/<stratum>/
```

Each Approach 1 stratum writes its figure inventory to
`results/01_set_count/<stratum>/figure_manifest.csv`. The runner verifies that every
listed PNG and PDF exists before reporting success.
