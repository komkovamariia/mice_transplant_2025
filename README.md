# Mouse TCRα repertoire analysis after allogeneic transplantation

This repository provides a reproducible, sample-resolved analysis of mouse T-cell receptor alpha-chain repertoires in the BALB/c → C57BL/6 transplantation model. Three complementary methods are run independently and compared only after each has produced stratum-specific results.

## Analysis structure

| Order | Analysis | Primary question | Main outputs |
|---:|---|---|---|
| 1 | [`01_set_count`](approaches/01_set_count/) | Which exact aaV clonotypes and V segments differ between g1 and the allogeneic g5+g6 reference? | Set overlap, edgeR, Fisher cross-check |
| 2 | [`02_sequence_embedding`](approaches/02_sequence_embedding/) | Which local sequence-space regions and whole repertoires shift between the same groups? | Exact density enrichment, RFF-MMD, PERMANOVA, witness scores |
| 3 | [`03_clone_alloreactivity`](approaches/03_clone_alloreactivity/) | Which individual clonotypes receive repeated mouse-level support? | Prevalence, normalized abundance, convergence, candidate ranking |
| 4 | [`04_cross_approach`](approaches/04_cross_approach/) | Which V-segment signals recur across methods and biological strata, with concordant direction? | Rank overlap, signed-effect agreement, consensus ranking |

Approach 4 is an integration stage, not a fourth biological method.

## Biological strata

Every method produces an independent result for:

- `cd4_thymus`
- `cd4_spleen`
- `cd8_thymus`
- `cd8_spleen`
- `cd4_combined`
- `cd8_combined`

The combined strata pool thymus and spleen counts within each mouse. CD4 and CD8 repertoires are never pooled. This mouse-level aggregation prevents two tissues from the same animal from being treated as independent replicates.

## Repository layout

```text
.
├── approaches/
│   ├── 01_set_count/
│   ├── 02_sequence_embedding/
│   ├── 03_clone_alloreactivity/
│   └── 04_cross_approach/
├── src/                       # Shared analytical functions
├── scripts/
│   ├── run_analysis.py        # Cell-aware pipeline runner
│   └── validate_repository.py
├── tests/                     # Lightweight regression tests
├── figures/                   # The only figure root
│   ├── 01_set_count/
│   ├── 02_sequence_embedding/
│   ├── 03_clone_alloreactivity/
│   └── 04_cross_approach/
├── outputs/
│   ├── tables/
│   ├── notebooks/
│   ├── logs/
│   └── cache/
├── docs/
│   └── notebook_input_provenance.md
├── environment.yml
├── requirements.txt
└── RUN_GUIDE.md
```

Generated figures are saved as 300-dpi PNG and vector PDF files under `figures/<approach>/<stratum>/`. A plot that exists only inside an executed notebook is not considered a pipeline output.

## Notebook-specific inputs

The historical notebooks do not share one universal input file. The verified first pass, `venn_original.ipynb`, reads the existing project metadata and MiXCR clonoset index:

```text
/projects/mice_transplant_2025/metadata_mice_transplant.csv
/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv
```

The clonoset index contains the paths to the exported MiXCR clonotype tables under the project working directory. `repseq.intersections.count_table()` reads those tables and constructs the aaV matrices in memory. The first pass contains no `read_parquet()` call and does not require `clean_clonotypes_aaV.parquet`.

Parquet files occur in later reanalysis notebooks. In particular, the historical `mirpy_analysis.ipynb` uses a derived flattened table named `clean_clonotypes_aaV.parquet`, and later sections of Study 2 consume derived Parquet result layers. These are analysis artifacts rather than primary inputs to `venn_original.ipynb`.

The modular notebooks under `approaches/01_set_count`, `02_sequence_embedding`, and `03_clone_alloreactivity` currently call `src.strata.load_repertoire()` and therefore expect the derived flattened Parquet table. The repository does not yet provide a verified MiXCR-to-Parquet materialization command. The modular command must consequently be used only when this derived table has already been prepared; it must not replace the verified first-pass command below.

The complete notebook-by-notebook audit is recorded in [docs/notebook_input_provenance.md](docs/notebook_input_provenance.md).

## Environment

Create the complete Python/R environment:

```bash
conda env create -f environment.yml
conda activate mice-transplant-2025
```

Update an existing environment:

```bash
conda env update -n mice-transplant-2025 -f environment.yml --prune
conda activate mice-transplant-2025
```

Approach 1 uses edgeR through rpy2. The active pipeline does not use DESeq2.

## Execution

### Verified HPC first pass

For the audit checkout that contains `venn_original.ipynb`, run Approach 1 exactly as follows:

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

This command uses the metadata CSV, the clonoset-index CSV, and the MiXCR exports referenced by that index. It does not use Parquet.

The DESeq2-free historical notebook can be restored from repository history if it is absent from a clean checkout:

```bash
git show 6678766a0e8913ec2feacb3efb2daebd6a27eda4:venn_original.ipynb \
  > venn_original.ipynb
```

### Modular pipeline

Validate the repository first:

```bash
python scripts/validate_repository.py
python -m pytest -q
```

Run modular Approach 1 only when the derived `clean_clonotypes_aaV.parquet` table already exists:

```bash
python scripts/run_analysis.py --approach 1
```

Run Approach 2:

```bash
python scripts/run_analysis.py --approach 2
```

Run Approach 3:

```bash
python scripts/run_analysis.py --approach 3
```

Compare completed outputs from Approaches 1–3:

```bash
python scripts/run_analysis.py --approach 4
```

Run the complete pipeline in publication order:

```bash
python scripts/run_analysis.py --approach all
```

Run selected approaches:

```bash
python scripts/run_analysis.py --approach 1,3
```

Run selected strata:

```bash
python scripts/run_analysis.py \
  --approach 2 \
  --strata cd4_thymus,cd8_thymus
```

Pass the input directory explicitly:

```bash
python scripts/run_analysis.py \
  --approach all \
  --data-dir /absolute/path/to/derived_data
```

This option points the modular runner to a directory containing the derived `clean_clonotypes_aaV.parquet`. It is unrelated to the existing `/projects/mice_transplant_2025/test_run/clonosets_mice_transplant_2025_df.csv` input used by the verified first pass.

Resume a failed long run:

```bash
python scripts/run_analysis.py --approach 2 --resume
```

The runner restores saved notebook outputs, re-executes bootstrap cells, and skips completed analytical cells. Approach 2 also validates and reuses its sequence-space cache when the clonotype fingerprint is unchanged.

## Progress logs

The runner writes a pipeline log and one log per approach:

```text
outputs/logs/pipeline.log
outputs/logs/<approach>.log
```

Each code cell records its position, ID, first executable line, status, elapsed time, and remaining cell count. Long-running cells emit a heartbeat every 30 seconds.

Follow the full pipeline:

```bash
tail -f outputs/logs/pipeline.log
```

Follow Approach 2 at cell level:

```bash
tail -f outputs/logs/02_sequence_embedding.log
```

## Outputs

Each approach writes to:

```text
outputs/tables/<approach>/
```

Every stratum produces:

- a standardized V-segment ranking with `v_gene`, `rank`, `score`, and `stratum`;
- approach-specific evidence tables;
- `<stratum>_conclusion.json`;
- `<stratum>_conclusion.md`;
- publication-oriented figures under `figures/<approach>/<stratum>/`.

The cross-approach stage also writes `overall_consensus_v_genes.csv` and `cross_stratum_conclusion.md`.

## Statistical interpretation

Approach 1 models g1 versus g5+g6 with edgeR quasi-likelihood inference and uses Fisher's exact test as an independent count-based check.

Approach 2 uses one shared technical TCREmp/PCA basis. Biological tests remain stratum-specific. Exact k-d-tree density enrichment is complemented by RFF-MMD, a deterministic 9,999-permutation two-group PERMANOVA, and witness scoring.

Approach 3 weights mice equally. Its abundance component is calculated from within-mouse relative abundance, preventing larger libraries or larger experimental groups from dominating the integrated score. The score prioritizes candidates for follow-up and is not a calibrated probability of antigen specificity.

Approach 4 separates rank overlap from effect direction. All method-specific effects are represented on the same g1-versus-allogeneic axis, so opposing directions remain explicit.

Historical PERMANOVA values from the earlier notebook implementation must not be reused because the previous helper ignored permuted labels. Only values regenerated by `src/runtime.py` are valid.

Further operational details are provided in [RUN_GUIDE.md](RUN_GUIDE.md) and [docs/parallel_execution.md](docs/parallel_execution.md).
