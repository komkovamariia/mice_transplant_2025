#!/usr/bin/env bash
#SBATCH --job-name=eq7_fisher_or
#SBATCH --partition=long
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00

set -euo pipefail

: "${EQ7_CACHE:?Set EQ7_CACHE to completed 1249 _inputs}"
: "${EQ7_PREVIOUS:?Set EQ7_PREVIOUS to previous complete aaV-fraction output}"
: "${EQ7_OUTPUT:?Set EQ7_OUTPUT to a NEW output folder}"

cd "${SLURM_SUBMIT_DIR}"
python_exe="${EQ7_PYTHON:-$(command -v python)}"
args=(--cache "${EQ7_CACHE}" --previous-results "${EQ7_PREVIOUS}"
      --output "${EQ7_OUTPUT}" --workers "${SLURM_CPUS_PER_TASK:-24}")
if [[ -n "${EQ7_ARTICLE_ZIP:-}" ]]; then
  [[ -f "${EQ7_ARTICLE_ZIP}" ]] || { echo "Missing article zip: ${EQ7_ARTICLE_ZIP}" >&2; exit 2; }
  args+=(--article-archive "${EQ7_ARTICLE_ZIP}")
fi

"${python_exe}" -m pytest -q tests/test_eq7_fisher_or.py
"${python_exe}" scripts/run_eq7_fisher_or_external.py "${args[@]}"
