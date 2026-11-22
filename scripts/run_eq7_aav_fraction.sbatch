#!/usr/bin/env bash
#SBATCH --job-name=eq7_aav_fraction
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00

set -euo pipefail

: "${EQ7_CACHE:?Set EQ7_CACHE to the existing eq7_benchmark_1249/_inputs directory}"
: "${EQ7_OUTPUT:?Set EQ7_OUTPUT to a NEW result directory}"
cd "${SLURM_SUBMIT_DIR}"
python_exe="${EQ7_PYTHON:-$(command -v python)}"

args=(--cache "${EQ7_CACHE}" --output "${EQ7_OUTPUT}" --workers "${SLURM_CPUS_PER_TASK:-4}")
if [[ -n "${EQ7_REF:-}" && -f "${EQ7_REF}/internal_predictions.csv" ]]; then
    args+=(--reference-results "${EQ7_REF}")
fi

"${python_exe}" -m pytest -q tests/test_eq7_aav_fraction.py
"${python_exe}" scripts/run_eq7_aav_fraction.py "${args[@]}"
