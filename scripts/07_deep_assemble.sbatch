#!/usr/bin/env bash
#SBATCH --job-name=deep-assemble
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%A_%a.out
# ADJUST for your cluster: --partition and --account.
# Submit via:  scripts/submit_array.sh scripts/07_deep_assemble.sbatch work/deep_accessions.txt 20 1
# The --array range is supplied there (C-10). rnaSPAdes -m is derived from
# SLURM_MEM_PER_NODE inside the worker, so --mem and -m cannot drift apart.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source config/config.sh

# C-11
if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null)" || CONDA_BASE=""
  [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ] \
    && . "$CONDA_BASE/etc/profile.d/conda.sh"
fi
conda activate obelisk 2>/dev/null || {
  echo "FATAL: cannot activate 'obelisk'. conda base: ${CONDA_BASE:-<not found>}" >&2; exit 2; }

: "${ACC_LIST:?ACC_LIST not set — submit via scripts/submit_array.sh}"
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACC_LIST")
[ -n "$ACC" ] && bash scripts/07_deep_assemble.sh "$ACC"
