#!/usr/bin/env bash
#SBATCH --job-name=logan-triage
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%A_%a.out
# ADJUST for your cluster: --partition and --account.
#
# The --array range is NOT set here. Submit via:
#     scripts/submit_array.sh scripts/05_logan_triage.sbatch work/accessions.txt 100
# which computes the bound from the accession file at submit time.
#
# CORRECTION C-10: PLAN.md Task 9 Step 5 patches the bound in with
#   sed -i "s/#SBATCH --array=1-100%20/.../"
# while its own Step 4 comment tells you to ADJUST that line by hand. Once you
# have, the literal pattern no longer matches, sed reports success, and you
# silently submit a truncated array.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source config/config.sh

# CORRECTION C-11: `source activate` / `conda activate` do not work in a
# non-interactive batch shell until conda's hook is sourced.
if [ -n "${CONDA_EXE:-}" ] || command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null)" || CONDA_BASE=""
  [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ] \
    && . "$CONDA_BASE/etc/profile.d/conda.sh"
fi
conda activate obelisk 2>/dev/null || {
  echo "FATAL: cannot activate the 'obelisk' conda environment." >&2
  echo "  conda base: ${CONDA_BASE:-<not found>}" >&2
  echo "  Create it with: conda env create -f env/environment.yml" >&2
  exit 2; }

: "${ACC_LIST:?ACC_LIST not set — submit via scripts/submit_array.sh}"
: "${CHUNK:=100}"

START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
echo "task $SLURM_ARRAY_TASK_ID: accessions $START-$END of $ACC_LIST"

sed -n "${START},${END}p" "$ACC_LIST" | while read -r A; do
  [ -n "$A" ] && bash scripts/05_logan_triage.sh "$A"
done
echo "task $SLURM_ARRAY_TASK_ID complete"
