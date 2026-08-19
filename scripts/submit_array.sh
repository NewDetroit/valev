#!/usr/bin/env bash
# Compute a SLURM array bound from an accession list and submit.
#
#   submit_array.sh <sbatch-file> <accession-list> [max-concurrent] [chunk]
#
# Exists because of CORRECTION C-10: PLAN.md rewrites the #SBATCH --array line
# inside the script with sed, which silently no-ops once the line has been edited
# by hand — and its own instructions tell you to edit it by hand. Passing the
# range on the command line cannot fail that way, and `sbatch --array` overrides
# any in-file directive.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

SBATCH_FILE="${1:?usage: submit_array.sh <sbatch-file> <accession-list> [max-concurrent] [chunk]}"
ACC_LIST="${2:?usage: submit_array.sh <sbatch-file> <accession-list> [max-concurrent] [chunk]}"
MAXC="${3:-100}"
CHUNK="${4:-1}"

[ -f "$SBATCH_FILE" ] || { echo "no such sbatch file: $SBATCH_FILE" >&2; exit 2; }
[ -s "$ACC_LIST" ]    || { echo "empty or missing accession list: $ACC_LIST" >&2; exit 2; }

N=$(grep -cve '^[[:space:]]*$' "$ACC_LIST")
TASKS=$(( (N + CHUNK - 1) / CHUNK ))
[ "$TASKS" -ge 1 ] || { echo "computed 0 array tasks from $N accessions" >&2; exit 1; }

echo "accessions : $N"
echo "chunk      : $CHUNK per task"
echo "array      : 1-${TASKS}%${MAXC}"
echo "script     : $SBATCH_FILE"

if [ "${DRY_RUN:-0}" = "1" ]; then
  echo "DRY_RUN=1 — would run:"
  echo "  sbatch --array=1-${TASKS}%${MAXC} --export=ALL,ACC_LIST=$(readlink -f "$ACC_LIST"),CHUNK=${CHUNK} $SBATCH_FILE"
  exit 0
fi

command -v sbatch >/dev/null 2>&1 || { echo "sbatch not found — not on a SLURM cluster" >&2; exit 2; }
exec sbatch --array="1-${TASKS}%${MAXC}" \
     --export="ALL,ACC_LIST=$(readlink -f "$ACC_LIST"),CHUNK=${CHUNK}" \
     "$SBATCH_FILE"
