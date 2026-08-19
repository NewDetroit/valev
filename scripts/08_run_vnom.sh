#!/usr/bin/env bash
# =============================================================================
# Task 10 — VNom circular-candidate calling, then pooling.
#
# Reads  : $WORK/deep/<ACC>/transcripts.fasta
# Writes : $WORK/vnom/<ACC>/                VNom output (MAY be several FASTAs)
#          $RESULTS/candidates.fna          pooled, length-filtered, deduplicated
#          $RESULTS/candidate_sources.tsv   candidate_id -> every source accession
#
# PROVENANCE. `seqkit rmdup -s` collapses identical sequences, which is what we
# want -- one element recovered from twenty libraries should count once -- but it
# discards which accessions those twenty were. Task 12's replication matrix is
# built entirely on that mapping, so it is captured BEFORE the dedup: every record
# is renamed <ACC>__<original-id> on the way in, and candidate_sources.tsv records
# the full set. PLAN.md pools with a bare `cat | rmdup` and loses it.
#
# VNom's CLI entry point could not be verified when this was written (its README
# was not reachable from the build environment; api.github.com is blocked there).
# The PARAMETERS below are from the published methods and are reliable. The
# INVOCATION is isolated in VNOM_CMD so there is exactly one line to correct.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

# The one line to fix if VNom's README says otherwise. Override without editing:
#   VNOM_CMD="python $VNOM/vnom/cli.py" bash scripts/08_run_vnom.sh
VNOM_CMD="${VNOM_CMD:-python $VNOM/VNom.py}"

# Published Zheludev et al. parameters. Do not change.
VNOM_ARGS=(-max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1)

command -v seqkit >/dev/null 2>&1 || { echo "missing seqkit — conda activate obelisk" >&2; exit 2; }
[ -d "$WORK/deep" ] || { echo "no $WORK/deep — run scripts/07_deep_assemble.sh first" >&2; exit 2; }

mkdir -p "$WORK/vnom"
n_in=0 n_ok=0 n_fail=0
for T in "$WORK"/deep/*/transcripts.fasta; do
  [ -e "$T" ] || continue
  ACC=$(basename "$(dirname "$T")")
  n_in=$((n_in+1))
  if [ -e "$WORK/vnom/$ACC/.done" ]; then continue; fi
  echo "=== VNom $ACC ==="
  mkdir -p "$WORK/vnom/$ACC"
  if $VNOM_CMD -i "$T" "${VNOM_ARGS[@]}" -o "$WORK/vnom/$ACC" 2>>"$LOGS/vnom.err"; then
    : > "$WORK/vnom/$ACC/.done"; n_ok=$((n_ok+1))
  else
    echo "  VNom failed on $ACC (see $LOGS/vnom.err)"; n_fail=$((n_fail+1))
  fi
done
echo "VNom: $n_in assemblies, $n_ok ok, $n_fail failed"

echo
echo "=== pooling candidates ==="
POOL="$WORK/pooled_raw.fna"
: > "$POOL"
shopt -s nullglob
n_files=0
for d in "$WORK"/vnom/*/; do
  ACC=$(basename "$d")
  # Every FASTA, not just the first — same class of defect as C-08.
  for fa in "$d"*.fasta "$d"*.fa "$d"*.fna; do
    [ -e "$fa" ] || continue
    n_files=$((n_files+1))
    seqkit replace -p '^' -r "${ACC}__" "$fa" >> "$POOL"
  done
done
shopt -u nullglob
echo "  read $n_files VNom FASTA(s)"

seqkit seq -m "$OBELISK_MIN_LEN" -M "$OBELISK_MAX_LEN" "$POOL" \
  | seqkit rmdup -s -D "$WORK/rmdup_duplicates.txt" \
  > "$RESULTS/candidates.fna"

# candidate_sources.tsv: the mapping rmdup would otherwise destroy.
python3 - <<'PY'
import os, re
from collections import defaultdict
WORK, RESULTS = os.environ["WORK"], os.environ["RESULTS"]

def ids(path):
    for line in open(path):
        if line.startswith(">"):
            yield line[1:].split()[0]

kept = list(ids(f"{RESULTS}/candidates.fna"))
sources = defaultdict(set)
for i in kept:
    sources[i].add(i.split("__", 1)[0])

# seqkit rmdup -D lists the ids collapsed into each representative.
dup = f"{WORK}/rmdup_duplicates.txt"
if os.path.exists(dup):
    for line in open(dup):
        parts = [p.strip() for p in re.split(r"[,\t]", line.strip()) if p.strip()]
        if len(parts) < 2:
            continue
        rep = next((p for p in parts if p in sources), None)
        if rep is None:
            continue
        for p in parts:
            sources[rep].add(p.split("__", 1)[0])

with open(f"{RESULTS}/candidate_sources.tsv", "w") as fh:
    fh.write("candidate\tn_accessions\taccessions\n")
    for c in sorted(sources):
        a = sorted(sources[c])
        fh.write(f"{c}\t{len(a)}\t{';'.join(a)}\n")
print(f"  wrote {RESULTS}/candidate_sources.tsv ({len(sources)} candidates)")
PY

seqkit stats "$RESULTS/candidates.fna"
echo
echo "Next: scripts/09_circularity_check.py"
