#!/usr/bin/env bash
# =============================================================================
# Task 10 — VNom circular-candidate calling, then pooling.
#
# Reads  : $WORK/deep/<ACC>/transcripts.fasta
# Writes : $WORK/vnom/<ACC>/4_final_clusters/    VNom output (may be several FASTAs)
#          $RESULTS/candidates.fna               pooled, length-filtered, deduplicated
#          $RESULTS/candidate_sources.tsv        candidate_id -> every source accession
#
# ---------------------------------------------------------------------------
# CORRECTION C-15 — PLAN.md's VNom invocation cannot work. It is:
#
#     python "$VNOM/VNom.py" -i "$T" -max 2000 ... -o "$WORK/vnom/$ACC"
#
# where $T is a full path ending in transcripts.fasta. Checked against VNom's
# README (fetched 2026-08-19), whose own worked example is:
#
#     python ../VNom.py -i peach_subset -max 2000 -CF_k 10 -CF_simple 0 \
#            -CF_tandem 1 -USG_vs_all 1 > peach_subset_VNom.log
#
# Four incompatibilities:
#   1. `-i` takes the file's BASENAME WITHOUT the .fasta extension, not a path.
#      The README: "you must specify this single underscore name without the
#      file ending for VNom".
#   2. There is NO -o flag. Outputs are written to `4_final_clusters` relative to
#      the working directory. The README: "outputs are stored to 4_final_clusters
#      (so if this dir wasn't written, VNom failed to nominate viroid-like
#      contigs)". Isolation therefore requires running each accession in its own
#      directory, which is what this script does.
#   3. The filename must contain EXACTLY ONE underscore and end in .fasta —
#      "X_Y.fasta is good, but XY.fasta is bad". `transcripts.fasta` has none.
#   4. seqIDs must keep the default rnaSPAdes layout; the README warns that
#      "adding more underscores will cause VNom to crash". Its example replaces
#      the literal NODE with the accession, which is what we do.
#
# Also from the README, and absent from PLAN.md entirely: VNom requires input
# derived from STRANDED RNA-seq, and depends on circuclust, USEARCH and MARS
# being installed under VNom/dependencies/. An unstranded library will produce
# nothing at the dual-polarity filter, which the author describes as where VNom
# most often stops.
#
# PROVENANCE. `seqkit rmdup -s` correctly collapses one element seen in twenty
# libraries into one record, but discards which twenty — and the replication
# matrix at GATE 3 is built entirely on that mapping. It is captured before the
# dedup and written to candidate_sources.tsv. PLAN.md pools with a bare
# `cat | rmdup` and loses it.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

VNOM_PY="${VNOM_PY:-$VNOM/VNom.py}"
# Published Zheludev et al. parameters. Scientific constants — do not change.
VNOM_ARGS=(-max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1)

command -v seqkit >/dev/null 2>&1 || { echo "missing seqkit — conda activate obelisk" >&2; exit 2; }
[ -f "$VNOM_PY" ] || { echo "missing $VNOM_PY — clone VNom (see docs/RUNBOOK.md)" >&2; exit 2; }
[ -d "$WORK/deep" ] || { echo "no $WORK/deep — run scripts/07_deep_assemble.sh first" >&2; exit 2; }

mkdir -p "$WORK/vnom"
n_in=0; n_ok=0; n_empty=0; n_fail=0
for T in "$WORK"/deep/*/transcripts.fasta; do
  [ -e "$T" ] || continue
  ACC=$(basename "$(dirname "$T")")
  n_in=$((n_in + 1))
  D="$WORK/vnom/$ACC"
  [ -e "$D/.done" ] && continue
  echo "=== VNom $ACC ==="
  mkdir -p "$D"

  # Exactly one underscore, .fasta extension (C-15 point 3), and NODE replaced
  # with the accession so cluster members stay traceable (point 4).
  STEM="${ACC}_contigs"
  seqkit grep -v -s -p 'N' "$T" 2>/dev/null | sed "s/NODE/${ACC}/g" > "$D/${STEM}.fasta"
  if [ ! -s "$D/${STEM}.fasta" ]; then
    echo "  no N-free contigs for $ACC"; n_empty=$((n_empty + 1)); : > "$D/.done"; continue
  fi

  # No -o flag exists, so isolate by running inside the accession's directory
  # (C-15 point 2) and pass the bare stem (point 1).
  if ( cd "$D" && python "$VNOM_PY" -i "$STEM" "${VNOM_ARGS[@]}" > "${STEM}_VNom.log" 2>&1 ); then
    if [ -d "$D/4_final_clusters" ]; then
      : > "$D/.done"; n_ok=$((n_ok + 1))
      echo "  ok — $(find "$D/4_final_clusters" -name '*.fasta' | wc -l) cluster FASTA(s)"
    else
      # Per the README this is the normal "nothing nominated" outcome, most often
      # at the dual-polarity filter. Not an error.
      : > "$D/.done"; n_empty=$((n_empty + 1))
      echo "  no 4_final_clusters — nothing nominated (see $D/${STEM}_VNom.log)"
    fi
  else
    n_fail=$((n_fail + 1))
    echo "  VNom FAILED on $ACC — see $D/${STEM}_VNom.log"
  fi
done
echo
echo "VNom: $n_in assemblies | $n_ok nominated | $n_empty nothing nominated | $n_fail failed"
[ "$n_ok" -gt 0 ] || { echo "No accession produced candidates. Nothing to pool."; exit 3; }

echo
echo "=== pooling ==="
POOL="$WORK/pooled_raw.fna"; : > "$POOL"
shopt -s nullglob
n_files=0
for d in "$WORK"/vnom/*/; do
  ACC=$(basename "${d%/}")
  # Every FASTA under 4_final_clusters, not just the first — same class of defect
  # as C-08 in the replication matrix.
  for fa in "$d"4_final_clusters/*.fasta "$d"4_final_clusters/*.fa "$d"4_final_clusters/*.fna; do
    [ -e "$fa" ] || continue
    n_files=$((n_files + 1))
    seqkit replace -p '^' -r "${ACC}__" "$fa" >> "$POOL"
  done
done
shopt -u nullglob
echo "  read $n_files VNom FASTA(s)"
[ -s "$POOL" ] || { echo "pooled file is empty" >&2; exit 3; }

seqkit seq -m "$OBELISK_MIN_LEN" -M "$OBELISK_MAX_LEN" "$POOL" \
  | seqkit rmdup -s -D "$WORK/rmdup_duplicates.txt" \
  > "$RESULTS/candidates.fna"

python3 - <<'PY'
import os, re
from collections import defaultdict
WORK, RESULTS = os.environ["WORK"], os.environ["RESULTS"]

def ids(path):
    for line in open(path):
        if line.startswith(">"):
            yield line[1:].split()[0]

sources = defaultdict(set)
for i in ids(f"{RESULTS}/candidates.fna"):
    sources[i].add(i.split("__", 1)[0])

dup = f"{WORK}/rmdup_duplicates.txt"
if os.path.exists(dup):
    for line in open(dup):
        parts = [p.strip() for p in re.split(r"[,\t]", line.strip()) if p.strip()]
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
echo "Next: python3 scripts/09_circularity_check.py"
