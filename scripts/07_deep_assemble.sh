#!/usr/bin/env bash
# =============================================================================
# Task 9 — deep assembly worker.  Usage: 07_deep_assemble.sh <ACCESSION>
#
# The published Zheludev et al. pipeline: fasterq-dump -> fastp -> rnaSPAdes.
# Runs ONLY on accessions that survived triage. That reduction -- millions down
# to hundreds -- is the entire point of the hybrid architecture.
#
# Writes : $WORK/deep/<ACC>/transcripts.fasta
#          $WORK/deep/<ACC>/.assembled     completion sentinel (C-05)
#
# `set -e` is omitted deliberately: one bad accession must not kill the array task.
# =============================================================================
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

ACC="${1:?usage: $(basename "$0") <ACCESSION>}"
D="$WORK/deep/$ACC"

# C-05: guard on the sentinel, not on file non-emptiness.
[ -e "$D/.assembled" ] && { echo "$ACC skip (already assembled)"; exit 0; }
mkdir -p "$D" || exit 1
cd "$D" || exit 1

# rnaSPAdes treats a non-empty output directory as a restart candidate, so the
# assembly gets its own subdirectory rather than sharing one with the reads.
# PLAN.md runs `rnaspades.py -o .` from inside the directory holding its own
# input FASTQs, which invites a spurious restart on any re-run.
ASM="$D/spades"
MEM_GB="${SLURM_MEM_PER_NODE:+$(( SLURM_MEM_PER_NODE / 1024 ))}"
MEM_GB="${MEM_GB:-${RNASPADES_MEM_GB:-64}}"   # derived from the allocation so the two cannot drift

echo "=== $ACC (mem ${MEM_GB}G, threads $THREADS) ==="
if ! fasterq-dump --split-3 -e "$THREADS" -O . "$ACC"; then
  echo "$ACC FASTERQ_FAIL"; exit 1
fi

# Published parameters. Scientific constants -- changing them forfeits the
# replication claim.
FASTP_ARGS=(--average_qual=30 --n_base_limit=0 --cut_front --cut_tail -j fastp.json -h fastp.html)

if [ -s "${ACC}_1.fastq" ] && [ -s "${ACC}_2.fastq" ]; then
  fastp -i "${ACC}_1.fastq" -I "${ACC}_2.fastq" -o r1.clean.fastq -O r2.clean.fastq \
        "${FASTP_ARGS[@]}" || { echo "$ACC FASTP_FAIL"; exit 1; }
  SPADES_IN=(-1 r1.clean.fastq -2 r2.clean.fastq)
  # --split-3 emits an orphan file for unpaired mates; PLAN.md ignores it.
  if [ -s "${ACC}.fastq" ]; then
    fastp -i "${ACC}.fastq" -o rs.clean.fastq \
          --average_qual=30 --n_base_limit=0 --cut_front --cut_tail \
          -j fastp.se.json -h fastp.se.html && SPADES_IN+=(-s rs.clean.fastq)
  fi
elif [ -s "${ACC}.fastq" ]; then
  fastp -i "${ACC}.fastq" -o r.clean.fastq "${FASTP_ARGS[@]}" \
    || { echo "$ACC FASTP_FAIL"; exit 1; }
  SPADES_IN=(-s r.clean.fastq)
else
  echo "$ACC NO_READS"; exit 1
fi

rnaspades.py "${SPADES_IN[@]}" -o "$ASM" -t "$THREADS" -m "$MEM_GB" \
  || { echo "$ACC ASSEMBLY_FAIL"; exit 1; }

if [ -s "$ASM/transcripts.fasta" ]; then
  cp "$ASM/transcripts.fasta" "$D/transcripts.fasta"
  # Only now are the reads expendable. PLAN.md deletes them unconditionally, so a
  # failed assembly forces a full re-download on retry.
  rm -f "$D"/*.fastq "$ASM"/*.fastq
  : > "$D/.assembled"
  echo "$ACC OK $(grep -c '^>' "$D/transcripts.fasta") transcripts"
else
  echo "$ACC ASSEMBLY_FAIL (no transcripts.fasta; reads kept for retry)"
  exit 1
fi
