#!/usr/bin/env bash
# =============================================================================
# Task 11 Step 4 — separate genuinely novel candidates from re-discoveries.
#
# Reads  : $RESULTS/candidates.fna, $REF_OBELISK_NT
# Writes : $RESULTS/candidates_novel.fna, $RESULTS/vs_known.tsv,
#          $RESULTS/rediscovered.txt
#
# Re-discovering known Obelisks is a GOOD sign: it independently validates the
# pipeline on data nobody curated for you. Report both numbers -- the ratio of
# re-discovered to novel is itself a credibility statistic, and a run that finds
# only novel things and no known ones should make you suspicious of your own
# pipeline rather than pleased.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

for t in makeblastdb blastn seqkit; do
  command -v "$t" >/dev/null 2>&1 || { echo "missing tool: $t" >&2; exit 2; }
done
[ -s "$RESULTS/candidates.fna" ] || { echo "missing $RESULTS/candidates.fna — run scripts/08_run_vnom.sh" >&2; exit 2; }
[ -s "$REF_OBELISK_NT" ]         || { echo "missing $REF_OBELISK_NT — run scripts/01_fetch_references.sh" >&2; exit 2; }

makeblastdb -in "$REF_OBELISK_NT" -dbtype nucl -out "$WORK/known_db" >/dev/null
blastn -query "$RESULTS/candidates.fna" -db "$WORK/known_db" \
       -outfmt "6 qseqid sseqid pident length evalue bitscore" \
       -evalue 1e-10 -num_threads "$THREADS" -out "$RESULTS/vs_known.tsv"

awk -F'\t' -v id="$KNOWN_IDENTITY_PCT" '$3 >= id {print $1}' "$RESULTS/vs_known.tsv" \
  | sort -u > "$RESULTS/rediscovered.txt"

TOTAL=$(grep -c '^>' "$RESULTS/candidates.fna")
KNOWN=$(wc -l < "$RESULTS/rediscovered.txt")

if [ "$KNOWN" -gt 0 ]; then
  seqkit grep -v -f "$RESULTS/rediscovered.txt" "$RESULTS/candidates.fna" > "$RESULTS/candidates_novel.fna"
else
  cp "$RESULTS/candidates.fna" "$RESULTS/candidates_novel.fna"
fi
NOVEL=$(grep -c '^>' "$RESULTS/candidates_novel.fna" || echo 0)

echo "candidates          : $TOTAL"
echo "re-discovered known : $KNOWN  (>= ${KNOWN_IDENTITY_PCT}% identity)"
echo "novel               : $NOVEL"
if [ "$KNOWN" -eq 0 ] && [ "$TOTAL" -gt 0 ]; then
  echo
  echo "NOTE: zero re-discoveries. Either this niche shares nothing with the published"
  echo "  set, or the reference in $REF_OBELISK_NT is not what you think it is."
  echo "  Check it before treating everything here as novel."
fi
seqkit stats "$RESULTS/candidates_novel.fna"
echo
echo "Next: scripts/11_replication_matrix.py   (GATE 3)"
