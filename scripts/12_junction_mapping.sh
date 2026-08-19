#!/usr/bin/env bash
# =============================================================================
# Task 12 Step 4 — physical evidence of circularity by junction read-mapping.
#
#   Usage: 12_junction_mapping.sh [--top N] [--candidate ID --accession ACC]
#
# THE STRONGEST EVIDENCE IN THE PROJECT. Terminal-repeat detection (scripts/09)
# is a heuristic about assembly output. This is a measurement on reads: if the
# molecule is genuinely circular, reads exist that cross the point where the
# assembler cut it open. Build a doubled reference S+S and those reads align
# across position L. On a linear molecule, none do.
#
# WHAT THE PLAN GETS WRONG HERE. PLAN.md ends this step with
#     samtools depth -a cand.bam | awk '{s+=$3; n++} END {print "mean depth:", s/n}'
# Mean depth over a doubled reference says nothing about circularity -- a purely
# linear molecule mapped to S+S also yields a perfectly good mean depth. The
# quantity that matters is the count of reads whose alignment SPANS position L.
# That is what this script reports.
#
# Reads  : $RESULTS/validated.fna, $RESULTS/replication_matrix.tsv
# Writes : $RESULTS/junction_support.tsv, $RESULTS/figures/junction/
# Exit: 0 at least one candidate has junction support · 2 missing input · 3 none do
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

TOP=3; ONE_CAND=""; ONE_ACC=""
while [ $# -gt 0 ]; do case "$1" in
  --top) TOP="${2:?}"; shift ;;
  --candidate) ONE_CAND="${2:?}"; shift ;;
  --accession) ONE_ACC="${2:?}"; shift ;;
  -h|--help) sed -n '2,26p' "$0"; exit 0 ;;
  *) echo "unknown argument: $1" >&2; exit 1 ;;
esac; shift; done

for t in bowtie2-build bowtie2 samtools seqkit fasterq-dump; do
  command -v "$t" >/dev/null 2>&1 || { echo "missing tool: $t — conda activate obelisk" >&2; exit 2; }
done
[ -s "$RESULTS/validated.fna" ] || { echo "missing $RESULTS/validated.fna — run scripts/11 then seqkit grep" >&2; exit 2; }
MATRIX="$RESULTS/replication_matrix.tsv"
[ -s "$MATRIX" ] || { echo "missing $MATRIX — run scripts/11_replication_matrix.py" >&2; exit 2; }

JD="$WORK/junction"; mkdir -p "$JD" "$RESULTS/figures/junction"
OUT="$RESULTS/junction_support.tsv"
printf 'candidate\taccession\tlength\tmapped_reads\tjunction_spanning_reads\tmean_depth\tverdict\n' > "$OUT"

# Resolve candidate -> accession from the matrix rather than by hand, which is how
# the plan leaves it (CAND=<candidate_id>; ACC=<an accession containing it>).
if [ -n "$ONE_CAND" ]; then
  PAIRS="$ONE_CAND	${ONE_ACC:-$(awk -F'\t' -v c="$ONE_CAND" 'NR>1 && $1==c {split($8,a,";"); print a[1]}' "$MATRIX")}"
else
  PAIRS=$(awk -F'\t' 'NR>1 {split($8,a,";"); print $1"\t"a[1]}' "$MATRIX" | head -"$TOP")
fi
[ -n "$PAIRS" ] || { echo "no candidate/accession pairs resolved from $MATRIX" >&2; exit 2; }

ANY=0
while IFS=$'\t' read -r CAND ACC; do
  [ -n "$CAND" ] && [ -n "$ACC" ] || continue
  echo "=== $CAND  (reads from $ACC) ==="
  cd "$JD"
  seqkit grep -p "$CAND" "$RESULTS/validated.fna" > "$CAND.fna"
  [ -s "$CAND.fna" ] || { echo "  $CAND not in validated.fna — skipping"; continue; }

  L=$(seqkit fx2tab -nl "$CAND.fna" | awk '{print $NF}' | head -1)
  python3 - "$CAND.fna" "$CAND.doubled.fna" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
name, seq = None, []
for line in open(src):
    if line.startswith(">"):
        name = line[1:].split()[0]
    else:
        seq.append(line.strip())
s = "".join(seq)
with open(dst, "w") as fh:
    fh.write(f">{name}_doubled\n{s*2}\n")
PY

  bowtie2-build -q "$CAND.doubled.fna" "${CAND}_idx"
  if [ ! -s "${ACC}_1.fastq" ] && [ ! -s "${ACC}.fastq" ]; then
    fasterq-dump --split-3 -e "$THREADS" -O . "$ACC" || { echo "  fasterq failed for $ACC"; continue; }
  fi
  if [ -s "${ACC}_1.fastq" ] && [ -s "${ACC}_2.fastq" ]; then
    bowtie2 -x "${CAND}_idx" -1 "${ACC}_1.fastq" -2 "${ACC}_2.fastq" -p "$THREADS" 2>"${CAND}.bt2.log"
  else
    bowtie2 -x "${CAND}_idx" -U "${ACC}.fastq" -p "$THREADS" 2>"${CAND}.bt2.log"
  fi | samtools sort -@ 2 -o "$CAND.bam" -
  samtools index "$CAND.bam"

  MAPPED=$(samtools view -c -F 4 "$CAND.bam")
  # The measurement the plan omits: alignments that start before L and end after it.
  SPAN=$(samtools view -F 4 "$CAND.bam" \
    | awk -v L="$L" '{
        pos=$4; cig=$6; len=0;
        while (match(cig, /^[0-9]+[MIDNSHP=X]/)) {
          n=substr(cig, RSTART, RLENGTH-1)+0; op=substr(cig, RSTART+RLENGTH-1, 1);
          if (op ~ /[MDN=X]/) len+=n;
          cig=substr(cig, RSTART+RLENGTH);
        }
        end=pos+len-1;
        if (pos <= L && end > L) c++
      } END {print c+0}')
  DEPTH=$(samtools depth -a "$CAND.bam" | awk '{s+=$3; n++} END {printf "%.2f", (n? s/n : 0)}')

  if [ "$SPAN" -ge 2 ]; then VERDICT="CIRCULAR_SUPPORTED"; ANY=1
  elif [ "$SPAN" -eq 1 ]; then VERDICT="WEAK_single_read"
  else VERDICT="NO_JUNCTION_SUPPORT"; fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$CAND" "$ACC" "$L" "$MAPPED" "$SPAN" "$DEPTH" "$VERDICT" >> "$OUT"
  echo "  length $L | mapped $MAPPED | junction-spanning $SPAN | mean depth $DEPTH | $VERDICT"

  samtools depth -a "$CAND.bam" > "$RESULTS/figures/junction/${CAND}.depth.tsv"
  cd - >/dev/null
done <<< "$PAIRS"

echo
column -t -s$'\t' "$OUT" 2>/dev/null || cat "$OUT"
echo
echo "wrote $OUT"
echo "Per-base depth for plotting: $RESULTS/figures/junction/<candidate>.depth.tsv"
echo "  Plot depth against position with a marker at L. Continuous coverage across L"
echo "  is the figure -- it is direct physical evidence, not a heuristic. Make it a poster panel."
[ "$ANY" -eq 1 ] || { echo; echo "No candidate has junction support (exit 3). Circularity is NOT demonstrated."; exit 3; }
