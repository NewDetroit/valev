#!/usr/bin/env bash
# =============================================================================
# Task 4 — build the search databases and calibrate the false-positive control.
#
# Reads  : $REF_OBLIN_FAA   (Oblin-1 proteins, from scripts/01_fetch_references.sh)
#          $REF_OBLIN_HMM   (published Obelisk HMM, ditto)
# Writes : $REF_OBLIN_DMND.dmnd   DIAMOND db for translated triage
#          $REF_DECOY_FAA         residue-shuffled negative control
#          $REF_DECOY_DMND.dmnd   DIAMOND db of the decoys
#          $REF/calibration.tsv   the decoy check's numbers, for the write-up
#
# ---------------------------------------------------------------------------
# CORRECTION C-03 — the plan's decoy control is inverted. Read this before
# "simplifying" the Python block back into a seqkit one-liner.
#
# PLAN.md Task 4 Step 1 builds the decoy as:
#
#     seqkit shuffle -s 42 oblin1_centroids.faa \
#       | seqkit mutate --any-point 100 -s 42 > decoy_shuffled.faa || python - <<EOF
#
# Two defects that compound into the opposite of a control:
#
#   1. `seqkit shuffle` shuffles the ORDER OF RECORDS IN A FILE. It does not
#      permute residues within a sequence. Its output is the identical protein
#      set in a different order.
#   2. `seqkit mutate` has no `--any-point` flag.
#
# So every "decoy" is a verbatim real Oblin. Step 4 then aligns the decoys
# against the real database and instructs:
#
#     "Expected: zero or very few lines. ... If the decoy hits heavily, your
#      E-value threshold is too permissive. Tighten DIAMOND_EVALUE in config.sh
#      until the decoy is clean ... This calibration IS your false-positive
#      control -- do not skip it."
#
# Every decoy hits itself at E~0, which is a maximal and unfixable hit rate.
# Following the instruction drives DIAMOND_EVALUE toward zero until genuine
# Oblins stop matching, destroying the sensitivity of the entire triage sweep.
# The step that exists to control false positives instead guarantees false
# negatives -- and it does so silently, because the sweep still runs and still
# reports a number.
#
# The plan's own Python fallback shuffles residues correctly, but it sits behind
# `||`, so it executes only if the seqkit pipeline FAILS. `seqkit shuffle`
# succeeds. The correct code never runs.
#
# Here the residue shuffle is unconditional, and a dirty decoy is a hard failure
# (exit 3) rather than a prompt to loosen the science.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

FROM_SCRATCH=0
for a in "$@"; do
  case "$a" in
    --from-scratch) FROM_SCRATCH=1 ;;
    -h|--help) sed -n '2,50p' "$0"; exit 0 ;;
    *) echo "unknown argument: $a" >&2; exit 1 ;;
  esac
done

need() { command -v "$1" >/dev/null 2>&1 || { echo "missing tool: $1 — conda activate obelisk" >&2; exit 2; }; }
need diamond
need hmmsearch
need python3

[ -s "$REF_OBLIN_FAA" ] || { echo "missing $REF_OBLIN_FAA — run scripts/01_fetch_references.sh" >&2; exit 2; }

echo "=== 1/4  DIAMOND database from Oblin-1 proteins ==="
diamond makedb --in "$REF_OBLIN_FAA" --db "$REF_OBLIN_DMND" --quiet
echo "  $REF_OBLIN_DMND.dmnd"

echo
echo "=== 2/4  profile HMM ==="
if [ "$FROM_SCRATCH" -eq 1 ]; then
  # Fallback path only. Aligning ~1,700 divergent Oblin-1 sequences with
  # `mafft --auto` and building one model over the result gives a markedly worse
  # profile than the published one, and its scores are not comparable to anything
  # in the literature. Kept because C-02's primary path depends on a public S3
  # object that could move.
  need mafft
  need hmmbuild
  echo "  --from-scratch: aligning and building (this is the WORSE model — see C-02)"
  mafft --auto --thread "$THREADS" "$REF_OBLIN_FAA" > "$REF/oblin1_aln.faa"
  hmmbuild --amino "$REF_OBLIN_HMM" "$REF/oblin1_aln.faa" >/dev/null
else
  [ -s "$REF_OBLIN_HMM" ] || {
    echo "missing $REF_OBLIN_HMM — run scripts/01_fetch_references.sh, or pass --from-scratch" >&2
    exit 2; }
  echo "  using the published Obelisk HMM (C-02) — keeps scores comparable to the literature"
fi

# CORRECTION C-09: hmmpress builds binary indices for hmmscan. This pipeline only
# ever calls hmmsearch, which reads the plain .hmm. Pressed anyway so hmmscan stays
# available, with -f because a bare hmmpress exits non-zero when the .h3* files
# already exist, which breaks this script's `set -e` on any re-run.
hmmpress -f "$REF_OBLIN_HMM" >/dev/null 2>&1 || true
grep -q '^HMMER3' "$REF_OBLIN_HMM" || { echo "  $REF_OBLIN_HMM is not a HMMER3 profile" >&2; exit 1; }
echo "  $(grep -m1 '^LENG' "$REF_OBLIN_HMM" || true)   $(grep -c '^HMMER3' "$REF_OBLIN_HMM") model(s)"

echo
echo "=== 3/4  residue-shuffled decoys (C-03) ==="
python3 - <<'PY'
import os, random, sys

src = os.environ["REF_OBLIN_FAA"]
dst = os.environ["REF_DECOY_FAA"]

def read_fasta(path):
    rec_id, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if rec_id is not None:
                    yield rec_id, "".join(buf)
                rec_id, buf = line[1:].split()[0], []
            elif line:
                buf.append(line.strip())
    if rec_id is not None:
        yield rec_id, "".join(buf)

# Fisher-Yates over the residues of each record, independently, with a fixed seed.
# Per-record shuffling preserves each sequence's exact length and amino-acid
# composition while destroying motif order -- which is precisely the null the
# calibration needs. Shuffling the concatenation instead would preserve only the
# global composition and would let residues migrate between records.
rng = random.Random(42)
n = 0
with open(dst, "w") as out:
    for rec_id, seq in read_fasta(src):
        if not seq:
            continue
        res = list(seq)
        for i in range(len(res) - 1, 0, -1):
            j = rng.randrange(i + 1)
            res[i], res[j] = res[j], res[i]
        shuf = "".join(res)
        assert sorted(shuf) == sorted(seq), f"composition changed for {rec_id}"
        out.write(f">decoy_{rec_id}\n{shuf}\n")
        n += 1

if n == 0:
    sys.exit(f"[02_build_search_dbs] no sequences read from {src}")
print(f"  wrote {n} decoys to {dst} (length and composition preserved, order destroyed)")
PY

diamond makedb --in "$REF_DECOY_FAA" --db "$REF_DECOY_DMND" --quiet
echo "  $REF_DECOY_DMND.dmnd"

echo
echo "=== 4/4  calibration: the decoys must NOT match the real Oblins ==="
CAL="$WORK/decoy_selfcheck.tsv"
diamond blastp \
  --query "$REF_DECOY_FAA" --db "$REF_OBLIN_DMND" \
  --evalue "$DIAMOND_EVALUE" $DIAMOND_SENS --threads "$THREADS" \
  --outfmt 6 qseqid sseqid pident length evalue bitscore \
  --out "$CAL" --quiet

N_DECOY=$(grep -c '^>' "$REF_DECOY_FAA")
N_HITROWS=$(wc -l < "$CAL")
N_HITSEQS=$(cut -f1 "$CAL" 2>/dev/null | sort -u | wc -l)
RATE=$(python3 -c "print(f'{($N_HITSEQS/$N_DECOY*100) if $N_DECOY else 0:.3f}')")

{
  printf 'metric\tvalue\n'
  printf 'diamond_evalue\t%s\n' "$DIAMOND_EVALUE"
  printf 'diamond_sensitivity\t%s\n' "$DIAMOND_SENS"
  printf 'decoys\t%s\n' "$N_DECOY"
  printf 'decoy_hit_rows\t%s\n' "$N_HITROWS"
  printf 'decoy_hit_sequences\t%s\n' "$N_HITSEQS"
  printf 'decoy_hit_rate_pct\t%s\n' "$RATE"
} > "$REF/calibration.tsv"

echo "  decoys                : $N_DECOY"
echo "  decoys with any hit   : $N_HITSEQS  (${RATE}%)"
echo "  at E <= $DIAMOND_EVALUE $DIAMOND_SENS"
echo "  wrote $REF/calibration.tsv"

# A handful of chance hits among thousands of shuffled sequences is expected and
# is what an E-value threshold is for. A large fraction is not: it means either
# the threshold is far too permissive or the shuffle did not happen. Either way
# the sweep must not run.
THRESH_PCT=1
BAD=$(python3 -c "print(1 if $RATE > $THRESH_PCT else 0)")
if [ "$BAD" -eq 1 ]; then
  cat >&2 <<MSG

FAIL (exit 3) — ${RATE}% of shuffled decoys match real Oblin-1 at E <= $DIAMOND_EVALUE.
Expected under ${THRESH_PCT}%.

Do NOT respond by tightening DIAMOND_EVALUE. That is what PLAN.md Task 4 Step 4
says, and it is how this control gets inverted (see C-03 at the top of this file):
if the decoys are not genuinely shuffled, no threshold makes them clean except one
that also rejects real Oblins, and the sweep then reports a confident zero.

Diagnose in this order:
  1. head -2 "$REF_DECOY_FAA" — do the ids start with 'decoy_' and does the
     sequence differ from the corresponding record in $REF_OBLIN_FAA?
  2. Is $REF_OBLIN_FAA actually protein? A nucleotide file shuffles into
     something that still aligns to itself.
  3. Are there many near-identical records in $REF_OBLIN_FAA? Shuffles of highly
     repetitive low-complexity sequence can still align. Cluster the reference
     first (cd-hit -c 0.9) and rebuild.
MSG
  exit 3
fi

echo
echo "Databases built and calibration clean."
echo "Next: scripts/03_positive_control.sh   (GATE 1 — do not skip)"
