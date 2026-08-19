#!/usr/bin/env bash
# =============================================================================
# Task 13 — RNA secondary structure of the validated candidates.
#
#   Usage: 13_rnafold.sh [--rfam]     (--rfam downloads/searches Rfam; ~1.5 GB)
#
# Obelisks are defined partly by a rod-like fold, so this is corroborating
# evidence, not decoration. PLAN.md says to "look for the rod-like fold"; looking
# is not a measurement, so this also computes paired fraction and MFE per
# nucleotide and writes them to a table you can defend.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

DO_RFAM=0
for a in "$@"; do case "$a" in
  --rfam) DO_RFAM=1 ;;
  -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
  *) echo "unknown argument: $a" >&2; exit 1 ;;
esac; done

command -v RNAfold >/dev/null 2>&1 || { echo "missing RNAfold (viennarna) — conda activate obelisk" >&2; exit 2; }
[ -s "$RESULTS/validated.fna" ] || { echo "missing $RESULTS/validated.fna — run scripts/11 first" >&2; exit 2; }

FD="$RESULTS/figures/rnafold"; mkdir -p "$FD"
( cd "$FD" && RNAfold --noPS < "$RESULTS/validated.fna" > validated_fold.txt )
echo "wrote $FD/validated_fold.txt"

python3 - <<'PY'
import os, re
R = os.environ["RESULTS"]
src = f"{R}/figures/rnafold/validated_fold.txt"
recs, cur = [], None
for line in open(src):
    line = line.rstrip("\n")
    if line.startswith(">"):
        cur = {"id": line[1:].split()[0], "seq": "", "db": "", "mfe": None}
        recs.append(cur)
    elif cur is not None and not cur["seq"]:
        cur["seq"] = line.strip()
    elif cur is not None and ("(" in line or "." in line):
        m = re.match(r"^([.()\[\]{}<>]+)\s*\(\s*(-?\d+\.\d+)\)", line.strip())
        if m:
            cur["db"], cur["mfe"] = m.group(1), float(m.group(2))

out = f"{R}/rna_structure.tsv"
with open(out, "w") as fh:
    fh.write("id\tlength\tmfe\tmfe_per_nt\tpaired_fraction\tlongest_helix\n")
    for r in recs:
        if not r["db"]:
            continue
        n = len(r["db"])
        paired = sum(1 for c in r["db"] if c in "()")
        longest = max((len(m) for m in re.findall(r"\(+|\)+", r["db"])), default=0)
        fh.write(f"{r['id']}\t{n}\t{r['mfe']}\t{r['mfe']/n:.4f}\t{paired/n:.4f}\t{longest}\n")
print(f"wrote {out} ({len(recs)} structures)")
print("  A rod-like Obelisk fold shows a HIGH paired fraction (typically >0.6) and one")
print("  dominant long helix rather than many short ones. Report both numbers rather")
print("  than asserting the morphology from a picture.")
PY

if [ "$DO_RFAM" -eq 1 ]; then
  command -v cmscan >/dev/null 2>&1 || { echo "missing cmscan (infernal)" >&2; exit 2; }
  CM="$REF/Rfam.cm"
  if [ ! -s "$CM.i1f" ]; then
    echo "downloading Rfam.cm (~1.5 GB, cached in $REF)"
    wget -q -O "$CM.gz" https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
    gunzip -f "$CM.gz"; cmpress -F "$CM"
  fi
  # -Z is the search-space size in Mb (2 x total residues, both strands). Without
  # it, cmscan's E-values are not comparable between runs.
  BP=$(seqkit stats -T "$RESULTS/validated.fna" | awk 'NR==2{print $5}')
  Z=$(python3 -c "print(f'{2*$BP/1e6:.6f}')")
  cmscan --tblout "$RESULTS/validated_rfam.tbl" --cut_ga -Z "$Z" --cpu "$THREADS" \
         "$CM" "$RESULTS/validated.fna" > /dev/null
  echo "wrote $RESULTS/validated_rfam.tbl (-Z $Z)"
  grep -v '^#' "$RESULTS/validated_rfam.tbl" | awk '{print $2, $3, $16}' | sort -u | head
  echo "  Hammerhead ribozyme hits are supporting evidence. Absence is NOT disqualifying —"
  echo "  many Obelisks carry no recognisable ribozyme."
fi
