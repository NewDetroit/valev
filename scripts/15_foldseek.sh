#!/usr/bin/env bash
# =============================================================================
# Task 14 Steps 3-4 — Foldseek structural homology, then the phylogeny.
#
#   Usage: 15_foldseek.sh [foldseek|tree|all]
#
# PLAN.md does Foldseek by hand on the web server. Scripted here when the CLI is
# available, with the manual route documented as the fallback rather than the
# default.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"
STAGE="${1:-all}"

CFD="$RESULTS/figures/colabfold"

do_foldseek() {
  if ! command -v foldseek >/dev/null 2>&1; then
    cat <<MSG
foldseek CLI not installed — use the web server instead:

  1. https://search.foldseek.com/search
  2. Upload each PDB from $CFD/
  3. Search against AFDB50, AFDB-proteome and PDB
  4. Record the top 5 hits and E-values per model in:
       $CFD/foldseek.tsv

Install the CLI with:  conda install -c conda-forge -c bioconda foldseek
MSG
    return 0
  fi
  shopt -s nullglob
  PDBS=( "$CFD"/*.pdb )
  shopt -u nullglob
  [ "${#PDBS[@]}" -gt 0 ] || { echo "no PDBs in $CFD — run ColabFold first (exit 2)" >&2; exit 2; }
  mkdir -p "$WORK/foldseek"
  : > "$RESULTS/foldseek_hits.tsv"
  printf 'query\ttarget\tfident\talnlen\tevalue\tbits\tprob\n' >> "$RESULTS/foldseek_hits.tsv"
  # Must word-split into separate database names; a quoted expansion iterates once
  # over the single joined word "afdb50 pdb" and foldseek then fails.
  read -r -a DBS <<< "${FOLDSEEK_DBS:-afdb50 pdb}"
  for db in "${DBS[@]}"; do
    # foldseek writes a file prefix, not a directory, so test the prefix itself.
    [ -s "$WORK/foldseek/$db" ] || foldseek databases "$db" "$WORK/foldseek/$db" "$WORK/foldseek/tmp"
    foldseek easy-search "${PDBS[@]}" "$WORK/foldseek/$db" \
      "$WORK/foldseek/hits_$db.tsv" "$WORK/foldseek/tmp" \
      --format-output "query,target,fident,alnlen,evalue,bits,prob" \
      -e 10 --max-seqs 50 >/dev/null
    cat "$WORK/foldseek/hits_$db.tsv" >> "$RESULTS/foldseek_hits.tsv"
  done
  echo "wrote $RESULTS/foldseek_hits.tsv"
  echo "  Both outcomes are good, and both need the numbers recorded:"
  echo "    no significant hits -> argues a genuinely novel fold"
  echo "    hits to known Oblin-1 -> confirms a true Obelisk relative"
}

do_tree() {
  for t in mafft iqtree; do
    command -v "$t" >/dev/null 2>&1 || { echo "missing $t — conda activate obelisk" >&2; exit 2; }
  done
  [ -s "$RESULTS/validated_oblins.faa" ] || { echo "missing $RESULTS/validated_oblins.faa — run scripts/14" >&2; exit 2; }
  [ -s "$REF_OBLIN_FAA" ] || { echo "missing $REF_OBLIN_FAA" >&2; exit 2; }
  cat "$RESULTS/validated_oblins.faa" "$REF_OBLIN_FAA" > "$WORK/tree_in.faa"
  mafft --auto --thread "$THREADS" "$WORK/tree_in.faa" > "$WORK/tree_aln.faa"
  iqtree -s "$WORK/tree_aln.faa" -m MFP -bb 1000 -nt "$THREADS" \
         -pre "$RESULTS/figures/oblin_tree" -redo
  echo "wrote $RESULTS/figures/oblin_tree.treefile"
  echo "  Where your sequences fall is the point. A well-supported clade separate from"
  echo "  all published Oblins is the phylogenetic version of 'this is new'. Report the"
  echo "  bootstrap value for that node — an unsupported clade claims nothing."
}

case "$STAGE" in
  foldseek) do_foldseek ;;
  tree)     do_tree ;;
  all)      do_foldseek; do_tree ;;
  *) echo "usage: $(basename "$0") [foldseek|tree|all]" >&2; exit 1 ;;
esac
