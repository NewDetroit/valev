#!/usr/bin/env bash
# =============================================================================
# Task 6 — enumerate the target niche via NCBI Entrez.
#
# Reads  : config/niche_query.txt
# Writes : $WORK/niche_runinfo.csv, $WORK/accessions.txt, $WORK/run_metadata.tsv
#
# REQUIRES NETWORK ACCESS TO eutils.ncbi.nlm.nih.gov. Some environments block it
# (see docs/CORRECTIONS.md OPEN-A). If yours does, use the offline equivalent:
#     scripts/04b_select_accessions_offline.sh
# which needs no NCBI and produces the same two output files.
#
# CORRECTION C-13: PLAN.md reads the BioProject column with `cut -d, -f22`.
# SRA runinfo column order is not contractual and position 22 silently yields the
# wrong field if NCBI reorders. Everything here selects by header name.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

QF="$ROOT/config/niche_query.txt"
[ -s "$QF" ] || { echo "missing $QF" >&2; exit 2; }
command -v esearch >/dev/null 2>&1 || {
  echo "missing esearch (entrez-direct) — conda activate obelisk" >&2
  echo "  If NCBI is unreachable here, use scripts/04b_select_accessions_offline.sh" >&2
  exit 2; }

# Exactly one active query line: comment lines, blanks and NICHE_* directives out.
QUERY=$(grep -vE '^[[:space:]]*(#|$)' "$QF" | grep -vE '^[[:space:]]*NICHE_[A-Z_]+=' || true)
N_Q=$(printf '%s\n' "$QUERY" | grep -cve '^[[:space:]]*$' || true)
if [ "$N_Q" -ne 1 ]; then
  echo "FATAL: $QF has $N_Q active query lines, expected exactly 1." >&2
  echo "  The usual mistake is uncommenting the alternative query without commenting" >&2
  echo "  the default, which silently sweeps a niche you did not choose." >&2
  exit 2
fi
echo "query: $QUERY"

esearch -db sra -query "$QUERY" | efetch -format runinfo > "$WORK/niche_runinfo.csv"
[ -s "$WORK/niche_runinfo.csv" ] || { echo "FATAL: empty runinfo from NCBI" >&2; exit 1; }

python3 - <<'PY'
import csv, os, sys
W = os.environ["WORK"]
WANT = ["Run", "BioProject", "LibraryStrategy", "Platform", "CenterName", "ScientificName"]
rows, hdr = [], None
with open(f"{W}/niche_runinfo.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        if not r.get("Run"):          # runinfo repeats its header between pages
            continue
        if hdr is None:
            hdr = list(r.keys())
        rows.append(r)
if not rows:
    sys.exit("[04_select_accessions] runinfo contained no Run rows (exit 1)")

missing = [c for c in WANT if c not in hdr]
if missing:
    sys.stderr.write(f"[04_select_accessions] WARNING: runinfo lacks {missing}; "
                     f"those columns will be NA. Present: {hdr}\n")

with open(f"{W}/accessions.txt", "w") as fh:
    for r in rows:
        fh.write(r["Run"].strip() + "\n")
with open(f"{W}/run_metadata.tsv", "w") as fh:
    fh.write("\t".join(WANT) + "\n")
    for r in rows:
        fh.write("\t".join((r.get(c) or "NA").strip() or "NA" for c in WANT) + "\n")

bps = {(r.get("BioProject") or "").strip() for r in rows} - {"", "NA"}
print(f"runs        : {len(rows):,}")
print(f"bioprojects : {len(bps):,}")

# PLAN.md's sanity thresholds, as checks rather than prose.
if len(rows) < 200:
    sys.stderr.write(f"\nWARNING: only {len(rows)} runs — the niche is probably too small.\n"
                     "  Broaden config/niche_query.txt.\n")
if len(rows) > 200_000:
    sys.stderr.write(f"\nWARNING: {len(rows):,} runs — you will not finish this sweep.\n"
                     "  Narrow config/niche_query.txt.\n")
if len(bps) < 3:
    sys.stderr.write(f"\nWARNING: only {len(bps)} BioProject(s). Independent replication across\n"
                     "  projects is the entire contamination defence; with fewer than 3 you\n"
                     "  cannot pass GATE 3 no matter what you find.\n")
PY

echo
echo "Next: python3 scripts/estimate_sweep_cost.py   (size the sweep before running it)"
