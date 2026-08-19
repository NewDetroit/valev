#!/usr/bin/env bash
# =============================================================================
# Task 6, offline — enumerate the target niche WITHOUT NCBI.
#
# Reads  : config/niche_query.txt   (the NICHE_* directives, not the Entrez line)
#          s3://logan-pub/stats/logan_accessions_v1.2_SRA2025.csv.zst   11.0 GiB
#          s3://logan-pub/stats/sra_taxid.csv.zst                        424 MiB
# Writes : $WORK/accessions.txt, $WORK/run_metadata.tsv
#          $WORK/logan_accessions_filtered.csv   (cache; --refresh to rebuild)
#
# WHY THIS EXISTS. PLAN.md's only route to an accession list is esearch/efetch.
# Some environments block eutils.ncbi.nlm.nih.gov outright (docs/CORRECTIONS.md
# OPEN-A), which makes the whole pipeline unrunnable there.
#
# Logan publishes its own accession metadata, and it carries every field this
# project needs. Verified 2026-08-19 by decompressing the first 3 MB of each:
#
#   logan_accessions_v1.2_SRA2025.csv.zst
#     acc, assay_type, center_name, consent, experiment, sample_name, instrument,
#     librarylayout, libraryselection, librarysource, platform, sample_acc,
#     biosample, organism, sra_study, releasedate, bioproject, mbytes, ...
#   sra_taxid.csv.zst
#     acc, assay_type, organism, tax_id, taxonomic_rank, scientific_name
#
# `bioproject`, `center_name` and `platform` are all present. That matters more
# than convenience: GATE 3 is built entirely on BioProject independence, so an
# offline path lacking it could find candidates but never validate them.
#
# Filter semantics, matching config/niche_query.txt:
#     (keyword match on organism  OR  tax_id in NICHE_TAXID)
#   AND assay_type     in NICHE_ASSAY_TYPE
#   AND librarysource  in NICHE_LIBRARY_SOURCE   (if any are set)
#
# The 11 GiB file is streamed and filtered, never stored whole. Expect a long
# single pass; the result is cached.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

REFRESH=0
for a in "$@"; do case "$a" in
  --refresh) REFRESH=1 ;;
  -h|--help) sed -n '2,36p' "$0"; exit 0 ;;
  *) echo "unknown argument: $a" >&2; exit 1 ;;
esac; done

QF="$ROOT/config/niche_query.txt"
[ -s "$QF" ] || { echo "missing $QF" >&2; exit 2; }
python3 -c "import zstandard" 2>/dev/null || {
  echo "missing python module 'zstandard' — pip install zstandard" >&2; exit 2; }

CACHE="$WORK/logan_accessions_filtered.csv"
if [ -s "$CACHE" ] && [ "$REFRESH" -eq 0 ]; then
  echo "using cached $CACHE (pass --refresh to rebuild)"
else
  python3 - "$QF" "$CACHE" "$LOGAN_HTTPS" <<'PY'
import csv, io, os, re, sys, urllib.request
import zstandard as zstd

qf, out, base = sys.argv[1:4]

def directives(name):
    vals = []
    for line in open(qf):
        s = line.strip()
        if s.startswith("#") or "=" not in s:
            continue
        k, _, v = s.partition("=")
        if k.strip() != name:
            continue
        v = v.split("#", 1)[0].strip()        # strip trailing comments
        if v:
            vals.append(v)
    return vals

keywords = [k.lower() for k in directives("NICHE_KEYWORD")]
taxids   = set(directives("NICHE_TAXID"))
assays   = {a.lower() for a in directives("NICHE_ASSAY_TYPE")}
sources  = {s.upper() for s in directives("NICHE_LIBRARY_SOURCE")}

if not keywords and not taxids:
    sys.exit("[04b] config/niche_query.txt defines no NICHE_KEYWORD and no NICHE_TAXID "
             "— nothing to select on (exit 2)")
if not assays:
    sys.exit("[04b] no NICHE_ASSAY_TYPE. RNA-Seq is mandatory: Obelisks are RNA elements "
             "and a DNA library cannot contain one (exit 2)")

print(f"keywords       : {keywords or '<none>'}")
print(f"taxids         : {sorted(taxids) or '<none>'}")
print(f"assay_type in  : {sorted(assays)}")
print(f"librarysource  : {sorted(sources) or '<any>'}")

def stream_csv(url, label):
    resp = urllib.request.urlopen(url, timeout=300)
    reader = zstd.ZstdDecompressor().stream_reader(resp)
    text = io.TextIOWrapper(reader, encoding="utf-8", errors="replace", newline="")
    sys.stderr.write(f"streaming {label}\n")
    return csv.DictReader(text)

# Pass 1 — accessions whose tax_id is in the niche. Only matches are retained, so
# memory is bounded by the niche, not by the 38M-row archive.
tax_hits = set()
if taxids:
    n = 0
    for row in stream_csv(f"{base}/stats/sra_taxid.csv.zst", "sra_taxid.csv.zst (424 MiB)"):
        n += 1
        if n % 5_000_000 == 0:
            sys.stderr.write(f"\r  taxid pass: {n:,} rows, {len(tax_hits):,} matched")
            sys.stderr.flush()
        if (row.get("tax_id") or "").strip() in taxids:
            tax_hits.add((row.get("acc") or "").strip())
    sys.stderr.write(f"\r  taxid pass: {n:,} rows, {len(tax_hits):,} matched\n")

# Pass 2 — the accession table, applying the full filter.
kept = seen = 0
with open(out, "w", newline="") as fh:
    w = None
    for row in stream_csv(f"{base}/stats/logan_accessions_v1.2_SRA2025.csv.zst",
                          "logan_accessions_v1.2_SRA2025.csv.zst (11.0 GiB)"):
        seen += 1
        if seen % 2_000_000 == 0:
            sys.stderr.write(f"\r  accession pass: {seen:,} rows, {kept:,} kept")
            sys.stderr.flush()
        acc = (row.get("acc") or "").strip()
        organism = (row.get("organism") or "").lower()
        if not (any(k in organism for k in keywords) or acc in tax_hits):
            continue
        if (row.get("assay_type") or "").strip().lower() not in assays:
            continue
        if sources and (row.get("librarysource") or "").strip().upper() not in sources:
            continue
        if w is None:
            w = csv.DictWriter(fh, fieldnames=row.keys())
            w.writeheader()
        w.writerow(row)
        kept += 1
    sys.stderr.write(f"\r  accession pass: {seen:,} rows, {kept:,} kept\n")

if kept == 0:
    sys.exit("[04b] no accessions matched. Widen NICHE_KEYWORD/NICHE_TAXID, or check "
             "NICHE_ASSAY_TYPE against the assay_type values Logan actually uses (exit 3)")
PY
fi

python3 - <<'PY'
import csv, os, sys
W = os.environ["WORK"]
rows = list(csv.DictReader(open(f"{W}/logan_accessions_filtered.csv", newline="")))
if not rows:
    sys.exit("[04b] cached filter result is empty (exit 3)")

# Logan's column names mapped onto the contract scripts/06 and /11 expect.
MAP = {"Run": "acc", "BioProject": "bioproject", "LibraryStrategy": "assay_type",
       "Platform": "platform", "CenterName": "center_name", "ScientificName": "organism"}
with open(f"{W}/accessions.txt", "w") as fh:
    for r in rows:
        fh.write((r.get("acc") or "").strip() + "\n")
with open(f"{W}/run_metadata.tsv", "w") as fh:
    fh.write("\t".join(MAP) + "\n")
    for r in rows:
        fh.write("\t".join(((r.get(s) or "NA").strip() or "NA") for s in MAP.values()) + "\n")

bps = {(r.get("bioproject") or "").strip() for r in rows} - {"", "NA"}
ctr = {(r.get("center_name") or "").strip() for r in rows} - {"", "NA"}
plt = {(r.get("platform") or "").strip() for r in rows} - {"", "NA"}
print(f"runs        : {len(rows):,}")
print(f"bioprojects : {len(bps):,}")
print(f"centers     : {len(ctr):,}")
print(f"platforms   : {len(plt):,}")
print(f"wrote {W}/accessions.txt and {W}/run_metadata.tsv")

if len(rows) < 200:
    sys.stderr.write(f"\nWARNING: only {len(rows)} runs — niche probably too small; broaden it.\n")
if len(rows) > 200_000:
    sys.stderr.write(f"\nWARNING: {len(rows):,} runs — you will not finish this sweep; narrow it.\n")
if len(bps) < 3:
    sys.stderr.write(f"\nWARNING: only {len(bps)} BioProject(s). Independent replication across\n"
                     "  projects is the entire contamination defence; with fewer than 3 you\n"
                     "  cannot pass GATE 3 no matter what you find.\n")
PY

echo
echo "NOTE: every accession selected here is by construction present in Logan v1.2,"
echo "  so the fast path covers 100% of this list. That is NOT true of the NCBI path,"
echo "  where coverage is whatever predates the freeze — check with estimate_sweep_cost.py."
echo
echo "Next: python3 scripts/estimate_sweep_cost.py"
