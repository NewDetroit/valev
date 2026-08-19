#!/usr/bin/env bash
# =============================================================================
# Re-run every automatable check behind docs/CORRECTIONS.md, so that record can be
# REFRESHED rather than re-derived. Run it before ISEF, and any time the pipeline
# has been idle for a while.
#
# It checks live network facts (they change) and local reference files (they can be
# rebuilt wrong). It cannot check the literature claims in OPEN-C — those need a
# human with journal access, and it says so rather than passing silently.
# =============================================================================
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

PASS=0; FAIL=0; SKIP=0
ok()   { printf '  PASS  %s\n' "$1"; PASS=$((PASS+1)); }
bad()  { printf '  FAIL  %s\n' "$1"; FAIL=$((FAIL+1)); }
skip() { printf '  SKIP  %s\n' "$1"; SKIP=$((SKIP+1)); }

echo "=== obelisk-hunt claim audit — $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo

echo "--- C-01: Logan release string ---"
RM=$(curl -s -m 30 https://raw.githubusercontent.com/IndexThePlanet/Logan/main/README.md || true)
if [ -z "$RM" ]; then
  skip "could not fetch Logan README (network?)"
elif printf '%s' "$RM" | grep -q '## v1\.2 release'; then
  ok "README still advertises v1.2 (config LOGAN_RELEASE=$LOGAN_RELEASE)"
else
  CUR=$(printf '%s' "$RM" | grep -oE '^## v[0-9.]+ release' | head -1)
  bad "README no longer says v1.2 — now '${CUR:-unknown}'. Update config/config.sh and docs/CORRECTIONS.md C-01."
fi

echo
echo "--- C-02: published Obelisk artifacts still present ---"
for k in 01diamond_hmm/Obelisk_paper.hmm \
         01diamond_hmm/Round1_Oblin1_sig.fasta \
         01diamond_hmm/Round1_Oblin2_sig.fasta \
         02_Building_Obelisk_DB/Obelisk_sig_nt.fasta \
         03_Obelisk_DB_QC/Obelisk_db_centroids_nt.fasta \
         03_Obelisk_DB_QC/Obelisk_cen_cen_orf.dmnd; do
  SZ=$(curl -s -I -m 25 "$LOGAN_HTTPS/paper/Obelisk/$k" | awk 'tolower($1)=="content-length:"{print $2}' | tr -d '\r')
  if [ -n "$SZ" ]; then ok "$k ($SZ bytes)"; else bad "$k MISSING — C-02's primary reference path is broken"; fi
done

echo
echo "--- C-07: stats parquet schema ---"
if python3 -c "import pyarrow" 2>/dev/null; then
  python3 - <<'PY' || echo "  FAIL  parquet schema probe errored"
import io, os, sys, urllib.request
import pyarrow.parquet as pq
class H(io.RawIOBase):
    def __init__(s,u):
        s.u,s.p=u,0
        s.size=int(urllib.request.urlopen(urllib.request.Request(u,method="HEAD"),timeout=60).headers["Content-Length"])
    def readable(s): return True
    def seekable(s): return True
    def tell(s): return s.p
    def seek(s,o,w=0):
        s.p = o if w==0 else (s.p+o if w==1 else s.size+o); return s.p
    def read(s,n=-1):
        if n is None or n<0: n=s.size-s.p
        n=min(n,s.size-s.p)
        if n<=0: return b""
        r=urllib.request.Request(s.u,headers={"Range":f"bytes={s.p}-{s.p+n-1}"})
        d=urllib.request.urlopen(r,timeout=120).read(); s.p+=len(d); return d
    def readinto(s,b):
        d=s.read(len(b)); b[:len(d)]=d; return len(d)
url=f"{os.environ['LOGAN_HTTPS']}/stats/{os.environ['LOGAN_STATS_PARQUET']}"
col=os.environ["LOGAN_SIZE_COL"]
pf=pq.ParquetFile(H(url)); names=pf.schema_arrow.names
print(f"  {'PASS' if col in names else 'FAIL'}  {os.environ['LOGAN_STATS_PARQUET']}: "
      f"{pf.metadata.num_rows:,} rows, LOGAN_SIZE_COL='{col}' {'present' if col in names else 'ABSENT'}")
PY
else
  skip "pyarrow not installed — cannot probe the parquet schema"
fi

echo
echo "--- C-14 / C-15: third-party CLIs ---"
for r in Zheludev/VNom rcedgar/circuclust; do
  got=0
  for b in main master; do
    if curl -s -o /dev/null -m 20 -w '%{http_code}' "https://raw.githubusercontent.com/$r/$b/README.md" | grep -q 200; then
      ok "$r README reachable ($b) — re-read it if the invocation stops working"; got=1; break
    fi
  done
  [ "$got" -eq 1 ] || skip "$r README unreachable from here"
done

echo
echo "--- reference files (rebuild with scripts/01 and 02 if these fail) ---"
if command -v seqkit >/dev/null 2>&1; then
  for f in "$REF_OBELISK_NT" "$REF_OBLIN_FAA"; do
    if [ -s "$f" ]; then ok "$(basename "$f"): $(seqkit stats -T "$f" | awk 'NR==2{print $4" seqs, avg len "$7}')"
    else skip "$(basename "$f") not built yet"; fi
  done
  if [ -s "$REF_OBELISK_NT" ]; then
    AVG=$(seqkit stats -T "$REF_OBELISK_NT" | awk 'NR==2{print $7}')
    OKLEN=$(python3 -c "print(1 if 500 <= float('$AVG') <= 2000 else 0)")
    [ "$OKLEN" = 1 ] && ok "obelisk_nt.fna mean length $AVG is in the ~1 kb range Obelisks occupy" \
                     || bad "obelisk_nt.fna mean length $AVG is NOT near 1 kb — likely the wrong column or object"
  fi
else
  skip "seqkit not installed — cannot stat reference files"
fi
if [ -s "$REF/calibration.tsv" ]; then
  RATE=$(awk -F'\t' '$1=="decoy_hit_rate_pct"{print $2}' "$REF/calibration.tsv")
  BADR=$(python3 -c "print(1 if float('${RATE:-100}') > 1 else 0)")
  [ "$BADR" = 0 ] && ok "decoy calibration clean (${RATE}% — C-03)" \
                  || bad "decoy hit rate ${RATE}% — C-03 control is not holding"
else
  skip "no ref/calibration.tsv — run scripts/02_build_search_dbs.sh"
fi

echo
echo "--- OPEN-C: literature figures this script CANNOT verify ---"
cat <<'MSG'
  These need a human with journal access. Do not put any of them on a poster
  until you have read them in the source yourself:
    - 1,744 stringent Obelisks; 788 with >=2 ORFs   (Zheludev et al., Cell 187:6521)
    - ~1,700 Oblin-1 centroid proteins              (Zheludev et al.)
    - 40 marine Obelisks                            (Lopez-Simon et al., ISME J 19:wraf033)
    - mean pLDDT 83.8 for Oblin-1                   (Zheludev et al.)
    - Foldseek best E-value 0.31                    (Urayama et al., Nat Commun 17:3041)
    - 17 SK36 runs across four labs                 (J Mol Evol 93:370)
MSG

echo
echo "======================================================================"
echo " $PASS passed, $FAIL failed, $SKIP skipped"
echo "======================================================================"
[ "$FAIL" -eq 0 ] || exit 1
