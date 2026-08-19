#!/usr/bin/env bash
# =============================================================================
# Task 5 — GATE 1. Recover the known Obelisk-S.s from S. sanguinis SK36 RNA-seq.
#
#   Usage: 03_positive_control.sh [runinfo|triage|deep|confirm|all] [--runs FILE]
#
# THIS IS THE MOST IMPORTANT SCRIPT IN THE PIPELINE. If it cannot re-find an
# Obelisk in data where one is known to exist, every downstream result is
# meaningless -- including a negative one. "We searched N accessions and found
# nothing" is a finding only if the pipeline demonstrably finds Obelisks when
# they are there. Without that, it is indistinguishable from a broken pipeline.
#
# Writes : $WORK/pc/                              working files
#          $RESULTS/positive_control_evidence.tsv the artifact GATE 3 checks for
#
# Exit: 0 passed · 2 missing input · 3 GATE 1 failed
#
# The runinfo stage needs NCBI (docs/CORRECTIONS.md OPEN-A). If eutils is blocked
# for you, supply accessions directly with --runs and start at the triage stage.
# =============================================================================
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

STAGE="all"; RUNS_FILE=""
while [ $# -gt 0 ]; do
  case "$1" in
    runinfo|triage|deep|confirm|all) STAGE="$1" ;;
    --runs) RUNS_FILE="${2:?--runs needs a file}"; shift ;;
    -h|--help) sed -n '2,22p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
  shift
done

PC="$WORK/pc"; mkdir -p "$PC"
RUNS_TSV="$PC/pc_runs.tsv"
SK36_QUERY='"Streptococcus sanguinis"[Organism] AND "RNA-Seq"[Strategy]'

# --------------------------------------------------------------------------
stage_runinfo() {
  if [ -n "$RUNS_FILE" ]; then
    [ -s "$RUNS_FILE" ] || { echo "no such runs file: $RUNS_FILE" >&2; exit 2; }
    # Accept a bare accession list or a Run<TAB>BioProject table.
    if head -1 "$RUNS_FILE" | grep -q $'\t'; then
      cp "$RUNS_FILE" "$RUNS_TSV"
    else
      { printf 'Run\tBioProject\n'; awk 'NF{print $1"\tUNKNOWN"}' "$RUNS_FILE"; } > "$RUNS_TSV"
      echo "NOTE: --runs had no BioProject column. The >=2 BioProject check cannot be"
      echo "  enforced, which weakens this control -- the point of using several projects"
      echo "  is to exercise the same independence logic GATE 3 depends on."
    fi
    echo "using $(( $(wc -l < "$RUNS_TSV") - 1 )) supplied runs"
    return 0
  fi

  command -v esearch >/dev/null 2>&1 || {
    echo "missing esearch and no --runs given." >&2
    echo "  NCBI may be blocked here (OPEN-A). Supply accessions with:" >&2
    echo "    $0 triage --runs my_sk36_runs.txt" >&2
    exit 2; }

  echo "querying SRA: $SK36_QUERY"
  esearch -db sra -query "$SK36_QUERY" | efetch -format runinfo > "$PC/sk36_runinfo.csv"
  [ -s "$PC/sk36_runinfo.csv" ] || { echo "FATAL: empty runinfo" >&2; exit 1; }

  python3 - <<'PY'
import csv, os, sys
from collections import defaultdict
PC = os.path.join(os.environ["WORK"], "pc")
rows = [r for r in csv.DictReader(open(f"{PC}/sk36_runinfo.csv", newline="")) if r.get("Run")]
if not rows:
    sys.exit("[03_positive_control] runinfo had no Run rows (exit 1)")

by_bp = defaultdict(list)
for r in rows:
    by_bp[(r.get("BioProject") or "UNKNOWN").strip()].append(r)

# Round-robin across BioProjects so five runs span as many projects as possible.
# The plan takes .head(3).head(5) of a groupby, which can return five runs from a
# single project and quietly destroys the independence this control is for.
sel, i = [], 0
while len(sel) < 5:
    added = False
    for bp in sorted(by_bp):
        if i < len(by_bp[bp]):
            sel.append((by_bp[bp][i]["Run"], bp)); added = True
            if len(sel) == 5: break
    if not added: break
    i += 1

n_bp = len({bp for _, bp in sel})
with open(f"{PC}/pc_runs.tsv", "w") as fh:
    fh.write("Run\tBioProject\n")
    for run, bp in sel:
        fh.write(f"{run}\t{bp}\n")

print(f"total SK36 RNA-Seq runs : {len(rows)}")
print(f"BioProjects available   : {len(by_bp)}")
print(f"selected                : {len(sel)} runs across {n_bp} BioProject(s)")
for run, bp in sel:
    print(f"   {run}\t{bp}")
if n_bp < 2:
    sys.exit("\n[03_positive_control] FATAL: selection spans <2 BioProjects.\n"
             "  Choosing runs from multiple projects is deliberate — it exercises the same\n"
             "  independence logic GATE 3 relies on. Broaden the query. (exit 3)")
PY
}

# --------------------------------------------------------------------------
stage_triage() {
  [ -s "$RUNS_TSV" ] || { echo "no $RUNS_TSV — run the runinfo stage first" >&2; exit 2; }
  [ -s "$REF_OBLIN_DMND.dmnd" ] || { echo "missing $REF_OBLIN_DMND.dmnd — run scripts/02" >&2; exit 2; }
  command -v diamond >/dev/null 2>&1 || { echo "missing diamond" >&2; exit 2; }

  : > "$PC/triage_summary.tsv"
  printf 'accession\tbioproject\tstatus\thits\tbest_evalue\n' >> "$PC/triage_summary.tsv"
  local any=0
  tail -n +2 "$RUNS_TSV" | while IFS=$'\t' read -r ACC BP; do
    [ -n "$ACC" ] || continue
    OUT="$PC/$ACC.oblin.tsv"
    if [ ! -e "$PC/$ACC.done" ]; then
      URL="$LOGAN_HTTPS/${LOGAN_C#s3://logan-pub/}/$ACC/$ACC.contigs.fa.zst"
      if command -v aws >/dev/null 2>&1; then
        FETCH=(aws s3 cp "$LOGAN_C/$ACC/$ACC.contigs.fa.zst" - --no-sign-request)
      else
        FETCH=(curl -sS --fail --max-time 3600 "$URL")
      fi
      # C-04: qseq_translated, not qseq.
      diamond blastx --db "$REF_OBLIN_DMND" \
        --query <("${FETCH[@]}" 2>/dev/null | zstd -dc 2>/dev/null) \
        --evalue "$DIAMOND_EVALUE" $DIAMOND_SENS --threads "$THREADS" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore qseq_translated \
        --out "$OUT" --quiet 2>>"$LOGS/pc.err" || : > "$OUT"
      : > "$PC/$ACC.done"
    fi
    N=$( [ -s "$OUT" ] && wc -l < "$OUT" || echo 0 )
    BEST=$( [ -s "$OUT" ] && sort -k5,5g "$OUT" | head -1 | cut -f5 || echo "NA" )
    ST=$( [ "$N" -gt 0 ] && echo HIT || echo NONE )
    printf '%s\t%s\t%s\t%s\t%s\n' "$ACC" "$BP" "$ST" "$N" "$BEST" >> "$PC/triage_summary.tsv"
    echo "  $ACC ($BP): $ST $N hits, best E=$BEST"
  done

  any=$(awk -F'\t' 'NR>1 && $3=="HIT"' "$PC/triage_summary.tsv" | wc -l)
  echo
  echo "accessions with hits: $any"
  if [ "$any" -eq 0 ]; then
    cat >&2 <<MSG

GATE 1 FAILED (exit 3) — zero Oblin hits in SK36 RNA-seq, where Obelisk-S.s is known
to be present. Do NOT continue: a silent zero here becomes a false "we found nothing
novel" later. Debug in this order:

  (a) Is $REF_OBLIN_FAA real protein?
        head -2 "$REF_OBLIN_FAA"
      If it is nucleotide, scripts/01 selected the wrong object.
  (b) Do these accessions exist in Logan $LOGAN_RELEASE?
        curl -sI "$LOGAN_HTTPS/${LOGAN_C#s3://logan-pub/}/<ACC>/<ACC>.contigs.fa.zst"
  (c) Retry at the relaxed threshold:
        DIAMOND_EVALUE=$DIAMOND_EVALUE_RELAXED $0 triage
  (d) The runs may postdate the freeze. NOTE: under correction C-01 the freeze is
      $LOGAN_FREEZE, not the December 2023 the original plan assumed, so this is a
      far less likely explanation than the plan suggests. Check the release date of
      the runs before believing it.
MSG
    exit 3
  fi
}

# --------------------------------------------------------------------------
stage_deep() {
  [ -s "$PC/triage_summary.tsv" ] || { echo "run the triage stage first" >&2; exit 2; }
  # Resolve the accession programmatically. The plan leaves ACC=<accession that
  # produced hits> as a manual fill-in.
  ACC=$(awk -F'\t' 'NR>1 && $3=="HIT"' "$PC/triage_summary.tsv" \
        | sort -t$'\t' -k5,5g | head -1 | cut -f1)
  [ -n "$ACC" ] || { echo "no hit accession found" >&2; exit 3; }
  echo "deep path on best-scoring accession: $ACC"
  bash "$SCRIPTS/07_deep_assemble.sh" "$ACC" || { echo "assembly failed for $ACC" >&2; exit 1; }

  # VNom's real interface: bare stem, no -o, output to 4_final_clusters in CWD,
  # filename needing exactly one underscore. See C-15 in scripts/08_run_vnom.sh.
  D="$PC/vnom/$ACC"; mkdir -p "$D"
  STEM="${ACC}_contigs"
  seqkit grep -v -s -p 'N' "$WORK/deep/$ACC/transcripts.fasta" \
    | sed "s/NODE/${ACC}/g" > "$D/${STEM}.fasta"
  ( cd "$D" && python "${VNOM_PY:-$VNOM/VNom.py}" -i "$STEM" \
      -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1 \
      > "${STEM}_VNom.log" 2>&1 ) || { echo "VNom failed on $ACC — see $D/${STEM}_VNom.log" >&2; exit 1; }
  [ -d "$D/4_final_clusters" ] || { echo "VNom nominated nothing for $ACC (no 4_final_clusters)" >&2; exit 3; }
  echo "$ACC" > "$PC/pc_accession.txt"
}

# --------------------------------------------------------------------------
stage_confirm() {
  ACC=$(cat "$PC/pc_accession.txt" 2>/dev/null || true)
  [ -n "$ACC" ] || { echo "run the deep stage first" >&2; exit 2; }
  for t in makeblastdb blastn; do
    command -v "$t" >/dev/null 2>&1 || { echo "missing $t" >&2; exit 2; }
  done
  [ -s "$REF_OBELISK_NT" ] || { echo "missing $REF_OBELISK_NT — run scripts/01" >&2; exit 2; }

  makeblastdb -in "$REF_OBELISK_NT" -dbtype nucl -out "$PC/obelisk_db" >/dev/null

  # Glob every VNom output rather than assuming one filename, as the plan does
  # with <vnom_output>.fasta.
  shopt -s nullglob
  QUERIES=( "$PC/vnom/$ACC"/4_final_clusters/*.fasta "$PC/vnom/$ACC"/4_final_clusters/*.fa "$PC/vnom/$ACC"/4_final_clusters/*.fna )
  shopt -u nullglob
  [ "${#QUERIES[@]}" -gt 0 ] || { echo "no VNom FASTA under $PC/vnom/$ACC" >&2; exit 3; }

  EV="$RESULTS/positive_control_evidence.tsv"
  { printf '# GATE 1 — recovery of a known Obelisk from S. sanguinis SK36\n'
    printf '# accession\t%s\n' "$ACC"
    printf '# logan_release\t%s (freeze %s)\n' "$LOGAN_RELEASE" "$LOGAN_FREEZE"
    printf '# generated\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'qseqid\tsseqid\tpident\tlength\tevalue\tbitscore\n'; } > "$EV"

  for q in "${QUERIES[@]}"; do
    blastn -query "$q" -db "$PC/obelisk_db" \
           -outfmt "6 qseqid sseqid pident length evalue bitscore" \
           -evalue 1e-20 -num_threads "$THREADS" 2>/dev/null
  done | sort -t$'\t' -k5,5g >> "$EV"

  BEST_ID=$(awk -F'\t' '!/^#/ && $3 ~ /^[0-9.]+$/ {print $3}' "$EV" | sort -rn | head -1)
  BEST_ID=${BEST_ID:-0}
  echo
  echo "best identity to a known Obelisk: ${BEST_ID}%"
  echo "evidence: $EV"

  PASS=$(python3 -c "print(1 if float('${BEST_ID}') >= 95 else 0)")
  echo
  echo "======================================================================"
  if [ "$PASS" -eq 1 ]; then
    echo "GATE 1 PASSED — recovered a known Obelisk at ${BEST_ID}% identity."
    echo "  The pipeline demonstrably finds Obelisks in data where they exist."
    echo "  Everything downstream, including a negative result, is now interpretable."
    echo "  Keep $EV: it is a poster figure and the most persuasive single artifact"
    echo "  you can show a skeptical judge."
    echo "======================================================================"
    return 0
  fi
  echo "GATE 1 FAILED (exit 3) — best identity ${BEST_ID}%, need >=95%."
  echo "  VNom produced circular candidates but none match a known Obelisk."
  echo "  Check that $REF_OBELISK_NT actually contains Obelisk-S.s, and that VNom's"
  echo "  invocation is correct (its CLI is unverified — see docs/CORRECTIONS.md OPEN-C)."
  echo "======================================================================"
  return 3
}

case "$STAGE" in
  runinfo) stage_runinfo ;;
  triage)  stage_triage ;;
  deep)    stage_deep ;;
  confirm) stage_confirm ;;
  all)     stage_runinfo && stage_triage && stage_deep && stage_confirm ;;
esac
