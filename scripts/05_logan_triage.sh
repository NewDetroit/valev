#!/usr/bin/env bash
# =============================================================================
# Task 7 — per-accession triage worker.  Usage: 05_logan_triage.sh <ACCESSION>
#
# Streams one accession's Logan contigs out of S3, through zstd, into DIAMOND,
# and writes any Oblin-like hits. Nothing touches disk: the contigs exist only
# in the pipe. This is what makes a sweep over hundreds of thousands of
# accessions affordable.
#
# Writes : $WORK/triage/<ACC>.tsv     DIAMOND hits (may be empty — that is a result)
#          $WORK/triage/<ACC>.done    completion sentinel
#          $WORK/triage/<ACC>.status  one word: HIT | NONE | ABSENT | FAIL
#
# `set -e` is deliberately NOT used. This runs once per accession inside a batch
# loop; one bad accession must not kill the other 99 in the array task.
# ---------------------------------------------------------------------------
# CORRECTIONS IMPLEMENTED HERE
#
# C-04  PLAN.md requests `--outfmt 6 ... qseq`. In blastx the query is nucleotide
#       and DIAMOND defines `qseq` as the aligned part of the QUERY, so it emits
#       DNA. Task 8 then writes that to a file named .faa and runs hmmsearch
#       against an amino-acid profile. Confirmed against DIAMOND's source
#       (src/output/blast_tab_format.cpp): line 92 defines `qseq_translated` as
#       "Aligned part of query sequence (translated)". That is the field we want.
#
# C-05  PLAN.md guards with `[ -s "$OUT" ] && exit 0`, but records no-hit
#       accessions with `touch`, so they are empty and `-s` is false. On re-run
#       every no-hit accession — the overwhelming majority — is re-downloaded and
#       re-searched. The guard skips almost nothing. Here the guard is the
#       explicit .done sentinel, so "ran, found nothing" is distinguishable from
#       "never ran".
#
# C-06  PLAN.md appends every hit to a single shared $RESULTS/triage_hits.tsv.
#       O_APPEND is atomic only up to PIPE_BUF (4096 bytes); hit rows carrying a
#       translated sequence exceed that, and --array=...%100 means 100 concurrent
#       writers interleaving partial lines into the file Task 8 then parses with
#       fixed column names. This worker writes ONLY its own file. scripts/06
#       merges them.
#
# Missing accessions are normal — not every SRA run is in Logan. PLAN.md sends
# all failure modes to /dev/null, which conflates "absent from Logan" with
# "download failed" and makes the final coverage statistic meaningless. They are
# recorded separately here.
# =============================================================================
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

ACC="${1:-}"
[ -n "$ACC" ] || { echo "usage: $(basename "$0") <ACCESSION>" >&2; exit 1; }

TDIR="$WORK/triage"
mkdir -p "$TDIR"
OUT="$TDIR/$ACC.tsv"
DONE="$TDIR/$ACC.done"
STATUS="$TDIR/$ACC.status"

# C-05: existence, not non-emptiness.
if [ -e "$DONE" ]; then
  echo "$ACC skip ($(cat "$STATUS" 2>/dev/null || echo done))"
  exit 0
fi

[ -s "$REF_OBLIN_DMND.dmnd" ] || {
  echo "missing $REF_OBLIN_DMND.dmnd — run scripts/02_build_search_dbs.sh" >&2; exit 2; }

URL_S3="$LOGAN_C/$ACC/$ACC.contigs.fa.zst"
URL_HTTP="$LOGAN_HTTPS/${LOGAN_C#s3://logan-pub/}/$ACC/$ACC.contigs.fa.zst"

# Probe before spending a transfer. An absent accession is a data fact worth
# recording, not an error worth retrying.
if command -v aws >/dev/null 2>&1; then
  if ! aws s3 ls "$URL_S3" --no-sign-request >/dev/null 2>&1; then
    : > "$OUT"; echo ABSENT > "$STATUS"; : > "$DONE"
    echo "$ACC ABSENT (not in Logan $LOGAN_RELEASE)"; exit 0
  fi
  FETCH=(aws s3 cp "$URL_S3" - --no-sign-request)
else
  CODE=$(curl -s -o /dev/null -I -m 60 -w '%{http_code}' "$URL_HTTP" 2>/dev/null || echo 000)
  if [ "$CODE" != "200" ]; then
    : > "$OUT"; echo ABSENT > "$STATUS"; : > "$DONE"
    echo "$ACC ABSENT (HTTP $CODE, not in Logan $LOGAN_RELEASE)"; exit 0
  fi
  FETCH=(curl -sS --fail --max-time 3600 "$URL_HTTP")
fi

# DIAMOND reads its query from a path. Process substitution hands it /dev/fd/N,
# so the stream still never lands on disk. Set TRIAGE_SPOOL=1 if a DIAMOND build
# refuses a non-seekable query — that spools to node-local scratch instead, which
# costs disk but is otherwise identical.
run_diamond() {
  diamond blastx \
    --db "$REF_OBLIN_DMND" --query "$1" \
    --evalue "$DIAMOND_EVALUE" $DIAMOND_SENS --threads "$TRIAGE_THREADS" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore qseq_translated \
    --out "$OUT" --quiet
}

rc=0
if [ "${TRIAGE_SPOOL:-0}" = "1" ]; then
  TMP="$(mktemp "${TMPDIR:-/tmp}/${ACC}.XXXXXX.fa")"
  trap 'rm -f "$TMP"' EXIT
  if "${FETCH[@]}" 2>>"$LOGS/triage.err" | zstd -dc > "$TMP" 2>>"$LOGS/triage.err"; then
    run_diamond "$TMP" 2>>"$LOGS/triage.err" || rc=$?
  else
    rc=1
  fi
else
  run_diamond <("${FETCH[@]}" 2>>"$LOGS/triage.err" | zstd -dc 2>>"$LOGS/triage.err") \
    2>>"$LOGS/triage.err" || rc=$?
fi

if [ "$rc" -ne 0 ]; then
  # Leave no .done: a genuine failure must be retried on the next pass, unlike
  # a legitimately empty result.
  echo FAIL > "$STATUS"
  echo "$ACC FAIL (rc=$rc — see $LOGS/triage.err)"
  exit 1
fi

if [ -s "$OUT" ]; then
  echo HIT > "$STATUS"; : > "$DONE"
  echo "$ACC HIT $(wc -l < "$OUT")"
else
  : > "$OUT"; echo NONE > "$STATUS"; : > "$DONE"
  echo "$ACC NONE"
fi
