#!/usr/bin/env bash
# =============================================================================
# 00_verify_logan.sh — Phase 0 gate.
#
# Proves Logan is reachable and that config.sh's hardcoded ground truth still
# holds, before any science runs. Refuses to proceed (non-zero exit) if it
# doesn't. This is the plan's Task 1, made real: every "Expected: ..." line
# in the original checklist becomes a live check here, not a thing a human
# eyeballs once and forgets to redo.
#
# Ground truth this script defends, verified live 2026-08-19
# (docs/CORRECTIONS.md C-01, C-02, C-07 — read those before touching this):
#
#   - Logan's current release is v1.2 (38,124,741 accessions at s3://logan-pub/c/).
#     v1.0 is parked at c1.0/ (27,269,310 accessions) and may be deleted —
#     Stats-v1.1.md said "will be deleted in ~1 year" and that window has
#     already passed. config.sh defaults LOGAN_RELEASE=v1.2 accordingly.
#     This script re-checks the *README's* claim independent of whatever
#     LOGAN_RELEASE is currently exported, because it is guarding config.sh's
#     DEFAULT (an operator may deliberately run with v1.0 for replication —
#     see docs/RUNBOOK.md — and that is not itself a staleness signal).
#   - The smoke-test accession DRR000016 was confirmed present under BOTH
#     c/ (189,571 bytes) and c1.0/ (190,293 bytes) via `curl -I` before being
#     hardcoded below. Contrary to the plan's claim, it is not the accession
#     used in Logan's own Accessions.md tutorial (that uses DRR000273) — it
#     is simply confirmed to exist, which is all the smoke test needs.
#
# Backend: prefers the `aws` CLI (what the rest of the pipeline uses on the
# cluster, per the plan) and falls back to anonymous HTTPS + curl when aws is
# absent — true of this repo's dev/CI container, and possibly true of a bare
# login node before env/environment.yml's conda env exists. Every check here
# is a read-only HEAD/list; nothing needs credentials (`--no-sign-request`).
#
# Exit codes (docs/../INTERFACES.md convention):
#   0  all checks passed — safe to proceed to Phase 1.
#   1  unexpected/infrastructure failure (network unreachable, an artifact
#      that should exist does not).
#   2  not used by this script — Phase 0 has no earlier script to depend on.
#   3  the C-01 guard tripped: Logan's advertised release has moved off v1.2.
#      This is a real result, not a crash — config.sh must be updated before
#      anything downstream is trustworthy.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

DOWNLOAD_STATS=0
for arg in "$@"; do
  case "$arg" in
    --download)
      DOWNLOAD_STATS=1
      ;;
    -h|--help)
      echo "Usage: $(basename "$0") [--download]"
      echo "  --download   also fetch \$LOGAN_STATS_PARQUET into \$REF (~1.7 GB)."
      echo "               Off by default — this script only HEADs it by default."
      exit 0
      ;;
    *)
      echo "unknown argument: $arg (use --help)" >&2
      exit 1
      ;;
  esac
done

# --- backend selection -------------------------------------------------------
BACKEND="curl"
command -v aws >/dev/null 2>&1 && BACKEND="aws"
LOGAN_BUCKET_NAME="${LOGAN_BUCKET#s3://}"

# --- bookkeeping --------------------------------------------------------------
# Every check is wrapped in an `if`, never run bare, so one network hiccup
# records a FAIL line and lets the rest of the checks still run — a partial
# report is far more useful here than an abort on the first failure.
SUMMARY=()
FAILS=0
RELEASE_STALE=0

pass() { SUMMARY+=("PASS  $1"); echo "PASS  $1"; }
fail() { SUMMARY+=("FAIL  $1"); echo "FAIL  $1" >&2; FAILS=$((FAILS + 1)); }
info() { echo "INFO  $1"; }

human_size() {
  awk -v b="$1" 'BEGIN {
    split("B K M G T P", u, " "); v = b; i = 1
    while (v >= 1024 && i < 6) { v /= 1024; i++ }
    printf "%.1f%s", v, u[i]
  }'
}

# s3_head_bytes KEY — print Content-Length for a bucket-relative key on
# stdout and return 0, or return 1 (nothing printed) if the object does not
# exist or the request failed. KEY has no leading slash.
s3_head_bytes() {
  local key="$1" bytes=""
  if [ "$BACKEND" = aws ]; then
    bytes="$(aws s3api head-object --bucket "$LOGAN_BUCKET_NAME" --key "$key" \
              --no-sign-request --query 'ContentLength' --output text 2>/dev/null || true)"
  else
    bytes="$(curl -fsSI "$LOGAN_HTTPS/$key" 2>/dev/null | tr -d '\r' \
              | awk -F': ' 'tolower($1) == "content-length" { print $2 }' || true)"
  fi
  case "$bytes" in
    ''|*[!0-9]*) return 1 ;;
    *) printf '%s' "$bytes" ;;
  esac
}

# s3_list_root — one top-level bucket entry per line; directories keep their
# trailing '/' (matching S3 CommonPrefixes convention) so entries are
# unambiguous against files of the same basename.
s3_list_root() {
  if [ "$BACKEND" = aws ]; then
    aws s3 ls "$LOGAN_BUCKET/" --no-sign-request 2>/dev/null | awk '{print $NF}'
  else
    curl -fsS "$LOGAN_HTTPS/?list-type=2&delimiter=/" 2>/dev/null \
      | grep -oE '<(Prefix|Key)>[^<]*</(Prefix|Key)>' \
      | sed -E 's#</?(Prefix|Key)>##g' \
      | grep -v '^$'
  fi
}

echo "=============================================================================="
echo " obelisk-hunt Phase 0 — Logan + toolchain verification"
echo " LOGAN_RELEASE=$LOGAN_RELEASE   LOGAN_C=$LOGAN_C   backend=$BACKEND"
echo "=============================================================================="

# --- 0. aws CLI presence — informational only. ------------------------------
# The plan's Task 1 Step 2 treated a missing aws CLI as a hard stop before
# continuing. This script does not: every check it performs is an anonymous
# read, which curl serves identically. The download-heavy scripts later in
# the pipeline still shell out to `aws s3 cp ... --no-sign-request` exactly
# as the plan specifies, so a real cluster deploy still needs awscli
# (env/environment.yml installs it) — this script just doesn't require it
# to verify itself.
if [ "$BACKEND" = aws ]; then
  pass "aws CLI found: $(aws --version 2>&1)"
else
  info "aws CLI not found on PATH — using curl against \$LOGAN_HTTPS instead."
  info "Install awscli (env/environment.yml) before running download-heavy scripts."
fi

# --- 1. Bucket root listing: report every top-level prefix, programmatically.
# The plan says "record every prefix you see" — this discovers them live
# rather than checking a hardcoded list, so a genuinely new prefix (e.g. a
# future v2 release directory) is surfaced instead of silently ignored.
echo
echo "--- bucket root: $LOGAN_BUCKET/ ---"
KNOWN_BASELINE=" c/ c1.0/ p/ paper/ stats/ u/ index.html "
mapfile -t ROOT_ENTRIES < <(s3_list_root | sort -u)
if [ "${#ROOT_ENTRIES[@]}" -eq 0 ]; then
  fail "could not list $LOGAN_BUCKET/ (network unreachable or listing blocked)"
else
  for e in "${ROOT_ENTRIES[@]}"; do
    case "$KNOWN_BASELINE" in
      *" $e "*) echo "      $e" ;;
      *)        echo "      $e   <-- not in the 2026-08-19 known baseline (c/ c1.0/ p/ paper/ stats/ u/ index.html); investigate before relying on it" ;;
    esac
  done
  pass "listed ${#ROOT_ENTRIES[@]} top-level entries: ${ROOT_ENTRIES[*]}"

  c_rel="${LOGAN_C#"$LOGAN_BUCKET"/}/"
  if printf '%s\n' "${ROOT_ENTRIES[@]}" | grep -qx "$c_rel"; then
    pass "configured LOGAN_C prefix '$c_rel' is present at bucket root"
  else
    fail "configured LOGAN_C prefix '$c_rel' NOT found at bucket root — LOGAN_RELEASE=$LOGAN_RELEASE is misconfigured, or that release was removed"
  fi
fi

# --- 2. C-01 guard: Logan's advertised current release must still be v1.2. -
echo
echo "--- Logan README release check (automated C-01 guard) ---"
README_URL="https://raw.githubusercontent.com/IndexThePlanet/Logan/main/README.md"
if readme="$(curl -fsS "$README_URL" 2>/dev/null)"; then
  # Exactly one "## vX.Y release" heading is expected today ("## v1.2
  # release"); take the max in case a future README keeps old headings
  # around alongside a new one.
  latest="$(printf '%s\n' "$readme" \
    | grep -oE '^## v[0-9]+\.[0-9]+ release' \
    | grep -oE '[0-9]+\.[0-9]+' \
    | sort -t. -k1,1n -k2,2n \
    | tail -1)"
  if [ -z "$latest" ]; then
    fail "no '## vX.Y release' heading found in Logan README — page format changed; the C-01 guard cannot confirm the release and must be updated by hand"
    RELEASE_STALE=1
  elif [ "$latest" != "1.2" ]; then
    fail "Logan README's newest release heading is v$latest, not v1.2 — config.sh's C-01 default has gone stale. Update LOGAN_RELEASE's default and the v1.2 case in config/config.sh (LOGAN_C, LOGAN_STATS_PARQUET, LOGAN_SIZE_COL, LOGAN_N_ACCESSIONS), log a new docs/CORRECTIONS.md entry, then re-run this script."
    RELEASE_STALE=1
  else
    pass "Logan README confirms v1.2 is the current release ('## v1.2 release' heading found)"
  fi
else
  fail "could not fetch $README_URL (network unreachable, or the repo/branch moved)"
fi

# --- 3. Smoke-test one known-present accession against $LOGAN_C. -----------
echo
echo "--- smoke-test accession ---"
SMOKE_ACC="DRR000016"
c_rel_bare="${LOGAN_C#"$LOGAN_BUCKET"/}"
smoke_key="$c_rel_bare/$SMOKE_ACC/$SMOKE_ACC.contigs.fa.zst"
if smoke_bytes="$(s3_head_bytes "$smoke_key")"; then
  pass "smoke-test accession $SMOKE_ACC present at \$LOGAN_C/$SMOKE_ACC/... ($smoke_bytes bytes, $(human_size "$smoke_bytes"))"
else
  fail "smoke-test accession $SMOKE_ACC NOT found at $LOGAN_HTTPS/$smoke_key"
fi

# --- 4. Stats parquet: confirm it exists and report size. Do NOT download --
# by default — it is ~1.7 GB. Pass --download to actually fetch it.
echo
echo "--- stats parquet: \$LOGAN_STATS_PARQUET ---"
stats_key="stats/$LOGAN_STATS_PARQUET"
if stats_bytes="$(s3_head_bytes "$stats_key")"; then
  pass "$LOGAN_STATS_PARQUET present ($stats_bytes bytes, $(human_size "$stats_bytes")) — not downloaded (pass --download to fetch)"
  if [ "$DOWNLOAD_STATS" = 1 ]; then
    mkdir -p "$REF"
    info "fetching to \$REF/$LOGAN_STATS_PARQUET ..."
    if [ "$BACKEND" = aws ]; then
      aws s3 cp "$LOGAN_BUCKET/$stats_key" "$REF/$LOGAN_STATS_PARQUET" --no-sign-request
    else
      curl -fSL "$LOGAN_HTTPS/$stats_key" -o "$REF/$LOGAN_STATS_PARQUET"
    fi
    info "downloaded: $(du -h "$REF/$LOGAN_STATS_PARQUET" | cut -f1)"
  fi
else
  fail "$LOGAN_STATS_PARQUET NOT found at $LOGAN_HTTPS/$stats_key"
fi

# --- 5. $LOGAN_OBELISK reference artifacts (C-02) — all must exist. --------
echo
echo "--- \$LOGAN_OBELISK reference artifacts (C-02) ---"
obelisk_rel="${LOGAN_OBELISK#"$LOGAN_BUCKET"/}"
OBELISK_ARTIFACTS=(
  "01diamond_hmm/Obelisk_paper.hmm"
  "01diamond_hmm/Round1_Oblin1_sig.fasta"
  "01diamond_hmm/Round1_Oblin2_sig.fasta"
  "02_Building_Obelisk_DB/Obelisk_sig_nt.fasta"
  "03_Obelisk_DB_QC/Obelisk_db_centroids_nt.fasta"
  "03_Obelisk_DB_QC/Obelisk_cen_cen_orf.dmnd"
  "02_Building_Obelisk_DB/case1/novel_sequences.fasta"
  "03_Obelisk_DB_QC/Circularity_presence.tsv"
  "03_Obelisk_DB_QC/Domain_A_presence.tsv"
)
for rel in "${OBELISK_ARTIFACTS[@]}"; do
  key="$obelisk_rel/$rel"
  if bytes="$(s3_head_bytes "$key")"; then
    pass "$rel ($bytes bytes, $(human_size "$bytes"))"
  else
    fail "$rel NOT found at $LOGAN_HTTPS/$key"
  fi
done

# --- final summary -----------------------------------------------------------
echo
echo "=============================================================================="
echo " SUMMARY — ${#SUMMARY[@]} checks run, $FAILS failed"
echo "=============================================================================="
for line in "${SUMMARY[@]}"; do echo " $line"; done

if [ "$RELEASE_STALE" = 1 ]; then
  echo
  echo "RESULT: FAIL (exit 3) — the C-01 guard tripped: Logan's advertised release"
  echo "        has moved off v1.2. Do not proceed until config.sh is updated."
  exit 3
elif [ "$FAILS" -gt 0 ]; then
  echo
  echo "RESULT: FAIL (exit 1) — $FAILS check(s) failed. Do not proceed to Phase 1"
  echo "        until every check above reads PASS."
  exit 1
else
  echo
  echo "RESULT: PASS — Logan is reachable and config.sh's ground truth holds."
  echo "        Safe to proceed to Phase 1 (scripts/01_fetch_references.sh)."
  exit 0
fi
