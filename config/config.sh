#!/usr/bin/env bash
# =============================================================================
# obelisk-hunt — single source of truth for all paths and thresholds.
# Every script sources this. No script hardcodes a path or a threshold.
#
# Ground truth here was re-verified against primary sources on 2026-08-19.
# See docs/CORRECTIONS.md before changing any value in the "Logan" block —
# several of these differ deliberately from the original PLAN.md.
# =============================================================================

# --- Project root -----------------------------------------------------------
# Resolves to the repository root regardless of where the project is cloned,
# so the tree works both in-repo and when deployed to ~/obelisk-hunt.
# Override by exporting ROOT before sourcing.
if [ -z "${ROOT:-}" ]; then
  _cfg_src="${BASH_SOURCE[0]:-$0}"
  ROOT="$(cd "$(dirname "$(readlink -f "$_cfg_src")")/.." && pwd)"
  unset _cfg_src
fi
export ROOT
export REF="$ROOT/ref"
export WORK="$ROOT/work"
export LOGS="$ROOT/logs"
export RESULTS="$ROOT/results"
export SCRIPTS="$ROOT/scripts"
export VNOM="$ROOT/VNom"
export CIRCUCLUST="$ROOT/circuclust"

mkdir -p "$REF" "$WORK" "$LOGS" "$RESULTS/figures"

# --- Logan release ----------------------------------------------------------
# CORRECTION C-01: PLAN.md asserted "Logan v1 freeze = December 2023" and
# explicitly retracted a (correct) earlier claim of a Dec-2025 / 87-Pbp release.
# The Logan README states the current release is v1.2 (published 21 Apr 2026),
# built over a 31 Dec 2025 SRA freeze. Verified 2026-08-19.
#
# LOGAN_RELEASE selects which contig set the pipeline streams:
#   v1.2  -> s3://logan-pub/c/     38,124,741 accessions   (current, DEFAULT)
#   v1.0  -> s3://logan-pub/c1.0/  27,269,310 accessions   (replication control)
#
# v1.0 is retained ONLY so you can reproduce the published Obelisk sweep, which
# was run against it. Stats-v1.1.md (29 Jan 2025) says c1.0/ "will be deleted in
# ~1 year" — that window has passed, so treat c1.0/ as liable to vanish and
# check scripts/00_verify_logan.sh output before relying on it.
export LOGAN_RELEASE="${LOGAN_RELEASE:-v1.2}"

case "$LOGAN_RELEASE" in
  v1.2)
    export LOGAN_FREEZE="2025-12-31"
    export LOGAN_C="s3://logan-pub/c"
    export LOGAN_STATS_PARQUET="logan-seqstats-contigs-v1.2.parquet"
    # Exact column name. CORRECTION C-07: PLAN.md selected this column by
    # substring match, which is ambiguous on the v1 parquet (two columns match)
    # and silently yields nothing if the schema shifts.
    export LOGAN_SIZE_COL="contigs_after_compression_bytes"
    export LOGAN_N_ACCESSIONS=38124741
    ;;
  v1.0|v1)
    export LOGAN_FREEZE="2023-12-10"
    export LOGAN_C="s3://logan-pub/c1.0"
    export LOGAN_STATS_PARQUET="logan-seqstats.parquet"
    export LOGAN_SIZE_COL="size_contigs_after_compression"
    export LOGAN_N_ACCESSIONS=27269310
    ;;
  *)
    echo "config.sh: unknown LOGAN_RELEASE='$LOGAN_RELEASE' (expected v1.2 or v1.0)" >&2
    return 1 2>/dev/null || exit 1
    ;;
esac

# Unitigs are release-independent in path terms; v1.2 unitigs superseded v1's.
export LOGAN_U="s3://logan-pub/u"
export LOGAN_HTTPS="https://s3.amazonaws.com/logan-pub"
export LOGAN_BUCKET="s3://logan-pub"

# Published Obelisk analysis by the Logan/Obelisk authors, in the same bucket.
# CORRECTION C-02: PLAN.md did not know this existed and built the reference
# set from scratch off a journal supplementary table instead.
export LOGAN_OBELISK="$LOGAN_BUCKET/paper/Obelisk"

# --- Search thresholds ------------------------------------------------------
export DIAMOND_EVALUE="1e-5"
export DIAMOND_SENS="--very-sensitive"
export HMM_EVALUE="1e-10"
# Relaxed values used only by the documented fallback branches at each GATE.
export DIAMOND_EVALUE_RELAXED="1e-3"
export HMM_EVALUE_RELAXED="1e-5"

# --- Validation thresholds --------------------------------------------------
export MIN_INDEPENDENT_BIOPROJECTS=3
export OBELISK_MIN_LEN=800
export OBELISK_MAX_LEN=1200
# Terminal-repeat length for the circularity heuristic (scripts/09).
export CIRC_MIN_OVERLAP=20
# Identity at or above which a candidate counts as an already-known Obelisk.
export KNOWN_IDENTITY_PCT=95

# --- Compute ----------------------------------------------------------------
export THREADS="${THREADS:-16}"
# Per-accession threads during the triage sweep (many small concurrent jobs).
export TRIAGE_THREADS="${TRIAGE_THREADS:-2}"

# --- Reference file names (built once by scripts/01 and 02, never edited) ----
export REF_OBLIN_HMM="$REF/oblin1.hmm"
export REF_OBLIN_FAA="$REF/oblin1_centroids.faa"
export REF_OBLIN_DMND="$REF/oblin1_dmnd"
export REF_OBELISK_NT="$REF/obelisk_nt.fna"
export REF_DECOY_FAA="$REF/decoy_shuffled.faa"
export REF_DECOY_DMND="$REF/decoy_dmnd"
