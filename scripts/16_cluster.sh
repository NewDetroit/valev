#!/usr/bin/env bash
# =============================================================================
# Task 14 Step 5 — circuclust nomenclature.
#
# The field's convention is Obelisk_X_Y_Z: X = cluster at 80% nt identity,
# Y = at 95%, Z = strain. Following it makes results directly comparable to
# published work. Real examples to validate against are in the bucket:
#   aws s3 ls s3://logan-pub/paper/Obelisk/03_Obelisk_DB_QC/ --no-sign-request
# e.g. Obelisk_000001_000001_000001.
#
# ---------------------------------------------------------------------------
# CORRECTION C-14 — every circuclust flag in PLAN.md is wrong. It writes:
#
#     "$CIRCUCLUST/bin/circuclust" --input "$WORK/all_obelisks.fna" \
#         --id 0.80 --output "$RESULTS/clusters80.tsv"
#
# Checked against circuclust's README (fetched 2026-08-19), whose clustering
# example is:
#
#     circuclust -cluster seqs.fa -id 0.9 -fastaout centroids.fa -tsvout hits.tsv
#
# Single dashes, not double. The input is the VALUE of -cluster, not --input.
# Output splits into -fastaout (centroids) and -tsvout (assignments); there is
# no --output. The plan's invocation would fail outright.
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

CC="${CIRCUCLUST_BIN:-$CIRCUCLUST/bin/circuclust}"
[ -x "$CC" ] || CC="$CIRCUCLUST/circuclust"
[ -x "$CC" ] || { echo "circuclust not found or not executable (tried $CIRCUCLUST/bin/circuclust and $CIRCUCLUST/circuclust)" >&2
                  echo "  Download a release binary: https://github.com/rcedgar/circuclust/releases" >&2
                  echo "  Or set CIRCUCLUST_BIN=/path/to/circuclust" >&2; exit 2; }
[ -s "$RESULTS/validated.fna" ] || { echo "missing $RESULTS/validated.fna" >&2; exit 2; }
[ -s "$REF_OBELISK_NT" ]        || { echo "missing $REF_OBELISK_NT" >&2; exit 2; }

ALL="$WORK/all_obelisks.fna"
cat "$RESULTS/validated.fna" "$REF_OBELISK_NT" > "$ALL"
echo "clustering $(grep -c '^>' "$ALL") sequences (new + published)"

# Both levels of the naming convention.
"$CC" -cluster "$ALL" -id 0.80 -fastaout "$RESULTS/centroids80.fa" -tsvout "$RESULTS/clusters80.tsv"
"$CC" -cluster "$ALL" -id 0.95 -fastaout "$RESULTS/centroids95.fa" -tsvout "$RESULTS/clusters95.tsv"

echo "wrote $RESULTS/clusters80.tsv  ($(wc -l < "$RESULTS/clusters80.tsv") rows)"
echo "wrote $RESULTS/clusters95.tsv  ($(wc -l < "$RESULTS/clusters95.tsv") rows)"
echo
echo "Assign Obelisk_X_Y_Z names by joining the two tables: X from the 80% cluster,"
echo "Y from the 95% cluster within it, Z per distinct strain. Cross-check the format"
echo "against the published names in s3://logan-pub/paper/Obelisk/03_Obelisk_DB_QC/"
echo "before adopting your numbering — colliding with an existing name is worse than"
echo "not naming at all."
