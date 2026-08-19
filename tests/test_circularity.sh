#!/usr/bin/env bash
# Regression test for scripts/09_circularity_check.py — corrections C-12 and the
# low-complexity terminal-repeat rejection. Builds two synthetic candidates:
#
#   cand_origin_spanning   circular, 30 bp head-to-tail repeat, ORF crossing the origin
#                          -> must be called circular AND the ORF must be found
#   cand_homopolymer_ends  25 bp poly-A at both ends, no real repeat
#                          -> must NOT be called circular
#
# The second case is the one the plan's bare first-20bp == last-20bp test gets wrong.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
source config/config.sh
RESULTS="$(mktemp -d)/results"; export RESULTS; mkdir -p "$RESULTS"

python3 - <<'PY'
import os, random
random.seed(7)
R = os.environ["RESULTS"]
cods = ["GCT","CGT","AAT","GAT","TGT","CAA","GAA","GGT","CAT","ATT",
        "TTA","AAA","ATG","TTT","CCT","TCT","ACT","TGG","TAT","GTT"]
orf    = "ATG" + "".join(random.choice(cods) for _ in range(158)) + "TAA"
filler = "".join(random.choice("ACGT") for _ in range(520))
linear = orf + filler
rot    = linear[300:] + linear[:300]      # rotate so the ORF straddles the origin
circ   = rot + rot[:30]                   # head-to-tail terminal repeat
neg    = "A"*25 + "".join(random.choice("ACGT") for _ in range(950)) + "A"*25
with open(f"{R}/candidates.fna", "w") as fh:
    fh.write(f">cand_origin_spanning\n{circ}\n>cand_homopolymer_ends\n{neg}\n")
PY

python3 scripts/09_circularity_check.py >/dev/null
T="$RESULTS/circularity.tsv"

fail=0
chk() { # chk <id> <column> <expected>
  local got; got=$(awk -F'\t' -v id="$1" -v c="$2" \
      'NR==1{for(i=1;i<=NF;i++)h[$i]=i} $1==id{print $h[c]}' "$T")
  if [ "$got" = "$3" ]; then echo "  ok   $1.$2 = $3"
  else echo "  FAIL $1.$2: expected '$3', got '$got'"; fail=1; fi
}
echo "test_circularity:"
chk cand_origin_spanning  circular            True
chk cand_origin_spanning  orf_spans_origin    True
chk cand_homopolymer_ends circular            False
chk cand_homopolymer_ends term_repeat_reject  homopolymer
[ "$fail" -eq 0 ] && echo "PASS" || { echo "FAILED"; exit 1; }
