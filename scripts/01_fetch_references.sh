#!/usr/bin/env bash
# =============================================================================
# PHASE 1 / Task 3 — build the reference set that everything downstream
# searches against.  A wrong reference set does not crash anything; it silently
# invalidates every result in the project.  Hence the assertions.
#
# CORRECTION C-02.  PLAN.md Task 3 Step 1 called a manual download of the
# Zheludev et al. *Cell* Table S1 "the backbone of the entire project".  It is
# not the best available source and it is not reproducible from a script.  The
# Obelisk authors published their own Logan-wide Oblin sweep at
# s3://logan-pub/paper/Obelisk/ — anonymously readable, already QC'd, already
# clustered, and already named in the field's Obelisk_X_Y_Z convention.  That is
# the PRIMARY path below.  The plan's from-scratch build survives behind
# --from-scratch for anyone who wants to reproduce it, or who needs a reference
# set independent of the incumbents'.
#
# Also see C-05 (idempotency guards test existence, not non-emptiness) and the
# exit-code convention in docs/INTERFACES:
#   0 ok · 1 unexpected error · 2 required input missing · 3 assertion failed
# =============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config/config.sh"

FROM_SCRATCH=0
FORCE=0
STRICT_SIZES=0
SKIP_TOOL_DOCS=0
SEQ_COL=""
ID_COL=""

usage() {
  cat <<'USAGE'
Usage: 01_fetch_references.sh [options]

  (default)          fetch the published paper/Obelisk artifacts  [PRIMARY, C-02]
  --from-scratch     rebuild from ref/obelisk_stringent.tsv instead (plan Task 3)
  --seq-col NAME     --from-scratch: sequence column in that table
  --id-col NAME      --from-scratch: identifier column in that table
  --force            re-download even if the local file already exists
  --strict-sizes     treat an upstream size change as fatal, not as a warning
  --skip-tool-docs   do not fetch the VNom / circuclust READMEs
  -h, --help         this text
USAGE
}

while [ $# -gt 0 ]; do
  case "$1" in
    --from-scratch)   FROM_SCRATCH=1 ;;
    --force)          FORCE=1 ;;
    --strict-sizes)   STRICT_SIZES=1 ;;
    --skip-tool-docs) SKIP_TOOL_DOCS=1 ;;
    --seq-col)        SEQ_COL="${2:?--seq-col needs a value}"; shift ;;
    --id-col)         ID_COL="${2:?--id-col needs a value}"; shift ;;
    -h|--help)        usage; exit 0 ;;
    *) echo "01_fetch_references.sh: unknown option '$1'" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

PROV="$REF/PROVENANCE.tsv"
NOW() { date -u +%Y-%m-%dT%H:%M:%SZ; }
have() { command -v "$1" >/dev/null 2>&1; }
warn() { echo "WARNING: $*" >&2; }

sha256_of() {
  if have sha256sum; then sha256sum "$1" | awk '{print $1}'
  elif have shasum;  then shasum -a 256 "$1" | awk '{print $1}'
  else echo "NA"; fi
}

# s3://logan-pub/... -> https://s3.amazonaws.com/logan-pub/...  Both endpoints
# were confirmed to serve these objects anonymously on 2026-08-19; the HTTPS one
# is the fallback for machines without the aws CLI (this build container, for one).
https_of() { printf '%s%s' "$LOGAN_HTTPS" "${1#"$LOGAN_BUCKET"}"; }

remote_size() {
  curl -sfI --retry 3 --retry-delay 2 "$1" 2>/dev/null \
    | tr -d '\r' \
    | awk 'tolower($1)=="content-length:"{print $2}' \
    | tail -1
}

# ref/PROVENANCE.tsv is mandatory (docs/INTERFACES): every reference file must be
# traceable to the URI it came from, because that is what lets the write-up cite
# its inputs.  Rows are upserted by local_path so re-running does not duplicate.
prov_init() {
  [ -e "$PROV" ] || printf 'local_path\tsource_uri\tbytes\tsha256\tfetched_utc\n' > "$PROV"
}
prov_record() {  # <path-relative-to-ROOT> <source-uri>
  local rel="$1" uri="$2" abs="$ROOT/$1" bytes sha
  [ -e "$abs" ] || { echo "prov_record: missing $abs" >&2; exit 1; }
  bytes="$(wc -c < "$abs" | tr -d ' ')"
  sha="$(sha256_of "$abs")"
  prov_init
  awk -F'\t' -v p="$rel" 'NR==1 || $1!=p' "$PROV" > "$PROV.part"
  printf '%s\t%s\t%s\t%s\t%s\n' "$rel" "$uri" "$bytes" "$sha" "$(NOW)" >> "$PROV.part"
  mv "$PROV.part" "$PROV"
}

# The one confusion this pipeline cannot tolerate is a nucleotide file in a
# protein slot or the reverse — it is the same class of error as C-04, and it
# produces plausible-looking garbage rather than a crash.  Every FASTA is
# checked before it is allowed to become a reference file.
assert_alphabet() {  # <file> <nt|aa> <label>
  python3 - "$1" "$2" "$3" <<'PY'
import collections, sys
path, want, label = sys.argv[1], sys.argv[2], sys.argv[3]
NT = set("ACGTUNRYSWKMBDHV-*.")
AA = set("ACDEFGHIKLMNPQRSTVWYBZXUO-*.")
# Every nucleotide letter is also a legal amino-acid letter, so "is this
# protein?" can only be answered by the letters that are amino-acid-only.
AA_ONLY = set("EFILPQZJ")
res = collections.Counter()
n = 0
with open(path, "r", errors="replace") as fh:
    for line in fh:
        if line.startswith(">"):
            n += 1
            if n > 500:
                break
            continue
        res.update(line.strip().upper())
tot = sum(res.values())
if tot == 0:
    sys.exit(f"{label}: {path} contains no residues")
f_nt = sum(v for k, v in res.items() if k in NT) / tot
f_aa = sum(v for k, v in res.items() if k in AA) / tot
f_only = sum(v for k, v in res.items() if k in AA_ONLY) / tot
print(f"  {label}: {min(n,500)} records sampled, "
      f"nt={f_nt:.4f} aa={f_aa:.4f} aa-only={f_only:.4f}")
if want == "nt" and (f_nt < 0.99 or f_only > 0.01):
    sys.exit(f"{label}: {path} is not nucleotide (nt={f_nt:.4f}, aa-only={f_only:.4f})")
if want == "aa" and (f_aa < 0.99 or f_only < 0.01):
    sys.exit(f"{label}: {path} is not protein (aa={f_aa:.4f}, aa-only={f_only:.4f})")
PY
}

# Downloads land on a .part file and are renamed only once complete.  C-05 makes
# every guard in this project test existence rather than non-emptiness, which
# means a half-written file from an interrupted run would be skipped forever.
fetch_object() {  # <s3-key-under-paper/Obelisk> <path-relative-to-ROOT> <expected-bytes>
  local key="$1" rel="$2" expect="$3"
  local uri="$LOGAN_OBELISK/$key" dest="$ROOT/$rel" url got remote
  url="$(https_of "$uri")"

  if [ -e "$dest" ] && [ "$FORCE" -eq 0 ]; then
    echo "  skip (present): $rel"
  else
    mkdir -p "$(dirname "$dest")"
    echo "  fetch: $uri"
    if have aws; then
      if ! aws s3 cp "$uri" "$dest.part" --no-sign-request --only-show-errors; then
        warn "aws s3 cp failed for $key; retrying over HTTPS"
        curl -fSL --retry 5 --retry-delay 3 -o "$dest.part" "$url"
      fi
    else
      curl -fSL --retry 5 --retry-delay 3 -o "$dest.part" "$url"
    fi
    mv "$dest.part" "$dest"
  fi

  got="$(wc -c < "$dest" | tr -d ' ')"
  remote="$(remote_size "$url" || true)"
  if [ -n "$remote" ] && [ "$got" != "$remote" ]; then
    echo "ERROR: $rel is $got bytes but the object is $remote — truncated download" >&2
    exit 1
  fi
  # A silent upstream replacement of a reference file is exactly the failure
  # that invalidates a project without breaking it, so say so loudly.  Not fatal
  # by default: the authors are entitled to update their own artifacts, and
  # PROVENANCE.tsv records the sha256 that was actually used.
  if [ "$got" != "$expect" ]; then
    if [ "$STRICT_SIZES" -eq 1 ]; then
      echo "ERROR: $rel is $got bytes, expected $expect (verified 2026-08-19)" >&2
      exit 1
    fi
    warn "$rel is $got bytes, expected $expect (verified 2026-08-19) — upstream changed"
  fi
  prov_record "$rel" "$uri"
}

# -----------------------------------------------------------------------------
# The manifest.  key | destination | alphabet | bytes-verified-2026-08-19
#
# Role assignments, and why each object won its slot:
#
#   Round1_Oblin1_sig.fasta -> $REF_OBLIN_FAA
#       The authors' Round-1 significant Oblin-1 ORFs: 38,742 protein records,
#       mean 131 aa.  This is the DIAMOND triage query set.  Taken verbatim and
#       unfiltered — the short partial ORFs in it cannot reach E<=1e-5 on their
#       own, and dropping records would be an unlogged change to a published
#       reference set.  (The config name says "centroids"; the file is the
#       significant-hit set, not a clustered centroid set.  Name kept because
#       config/config.sh and docs/INTERFACES are binding.)
#
#   Round1_Oblin2_sig.fasta -> ref/oblin2.faa
#       194 Oblin-2 protein records.  PLAN.md ignores Oblin-2 entirely; an
#       Obelisk with a divergent Oblin-1 and a recognisable Oblin-2 is a real
#       and reachable case, so the set is fetched even though the plan's
#       pipeline does not yet query it.
#
#   Obelisk_db_centroids_nt.fasta -> $REF_OBELISK_NT
#       10,857 nucleotide centroids under the canonical Obelisk_X_Y_Z names.
#       Chosen over Obelisk_sig_nt.fasta (25.7 MB, 'known-Obelisk' screening at
#       $KNOWN_IDENTITY_PCT wants deduplicated representatives, and a hit here
#       reports as "Obelisk_000123" rather than as a raw Logan contig id) and
#       over Obelisk_cen_cen.fasta (a further-collapsed 1.8 MB set that would
#       under-report re-discoveries).
#
#   Obelisk_sig_nt.fasta, case1/novel_sequences.fasta, Obelisk_cen_cen_orf.fa
#       Fetched as cross-checks, not as query sets.  novel_sequences.fasta is
#       the authors' own novel calls — anything this project rediscovers that is
#       in there is theirs, not ours, which is the core of C-02's warning.
#       Obelisk_cen_cen_orf.fa is every ORF of every centroid including 15-aa
#       fragments; it is the input to their published .dmnd and is useless as a
#       query set, so it is deliberately NOT $REF_OBLIN_FAA.
#
#   Obelisk_paper.hmm
#       Fetched here, turned into $REF_OBLIN_HMM by scripts/02 (which owns that
#       path per docs/INTERFACES).  Contains TWO profiles: COR_ORF_1 (Oblin-1,
#       LENG 198, NSEQ 575) and COR_ORF_2 (Oblin-2, LENG 180, NSEQ 134).
#
#   obelisk_name_mapping.tsv
#       1,743 rows mapping Logan contigs to Obelisk_X_Y_Z clusters — the
#       published stringent set.  It is also the evidence for the positive
#       control's identity; see PC_ID below.
# -----------------------------------------------------------------------------
MANIFEST=(
  "01diamond_hmm/Obelisk_paper.hmm|ref/obelisk_paper.hmm|none|177803"
  "01diamond_hmm/Round1_Oblin1_sig.fasta|ref/oblin1_centroids.faa|aa|9704316"
  "01diamond_hmm/Round1_Oblin2_sig.fasta|ref/oblin2.faa|aa|37160"
  "01diamond_hmm/obelisk_name_mapping.tsv|ref/obelisk_name_mapping.tsv|none|104471"
  "02_Building_Obelisk_DB/Obelisk_sig_nt.fasta|ref/obelisk_sig_nt.fna|nt|25653784"
  "03_Obelisk_DB_QC/Obelisk_db_centroids_nt.fasta|ref/obelisk_nt.fna|nt|7698335"
  "03_Obelisk_DB_QC/Obelisk_cen_cen_orf.fa|ref/obelisk_cen_cen_orf.faa|aa|8484393"
  "02_Building_Obelisk_DB/case1/novel_sequences.fasta|ref/obelisk_authors_novel_nt.fna|nt|14340693"
  "03_Obelisk_DB_QC/Circularity_presence.tsv|ref/Circularity_presence.tsv|none|377268"
  "03_Obelisk_DB_QC/Domain_A_presence.tsv|ref/Domain_A_presence.tsv|none|370921"
)

# Obelisk-S.s, the GATE 1 positive control (scripts/03_positive_control.sh).
# Identity established from two independent artifacts in the bucket:
#   1. obelisk_name_mapping.tsv maps the contig 'var.obelisk-1.Ss_' — the only
#      strain-named, non-SRA entry in all 1,743 rows — to this cluster.
#   2. the record is exactly 1,137 nt, the published length of Obelisk-S.s.
# Both were checked live on 2026-08-19.  PC_LEN is asserted below so a change to
# either the id scheme or the sequence fails here rather than at GATE 1.
PC_ID="Obelisk_000003_000001_000001"
PC_LEN=1137

echo "=== 01_fetch_references.sh — $(NOW) ==="
echo "ROOT=$ROOT"
echo "LOGAN_OBELISK=$LOGAN_OBELISK"
have aws || warn "aws CLI not found — falling back to anonymous HTTPS reads"
prov_init

if [ "$FROM_SCRATCH" -eq 0 ]; then
  # ---------------------------------------------------------------------------
  # PRIMARY — the published artifacts (C-02)
  # ---------------------------------------------------------------------------
  echo "--- primary path: published paper/Obelisk artifacts ---"
  for entry in "${MANIFEST[@]}"; do
    IFS='|' read -r key rel kind expect <<<"$entry"
    fetch_object "$key" "$rel" "$expect"
    case "$kind" in
      nt|aa) assert_alphabet "$ROOT/$rel" "$kind" "$(basename "$rel")" ;;
    esac
  done
else
  # ---------------------------------------------------------------------------
  # FALLBACK — PLAN.md Task 3, rebuilt from the Cell supplementary table.
  # Kept because it is the only path that does not depend on the incumbents'
  # own artifacts.  Step 1 is unavoidably manual: the table is behind the
  # publisher's site and cannot be fetched non-interactively.
  # ---------------------------------------------------------------------------
  echo "--- fallback path: --from-scratch from the Cell supplementary table ---"
  SUPP="$REF/obelisk_stringent.tsv"
  if [ ! -e "$SUPP" ]; then
    cat >&2 <<EOF
ERROR: $SUPP not found.

  --from-scratch rebuilds the reference set from Zheludev et al. Table S1
  (doi:10.1016/j.cell.2024.09.033; preprint doi:10.1101/2024.01.20.576352),
  file 'Supplementary_table_1_stringent_Obelisk_clustering_011724.tsv'.
  Download it by hand and save it as:
      $SUPP
  Then re-run with --from-scratch --seq-col <NAME> --id-col <NAME>.

  Or drop --from-scratch and use the published artifacts, which need no manual
  step and are what the rest of this pipeline is calibrated against.
EOF
    exit 2
  fi

  if [ -z "$SEQ_COL" ] || [ -z "$ID_COL" ]; then
    echo "ERROR: --from-scratch needs --seq-col and --id-col." >&2
    echo "Columns present in $SUPP:" >&2
    head -1 "$SUPP" | tr '\t' '\n' | cat -n >&2
    exit 2
  fi

  python3 "$SCRIPTS/tsv_to_fasta.py" "$SUPP" "$SEQ_COL" "$ID_COL" \
    "$REF_OBELISK_NT" --alphabet nt --rna-to-dna
  prov_record "ref/obelisk_nt.fna" "file://$SUPP (Zheludev et al. Cell Table S1, manual download)"
  prov_record "ref/obelisk_stringent.tsv" "manual: Zheludev et al. Cell Table S1"

  # Plan Task 3 Step 5: proteins either come out of the table directly or have
  # to be translated from ORF coordinates.  Both need input only the operator
  # has, so fail with the instruction rather than guessing an ORF caller.
  if [ ! -e "$REF_OBLIN_FAA" ]; then
    cat >&2 <<EOF
ERROR: $REF_OBLIN_FAA not built.

  --from-scratch produced the nucleotide set but cannot produce the Oblin-1
  protein set without one of:
    (a) a protein column in $SUPP — then run, by hand:
          python3 $SCRIPTS/tsv_to_fasta.py $SUPP <PROT_COL> $ID_COL \\
            $REF_OBLIN_FAA --alphabet aa
    (b) ORF coordinates as BED — then run, by hand (plan Task 3 Step 5):
          seqkit subseq --bed <orf_coords.bed> $REF_OBELISK_NT \\
            | seqkit translate --frame 1 --trim > $REF_OBLIN_FAA
EOF
    exit 2
  fi
  assert_alphabet "$REF_OBELISK_NT" nt "obelisk_nt.fna"
  assert_alphabet "$REF_OBLIN_FAA" aa "oblin1_centroids.faa"
fi

# -----------------------------------------------------------------------------
# Positive control — ref/positive_control.fna (docs/INTERFACES, producer T3)
# -----------------------------------------------------------------------------
echo "--- positive control: $PC_ID ---"
if [ -e "$REF/positive_control.fna" ] && [ "$FORCE" -eq 0 ]; then
  echo "  skip (present): ref/positive_control.fna"
else
  python3 - "$REF_OBELISK_NT" "$REF/positive_control.fna" "$PC_ID" "$PC_LEN" <<'PY'
import sys
src, dst, want, want_len = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
seq, hit = [], False
with open(src, "r", errors="replace") as fh:
    for line in fh:
        if line.startswith(">"):
            if hit:
                break
            hit = line[1:].split()[0] == want
            continue
        if hit:
            seq.append(line.strip())
s = "".join(seq).upper()
if not s:
    sys.exit(f"positive control: {want} not found in {src}")
if len(s) != want_len:
    sys.exit(f"positive control: {want} is {len(s)} nt, expected {want_len} — "
             f"the reference set changed; re-verify the Obelisk-S.s identity "
             f"before trusting GATE 1")
bad = set(s) - set("ACGT")
if bad:
    sys.exit(f"positive control: {want} contains non-DNA characters {sorted(bad)}")
with open(dst + ".part", "w") as out:
    out.write(f">{want}\n")
    for i in range(0, len(s), 60):
        out.write(s[i:i+60] + "\n")
import os
os.replace(dst + ".part", dst)
print(f"  wrote {dst}: {want}, {len(s)} nt")
PY
fi
prov_record "ref/positive_control.fna" "$LOGAN_OBELISK/03_Obelisk_DB_QC/Obelisk_db_centroids_nt.fasta#$PC_ID"

# -----------------------------------------------------------------------------
# Tool documentation.  CORRECTIONS OPEN-C: api.github.com is blocked from the
# build environment but raw.githubusercontent.com is not, so the two READMEs
# that hold the CLI details the plan leaves as placeholders (VNom's entry point,
# circuclust's flags) are fetched here rather than left to be looked up later.
# Cloning the repos themselves is plan Task 2 (scripts/environment setup).
# -----------------------------------------------------------------------------
if [ "$SKIP_TOOL_DOCS" -eq 0 ]; then
  echo "--- tool READMEs ---"
  mkdir -p "$REF/tools"
  # VNom is on 'main', circuclust on 'master'; try both rather than hardcode.
  for spec in "Zheludev/VNom|VNom_README.md" "rcedgar/circuclust|circuclust_README.md"; do
    IFS='|' read -r repo out <<<"$spec"
    dest="$REF/tools/$out"
    if [ -e "$dest" ] && [ "$FORCE" -eq 0 ]; then
      echo "  skip (present): ref/tools/$out"
      continue
    fi
    ok=0
    for branch in main master; do
      url="https://raw.githubusercontent.com/$repo/$branch/README.md"
      if curl -fsSL --retry 3 --retry-delay 2 -o "$dest.part" "$url"; then
        mv "$dest.part" "$dest"
        prov_record "ref/tools/$out" "$url"
        echo "  fetched: ref/tools/$out ($branch)"
        ok=1
        break
      fi
    done
    # Non-fatal: these are documentation, not pipeline inputs.  Everything that
    # IS a pipeline input fails hard above.
    [ "$ok" -eq 1 ] || { rm -f "$dest.part"; warn "could not fetch README for $repo"; }
  done
fi

# -----------------------------------------------------------------------------
echo "--- reference set ---"
for f in "$REF_OBELISK_NT" "$REF_OBLIN_FAA" "$REF/oblin2.faa" "$REF/positive_control.fna"; do
  if [ ! -e "$f" ]; then
    echo "ERROR: contract file missing after run: $f" >&2
    exit 1
  fi
  printf '  %-34s %10s bytes  %6s records\n' \
    "$(basename "$f")" "$(wc -c < "$f" | tr -d ' ')" "$(grep -c '^>' "$f" || true)"
done
echo "  PROVENANCE.tsv rows: $(( $(wc -l < "$PROV") - 1 ))"
echo "References ready.  Next: scripts/02_build_search_dbs.sh"
