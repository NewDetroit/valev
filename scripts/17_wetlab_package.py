#!/usr/bin/env python3
"""Task 15 — rank candidates for wet-lab confirmation and write the request brief.

Reads  : $RESULTS/replication_matrix.tsv, $RESULTS/validated.fna
         $RESULTS/junction_support.tsv, $RESULTS/circularity.tsv   (optional, used if present)
Writes : $RESULTS/wetlab_priority.tsv
         $RESULTS/primer3_input.txt      -> feed to Primer3 for divergent primers
         $RESULTS/WETLAB_BRIEF.md

Ranking is by INDEPENDENCE, not abundance: n_bioprojects*3 + n_centers*2 + n_platforms.
The most-replicated element is the least likely to be an artifact, and abundance within a
single project is exactly what a contaminant looks like.

ON PRIMERS. The plan asks for outward-facing (divergent) primers, which amplify only across
the circular junction and give no product from a linear template — the standard RT-PCR test
for circularity. Primer3 is not assumed to be installed, and hand-rolled primers with no Tm
or secondary-structure checking are worse than none: they waste bench time and can produce a
confident wrong answer. This emits correctly formatted Primer3 input with the divergent
design already set up, rather than guessing at sequences.
"""
from __future__ import annotations
import os, sys
from datetime import datetime, timezone

PRODUCT_MIN, PRODUCT_MAX = 100, 300


def env(name: str) -> str:
    v = os.environ.get(name)
    if not v:
        sys.stderr.write(f"[17_wetlab_package] {name} unset — run: source config/config.sh\n")
        raise SystemExit(2)
    return v


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]


def read_fasta(path):
    out, rid, buf = {}, None, []
    if not os.path.exists(path):
        return out
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if rid:
                out[rid] = "".join(buf)
            rid, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.strip())
    if rid:
        out[rid] = "".join(buf)
    return out


def main() -> int:
    RESULTS = env("RESULTS")
    min_bp = int(os.environ.get("MIN_INDEPENDENT_BIOPROJECTS", "3"))
    freeze = os.environ.get("LOGAN_FREEZE", "unknown")
    release = os.environ.get("LOGAN_RELEASE", "unknown")

    matrix = read_tsv(os.path.join(RESULTS, "replication_matrix.tsv"))
    if not matrix:
        sys.stderr.write("[17_wetlab_package] missing replication_matrix.tsv — run scripts/11 (exit 2)\n")
        return 2
    seqs = read_fasta(os.path.join(RESULTS, "validated.fna"))
    junction = {r["candidate"]: r for r in read_tsv(os.path.join(RESULTS, "junction_support.tsv"))}
    circ = {r["id"]: r for r in read_tsv(os.path.join(RESULTS, "circularity.tsv"))}

    def num(r, k):
        try:
            return int(r.get(k, 0) or 0)
        except ValueError:
            return 0

    for r in matrix:
        r["score"] = num(r, "n_bioprojects") * 3 + num(r, "n_centers") * 2 + num(r, "n_platforms")
    ranked = sorted(matrix, key=lambda r: -r["score"])
    top = ranked[:5]

    cols = ["candidate", "score", "n_accessions", "n_bioprojects", "n_centers",
            "n_platforms", "bioprojects"]
    out = os.path.join(RESULTS, "wetlab_priority.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in top:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"wrote {out}")
    for r in top:
        print(f"  {r['candidate'][:44]:<44} score {r['score']:>3}  "
              f"BP {r.get('n_bioprojects')}  ctr {r.get('n_centers')}  plat {r.get('n_platforms')}")

    # Primer3 input. Divergent primers are obtained by rotating the sequence so the
    # junction sits in the middle: primers flanking that midpoint point outward on the
    # original molecule and only close across a circular template.
    p3 = os.path.join(RESULTS, "primer3_input.txt")
    n_p3 = 0
    with open(p3, "w") as fh:
        for r in top[:3]:
            s = seqs.get(r["candidate"])
            if not s:
                continue
            half = len(s) // 2
            rotated = s[half:] + s[:half]          # junction now at position len(s)-half
            fh.write(
                f"SEQUENCE_ID={r['candidate']}_junction\n"
                f"SEQUENCE_TEMPLATE={rotated}\n"
                f"SEQUENCE_TARGET={len(s)-half-10},20\n"
                f"PRIMER_PRODUCT_SIZE_RANGE={PRODUCT_MIN}-{PRODUCT_MAX}\n"
                f"PRIMER_NUM_RETURN=3\nPRIMER_OPT_TM=60.0\nPRIMER_MIN_TM=57.0\n"
                f"PRIMER_MAX_TM=63.0\nPRIMER_EXPLAIN_FLAG=1\n=\n")
            n_p3 += 1
    print(f"wrote {p3} ({n_p3} templates, junction centred)")

    def tier(r):
        c = r["candidate"]
        j = junction.get(c, {}).get("verdict", "")
        bp = num(r, "n_bioprojects")
        if j == "CIRCULAR_SUPPORTED" and bp >= min_bp:
            return "A — replicated and circular by read mapping"
        if bp >= min_bp:
            return "B — replicated; circularity not yet confirmed on reads"
        return "C — tentative; insufficient independent replication"

    brief = os.path.join(RESULTS, "WETLAB_BRIEF.md")
    with open(brief, "w") as fh:
        fh.write(f"""# Wet-lab confirmation request — Obelisk-like elements

Generated {datetime.now(timezone.utc):%Y-%m-%d} from Logan {release} (SRA freeze {freeze}).

## What an Obelisk is

Obelisks are ~1 kb circular RNA elements, first described in 2024, that encode a protein
family called Oblin-1 and fold into a characteristic rod-like secondary structure. They
carry no capsid and have no known DNA intermediate, so they exist only in RNA sequencing
data. They appear to be widespread in host-associated microbiomes and almost nothing is
known about what they do.

## What we found

A computational search of this niche recovered {len(matrix)} candidate elements, of which
{sum(1 for r in matrix if num(r,'n_bioprojects') >= min_bp)} are supported by independent
replication across at least {min_bp} unrelated BioProjects.

The pipeline was validated first on a positive control: it recovers the known Obelisk-S.s
from *Streptococcus sanguinis* SK36 RNA-seq. Evidence in `results/positive_control_evidence.tsv`.

### Priority candidates

| # | Candidate | Evidence tier | BioProjects | Centers | Platforms | Length |
|---|---|---|---|---|---|---|
""")
        for n, r in enumerate(top, 1):
            c = r["candidate"]
            L = circ.get(c, {}).get("length", len(seqs.get(c, "")) or "?")
            fh.write(f"| {n} | `{c}` | {tier(r)} | {r.get('n_bioprojects')} | "
                     f"{r.get('n_centers')} | {r.get('n_platforms')} | {L} |\n")
        fh.write(f"""
Candidates are ranked by **independence**, not abundance. An element seen many times inside
one project is what a contaminant looks like; an element seen across unrelated projects,
sequencing centers and platforms is not.

## What we are asking for

**RT-PCR with outward-facing (divergent) primers**, on RNA from the sample types listed
below.

Divergent primers point away from each other on the linear sequence. On a linear template
they produce no product. On a circular template they amplify across the junction and give a
band of {PRODUCT_MIN}-{PRODUCT_MAX} bp. That single experiment distinguishes a genuinely
circular RNA from an assembly artifact, which is the main alternative explanation for
everything above.

Primer3 input with the junction already centred is in `results/primer3_input.txt`.
Run it through Primer3 and record the chosen pairs in `results/primers.tsv` before ordering.

## Sample sources

Source BioProjects per candidate are in `results/wetlab_priority.tsv`. Look each up in the
SRA to identify the sample type and origin, and note which are plausibly obtainable locally
— a rumen sample from a veterinary school, a poultry sample from an agricultural program.

## What a positive result would establish

A junction-spanning RT-PCR product, Sanger-confirmed, would demonstrate that the element is
a real circular RNA present in biological material rather than an artifact of assembly or
of the sequencing archive. Combined with the Oblin-1 homology and the rod-like fold, that is
sufficient to describe it as a new Obelisk-like element.

A negative result is also informative: it would indicate the element is either absent from
the available sample or not circular, and would stop us reporting it as a discovery.
""")
    print(f"wrote {brief}")
    print("\nBring the brief, not an idea. A one-page request with validated targets and")
    print("designed primers is a far easier yes than a conversation about a possibility.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
