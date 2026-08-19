#!/usr/bin/env python3
"""Task 11 — confirm head-to-tail circularity and find Oblin-like ORFs.

Reads  : $RESULTS/candidates.fna          (pooled VNom candidates, from scripts/08)
Writes : $RESULTS/circularity.tsv         (per-candidate evidence table)
         $WORK/candidate_orfs.faa         (longest ORF per candidate, >= --min-orf-aa)

Two deliberate departures from PLAN.md Task 11 Step 1:

C-12 — the plan translates only the LINEAR sequence. Obelisk genomes are circular and
Oblin-1 ORFs frequently cross the point at which the assembler linearised the molecule, so
a linear-only six-frame search truncates exactly the ORFs this project most wants to find.
Here the search runs on the DOUBLED sequence (S+S), ORF length is capped at one full turn
(an ORF cannot be longer than the molecule), ORFs starting at or beyond position L are
discarded as duplicates of the copy-1 hit, and origin-spanning ORFs are reported in the
separate `orf_spans_origin` column that the interface contract requires.

Terminal-repeat heuristic — the plan compares the first 20 bp to the last 20 bp and returns
a bare boolean. Two things are wrong with that as evidence. (1) A boolean throws away the
size of the observation: a 3 bp and a 300 bp terminal repeat are not the same claim, so this
script reports the length of the MAXIMAL exact terminal repeat. (2) On ~1 kb sequences a
20 bp homopolymer or dinucleotide run at both ends arises readily and is not evidence of
anything, so low-complexity repeats are rejected with the reason recorded.

None of this proves circularity. A terminal repeat is consistent with a circular molecule
that the assembler opened, and also with a tandem duplication, a chimeric join, or a
repetitive element. Physical proof is reads that align ACROSS the junction on a doubled
reference — scripts/12_junction_mapping.sh. Treat this table as supporting evidence only.
"""

import argparse
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq

# Rejection thresholds for the terminal-repeat complexity filter. These guard the heuristic
# against the failure mode it is most prone to on short sequences: a run of one or two bases
# repeated at both ends of an unrelated contig.
MAX_HOMOPOLYMER_FRAC = 0.60   # longest single-base run, as a fraction of repeat length
MIN_SHANNON_BITS = 1.20       # mononucleotide entropy, of 2.0 bits maximum
MIN_KMER_DIVERSITY = 0.50     # distinct 4-mers / (k - 3); kills (AT)n and similar
MIN_ACGT_FRAC = 0.90          # a repeat made of Ns is not an observation
_COMPLEXITY_K = 4


def env(name):
    v = os.environ.get(name)
    if not v:
        sys.exit(f"{name} unset — run: source config/config.sh")
    return v


# ---------------------------------------------------------------------------- circularity

def maximal_terminal_repeat(seq, max_k=None):
    """Length of the longest prefix that is also the suffix, ignoring the trivial whole
    sequence. Returns 0 if there is none."""
    s = str(seq).upper()
    n = len(s)
    if n < 2:
        return 0
    hi = n // 2 if max_k is None else min(max_k, n // 2)
    for k in range(hi, 0, -1):
        if s[:k] == s[-k:]:
            return k
    return 0


def _longest_run(s):
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best if s else 0


def _shannon_bits(s):
    counts = {}
    for c in s:
        counts[c] = counts.get(c, 0) + 1
    n = len(s)
    if n == 0:
        return 0.0
    import math
    return -sum((c / n) * math.log2(c / n) for c in counts.values())


def repeat_rejection(rep):
    """Why this terminal repeat is not usable evidence, or '' if it is.

    Order matters only for which reason gets reported first; any one of them disqualifies.
    """
    if not rep:
        return "no_repeat"
    r = rep.upper()
    acgt = sum(1 for c in r if c in "ACGT")
    if acgt / len(r) < MIN_ACGT_FRAC:
        return "non_acgt"
    if _longest_run(r) / len(r) >= MAX_HOMOPOLYMER_FRAC:
        return "homopolymer"
    if _shannon_bits(r) < MIN_SHANNON_BITS:
        return "low_entropy"
    if len(r) > _COMPLEXITY_K:
        kmers = {r[i:i + _COMPLEXITY_K] for i in range(len(r) - _COMPLEXITY_K + 1)}
        if len(kmers) / (len(r) - _COMPLEXITY_K + 1) < MIN_KMER_DIVERSITY:
            return "low_kmer_diversity"
    return ""


def circularity_call(seq, min_overlap):
    """(circular, repeat_len, rejection_reason). `circular` requires a terminal repeat of at
    least `min_overlap` bp that survives the complexity filter."""
    k = maximal_terminal_repeat(seq)
    if k == 0:
        return False, 0, "no_repeat"
    reason = repeat_rejection(str(seq).upper()[:k])
    if reason:
        return False, k, reason
    return (k >= min_overlap), k, ""


# ------------------------------------------------------------------------------------ ORFs

class Orf(object):
    __slots__ = ("aa", "strand", "start_nt", "spans_origin", "full_turn")

    def __init__(self, aa="", strand="", start_nt=-1, spans_origin=False, full_turn=False):
        self.aa = aa
        self.strand = strand
        self.start_nt = start_nt
        self.spans_origin = spans_origin
        self.full_turn = full_turn

    @property
    def aa_len(self):
        return len(self.aa)


def _translate(nt):
    """Trim to a whole number of codons before translating: a partial trailing codon makes
    Biopython warn, and it can never contribute a residue anyway."""
    usable = len(nt) - (len(nt) % 3)
    if usable <= 0:
        return ""
    return str(Seq(str(nt)[:usable]).translate(to_stop=False))


def _scan_frame(prot, frame, strand, mol_len, out):
    """Collect the best M-initiated ORF of every stop-free segment of one frame.

    `prot` is the translation of the DOUBLED molecule, so nucleotide offsets run 0..2L.
    Only ORFs whose start lies in the first copy are kept; a start at >= L is the same ORF
    found one turn later.
    """
    max_aa = mol_len // 3
    aa_pos = 0
    for seg in prot.split("*"):
        m = seg.find("M")
        if m >= 0:
            aa_start = aa_pos + m
            nt_start = frame + 3 * aa_start
            if nt_start < mol_len:
                aa = seg[m:]
                full_turn = len(aa) > max_aa
                if full_turn:
                    aa = aa[:max_aa]
                out.append(Orf(
                    aa=aa,
                    strand=strand,
                    start_nt=nt_start,
                    spans_origin=(nt_start + 3 * len(aa)) > mol_len,
                    full_turn=full_turn,
                ))
        aa_pos += len(seg) + 1
    return out


def longest_orf_circular(seq):
    """Longest ORF of a circular molecule, six frames, origin-crossing included."""
    s = str(seq).upper()
    L = len(s)
    if L < 3:
        return Orf()
    doubled = {"+": s + s, "-": str(Seq(s).reverse_complement()) * 2}
    found = []
    for strand, d in doubled.items():
        for frame in range(3):
            _scan_frame(_translate(d[frame:]), frame, strand, L, found)
    if not found:
        return Orf()
    # Deterministic tie-break: longest, then + strand, then earliest start.
    found.sort(key=lambda o: (-o.aa_len, o.strand != "+", o.start_nt))
    return found[0]


def longest_orf_linear(seq):
    """The plan's linear-only search, kept so the cost of C-12 can be measured on real data
    (and so tests can demonstrate the truncation). Not used for the reported ORF."""
    s = str(seq).upper()
    L = len(s)
    if L < 3:
        return Orf()
    found = []
    for strand, nt in (("+", s), ("-", str(Seq(s).reverse_complement()))):
        for frame in range(3):
            aa_pos = 0
            for seg in _translate(nt[frame:]).split("*"):
                m = seg.find("M")
                if m >= 0:
                    found.append(Orf(aa=seg[m:], strand=strand,
                                     start_nt=frame + 3 * (aa_pos + m)))
                aa_pos += len(seg) + 1
    if not found:
        return Orf()
    found.sort(key=lambda o: (-o.aa_len, o.strand != "+", o.start_nt))
    return found[0]


# ------------------------------------------------------------------------------------ main

COLUMNS = [
    # The first five are fixed by the interface contract; the rest are appended so the
    # heuristic can be audited rather than taken on faith.
    "id", "length", "circular", "orf_aa", "orf_spans_origin",
    "term_repeat_len", "term_repeat_reject", "orf_strand", "orf_start_nt", "orf_full_turn",
]


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", default=None, help="default $RESULTS/candidates.fna")
    p.add_argument("--tsv-out", default=None, help="default $RESULTS/circularity.tsv")
    p.add_argument("--faa-out", default=None, help="default $WORK/candidate_orfs.faa")
    p.add_argument("--min-overlap", type=int, default=None,
                   help="terminal repeat length required to call circular; default $CIRC_MIN_OVERLAP")
    p.add_argument("--min-orf-aa", type=int, default=100,
                   help="ORF length written to the .faa (plan value: 100)")
    p.add_argument("--force", action="store_true",
                   help="recompute even if outputs are newer than the input")
    return p.parse_args(argv)


def write_atomic(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write(text)
    os.replace(tmp, path)


def main(argv=None):
    args = parse_args(argv)
    results = env("RESULTS")
    work = env("WORK")
    fasta = args.fasta or os.path.join(results, "candidates.fna")
    tsv_out = args.tsv_out or os.path.join(results, "circularity.tsv")
    faa_out = args.faa_out or os.path.join(work, "candidate_orfs.faa")
    min_overlap = args.min_overlap if args.min_overlap is not None else int(env("CIRC_MIN_OVERLAP"))

    if not os.path.exists(fasta):
        sys.exit(f"missing {fasta} — run scripts/08_run_vnom.sh first (exit 2)" and 2)

    # Idempotent, but not stale: skip only when both outputs postdate the input.
    if not args.force and os.path.exists(tsv_out) and os.path.exists(faa_out):
        src = os.path.getmtime(fasta)
        if os.path.getmtime(tsv_out) >= src and os.path.getmtime(faa_out) >= src:
            print(f"circularity.tsv and candidate_orfs.faa are up to date — "
                  f"use --force to recompute")
            return 0

    rows, faa, seen = [], [], {}
    n_gained = 0
    gained_aa = 0
    for rec in SeqIO.parse(fasta, "fasta"):
        if rec.id in seen:
            # Duplicate ids break `seqkit grep -f` and would silently merge two different
            # elements in the replication matrix. See "Notes for downstream" in the report.
            sys.exit(f"duplicate candidate id {rec.id!r} in {fasta}: pooled candidate ids "
                     f"must be unique (scripts/08 must namespace them per accession)")
        seen[rec.id] = True
        if any(c.isspace() or c in "'\"" for c in rec.id):
            sys.exit(f"candidate id {rec.id!r} contains whitespace or a quote; "
                     f"$WORK/passed_ids.txt must stay consumable by `seqkit grep -f`")

        circ, rep_len, reason = circularity_call(rec.seq, min_overlap)
        orf = longest_orf_circular(rec.seq)
        lin = longest_orf_linear(rec.seq)
        if orf.aa_len > lin.aa_len:
            n_gained += 1
            gained_aa += orf.aa_len - lin.aa_len

        rows.append((rec.id, len(rec.seq), circ, orf.aa_len, orf.spans_origin,
                     rep_len, reason, orf.strand or "NA", orf.start_nt, orf.full_turn))
        if orf.aa_len >= args.min_orf_aa:
            faa.append(f">{rec.id}\n{orf.aa}\n")

    if not rows:
        sys.stderr.write(f"{fasta} contains 0 records — nothing to validate; "
                         f"re-run scripts/08_run_vnom.sh\n")
        return 2

    write_atomic(faa_out, "".join(faa))
    write_atomic(tsv_out,
                 "\t".join(COLUMNS) + "\n" +
                 "".join("\t".join(map(str, r)) + "\n" for r in rows))

    n_circ = sum(1 for r in rows if r[2])
    n_span = sum(1 for r in rows if r[4])
    print(f"candidates:            {len(rows)}")
    print(f"circular (heuristic):  {n_circ}")
    print(f"ORF>={args.min_orf_aa}aa:            {len(faa)}")
    print(f"ORFs crossing origin:  {n_span}   <- invisible to the plan's linear-only search")
    print(f"ORFs longer than the linear search would report: {n_gained} "
          f"(+{gained_aa} aa total)")
    print("terminal repeats rejected as low-complexity: "
          f"{sum(1 for r in rows if r[6] not in ('', 'no_repeat'))}")
    print("NOTE: a terminal repeat is supporting evidence only. Confirm circularity with "
          "scripts/12_junction_mapping.sh.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
