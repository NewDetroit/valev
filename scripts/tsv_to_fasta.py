#!/usr/bin/env python3
"""Extract a sequence column from a delimited table into FASTA.

Used only by the ``--from-scratch`` branch of ``scripts/01_fetch_references.sh``,
which rebuilds the reference set from the Zheludev et al. *Cell* supplementary
table instead of the published ``s3://logan-pub/paper/Obelisk/`` artifacts
(CORRECTION C-02).  The primary path never calls this.

Positional arguments are the plan's, in the plan's order, so an operator
following PLAN.md verbatim gets the same behaviour:

    tsv_to_fasta.py <table.tsv> <SEQ_COLUMN> <ID_COLUMN> <out.fasta>

The plan then told the operator to eyeball ``seqkit stats`` output and, if the
mean length was "far from ~1,000", to go back and pick a different column.
That check is encoded here instead: a column that is not sequence fails the
alphabet assertion and the script exits non-zero.  Guessing wrong about which
column holds the sequence is the single failure mode that would silently
invalidate every downstream search, so it is not left to the eye.

stdlib only: this runs before the conda environment is guaranteed to exist, and
``csv`` round-trips raw sequence strings without the dtype inference and NaN
coercion that ``pandas.read_csv`` applies to a column of A/C/G/T text.

Exit codes follow docs/INTERFACES: 0 ok, 1 unexpected, 2 required input
missing or malformed, 3 a content assertion failed.
"""

import argparse
import collections
import csv
import gzip
import io
import os
import statistics
import sys

# IUPAC nucleotide + gap/unknown.  'U' is accepted so an RNA-alphabet table is
# recognised as nucleotide rather than rejected; --rna-to-dna converts it.
NT_ALPHABET = set("ACGTUNRYSWKMBDHV-*.")
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYBZXUO-*.")
# Letters that are amino-acid-only.  Their absence over a large sample is what
# distinguishes a protein column from a nucleotide column, because every
# nucleotide letter is also a legal amino-acid letter.
AA_ONLY = set("EFILPQZJ")

NULLS = {"", "NA", "N/A", "NAN", "NONE", "NULL", "-", "."}


def die(msg, code):
    print(f"tsv_to_fasta: {msg}", file=sys.stderr)
    raise SystemExit(code)


def opener(path):
    """Supplementary tables are routinely distributed gzipped."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace", newline="")


def sniff_delimiter(path, override):
    if override:
        return {"tab": "\t", "comma": ",", "semicolon": ";"}[override]
    with opener(path) as fh:
        head = fh.readline()
    # A journal supplementary "tsv" is sometimes a csv.  Pick whichever
    # delimiter actually splits the header into more than one field.
    return "\t" if head.count("\t") >= head.count(",") else ","


def wrap(seq, width):
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main(argv):
    ap = argparse.ArgumentParser(
        prog="tsv_to_fasta.py",
        description="Extract a sequence column from a table into FASTA.",
    )
    ap.add_argument("table")
    ap.add_argument("seq_col")
    ap.add_argument("id_col")
    ap.add_argument("out")
    ap.add_argument("--alphabet", choices=("nt", "aa", "any"), default="nt",
                    help="assert the extracted residues are of this type (default: nt)")
    ap.add_argument("--delimiter", choices=("tab", "comma", "semicolon"), default=None)
    ap.add_argument("--wrap", type=int, default=60)
    ap.add_argument("--min-len", type=int, default=1,
                    help="skip records shorter than this")
    ap.add_argument("--rna-to-dna", action="store_true",
                    help="rewrite U->T (blastn and makeblastdb expect DNA)")
    ap.add_argument("--allow-duplicate-ids", action="store_true")
    args = ap.parse_args(argv)

    if not os.path.exists(args.table):
        die(f"no such table: {args.table}", 2)

    delim = sniff_delimiter(args.table, args.delimiter)

    # QUOTE_NONE: sequence tables are unquoted, and honouring quotes would let a
    # stray '"' in a description field swallow the rest of the row.
    with opener(args.table) as fh:
        rdr = csv.DictReader(fh, delimiter=delim, quoting=csv.QUOTE_NONE)
        if rdr.fieldnames is None:
            die(f"{args.table} is empty", 2)
        cols = [c.strip() for c in rdr.fieldnames]
        for want in (args.seq_col, args.id_col):
            if want not in cols:
                print("tsv_to_fasta: available columns:", file=sys.stderr)
                for i, c in enumerate(cols, 1):
                    print(f"  {i:3d}  {c}", file=sys.stderr)
                die(f"column {want!r} not present in {args.table}", 2)

        seen = {}
        lengths = []
        residues = collections.Counter()
        n_null = 0
        n_short = 0
        n_renamed = 0
        rows = 0

        tmp = args.out + ".part"
        with open(tmp, "w", encoding="utf-8") as out:
            for raw in rdr:
                rows += 1
                row = {(k.strip() if k else k): v for k, v in raw.items()}
                s = (row.get(args.seq_col) or "")
                s = "".join(s.split()).upper()
                if s in NULLS:
                    n_null += 1
                    continue
                if args.rna_to_dna:
                    s = s.replace("U", "T")
                if len(s) < args.min_len:
                    n_short += 1
                    continue

                rid = (row.get(args.id_col) or "").strip()
                if not rid:
                    die(f"row {rows}: empty id in column {args.id_col!r}", 3)
                # A FASTA id ends at the first space.  seqkit grep -f, DIAMOND's
                # sseqid and hmmsearch's target name all key on that token, so an
                # id containing whitespace would silently truncate downstream.
                if any(ch.isspace() for ch in rid):
                    rid = "_".join(rid.split())
                    n_renamed += 1
                if rid in seen and not args.allow_duplicate_ids:
                    die(f"duplicate id {rid!r} (rows {seen[rid]} and {rows}); "
                        f"pick a unique id column or pass --allow-duplicate-ids", 3)
                seen[rid] = rows

                lengths.append(len(s))
                residues.update(s)
                out.write(f">{rid}\n{wrap(s, args.wrap)}\n")

    n = len(lengths)
    if n == 0:
        os.unlink(tmp)
        die(f"0 sequences extracted from {rows} rows — column {args.seq_col!r} "
            f"holds no usable sequence", 3)

    total = sum(residues.values())
    frac_nt = sum(v for k, v in residues.items() if k in NT_ALPHABET) / total
    frac_aa = sum(v for k, v in residues.items() if k in AA_ALPHABET) / total
    frac_aa_only = sum(v for k, v in residues.items() if k in AA_ONLY) / total

    print(f"wrote {n} sequences to {args.out}")
    print(f"  rows read      : {rows}  (null/blank {n_null}, "
          f"shorter than {args.min_len} {n_short}, ids despaced {n_renamed})")
    print(f"  length         : min {min(lengths)}  max {max(lengths)}  "
          f"mean {statistics.mean(lengths):.1f}  median {statistics.median(lengths):.1f}")
    print(f"  residues       : {total}  nt-alphabet {frac_nt:.4f}  "
          f"aa-alphabet {frac_aa:.4f}  aa-only-letters {frac_aa_only:.4f}")
    top = "".join(k for k, _ in residues.most_common(8))
    print(f"  commonest      : {top}")

    # Nucleotide letters are a subset of amino-acid letters, so "looks like
    # protein" has to be decided on the presence of aa-only letters, not on the
    # aa fraction, which is ~1.0 for DNA too.
    if args.alphabet == "nt":
        if frac_nt < 0.99 or frac_aa_only > 0.01:
            os.unlink(tmp)
            die(f"column {args.seq_col!r} is not nucleotide "
                f"(nt-alphabet {frac_nt:.4f}, aa-only letters {frac_aa_only:.4f}). "
                f"Re-run against the column that holds the sequence.", 3)
    elif args.alphabet == "aa":
        if frac_aa < 0.99 or frac_aa_only < 0.01:
            os.unlink(tmp)
            die(f"column {args.seq_col!r} is not protein "
                f"(aa-alphabet {frac_aa:.4f}, aa-only letters {frac_aa_only:.4f}). "
                f"Re-run against the column that holds the sequence.", 3)

    os.replace(tmp, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
