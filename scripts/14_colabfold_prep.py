#!/usr/bin/env python3
"""Task 14 Steps 1-2 — prepare Oblin-like ORFs for ColabFold, and validate what comes back.

  Usage: 14_colabfold_prep.py [--validate results/figures/colabfold/plddt.tsv]

PLAN.md leaves this step entirely manual ("Open ColabFold in Google Colab, paste each
sequence"). The prediction itself has to be manual — ColabFold is a notebook — but the
preparation and the checking do not, and the checking is where mistakes get onto posters.

Writes : $RESULTS/validated_oblins.faa            all Oblin-like ORFs
         $RESULTS/figures/colabfold/<id>.fasta    one per sequence, for pasting
         $RESULTS/figures/colabfold/batch.fasta   combined, for batch mode
         $RESULTS/figures/colabfold/plddt.tsv     TEMPLATE — fill in from ColabFold
         $RESULTS/figures/colabfold/foldseek.tsv  TEMPLATE — fill in from Foldseek

UNVERIFIED FIGURE. PLAN.md states Zheludev et al. reported a mean pLDDT of 83.8 for
Oblin-1 and uses it as a comparison point. That number could not be verified from the
build environment (docs/CORRECTIONS.md OPEN-C) and is NOT baked in here. The only
threshold applied is the generic structural-biology convention that pLDDT < 70 is low
confidence. Do not put 83.8 on a poster until you have read it in the paper yourself.
"""
from __future__ import annotations
import argparse, os, sys

LOW_PLDDT = 70.0


def env(name: str) -> str:
    v = os.environ.get(name)
    if not v:
        sys.stderr.write(f"[14_colabfold_prep] {name} unset — run: source config/config.sh\n")
        raise SystemExit(2)
    return v


def read_fasta(path):
    rid, buf = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if rid:
                yield rid, "".join(buf)
            rid, buf = line[1:].split()[0], []
        elif line:
            buf.append(line.strip())
    if rid:
        yield rid, "".join(buf)


def validate(path: str) -> int:
    if not os.path.exists(path):
        sys.stderr.write(f"[14_colabfold_prep] no {path} to validate (exit 2)\n")
        return 2
    rows, bad = [], 0
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 2 or f[1].strip() in ("", "FILL_IN"):
                continue
            try:
                v = float(f[1])
            except ValueError:
                continue
            rows.append((f[0], v))
            if v < LOW_PLDDT:
                bad += 1
    if not rows:
        sys.stderr.write(f"[14_colabfold_prep] {path} has no filled-in pLDDT values yet\n")
        return 2
    mean = sum(v for _, v in rows) / len(rows)
    print(f"models with pLDDT : {len(rows)}")
    print(f"mean pLDDT        : {mean:.1f}")
    print(f"below {LOW_PLDDT:.0f}         : {bad}")
    for i, v in sorted(rows, key=lambda x: x[1]):
        flag = "  LOW CONFIDENCE — report as such" if v < LOW_PLDDT else ""
        print(f"  {i:<44} {v:6.1f}{flag}")
    if bad:
        print(f"\n{bad} model(s) below {LOW_PLDDT:.0f}. Say so explicitly wherever they appear.")
        print("A low-pLDDT model is not evidence of a novel fold — it is evidence of an")
        print("uncertain prediction, and a judge who knows the difference will ask.")
    return 0


def main() -> int:
    RESULTS, WORK = env("RESULTS"), env("WORK")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--validate", nargs="?", const="", help="check a filled-in plddt.tsv")
    args = ap.parse_args()

    outdir = os.path.join(RESULTS, "figures", "colabfold")
    if args.validate is not None:
        return validate(args.validate or os.path.join(outdir, "plddt.tsv"))

    orfs = os.path.join(WORK, "candidate_orfs.faa")
    ids_file = os.path.join(WORK, "passed_ids.txt")
    for p in (orfs, ids_file):
        if not os.path.exists(p):
            sys.stderr.write(f"[14_colabfold_prep] missing {p} — run scripts/09 and scripts/11 (exit 2)\n")
            return 2

    keep = {ln.strip() for ln in open(ids_file) if ln.strip()}
    os.makedirs(outdir, exist_ok=True)
    kept = [(i, s) for i, s in read_fasta(orfs) if i in keep]
    if not kept:
        sys.stderr.write("[14_colabfold_prep] no validated ORFs matched passed_ids.txt (exit 3)\n")
        return 3

    with open(os.path.join(RESULTS, "validated_oblins.faa"), "w") as fh, \
         open(os.path.join(outdir, "batch.fasta"), "w") as bfh:
        for i, s in kept:
            fh.write(f">{i}\n{s}\n")
            bfh.write(f">{i}\n{s}\n")
            safe = i.replace("/", "_")
            with open(os.path.join(outdir, f"{safe}.fasta"), "w") as one:
                one.write(f">{i}\n{s}\n")

    for name, cols in [("plddt.tsv", ["id", "mean_plddt", "model_file", "notes"]),
                       ("foldseek.tsv", ["id", "rank", "target", "database",
                                         "evalue", "prob", "description"])]:
        p = os.path.join(outdir, name)
        if os.path.exists(p):
            continue
        with open(p, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for i, _ in kept:
                if name == "plddt.tsv":
                    fh.write(f"{i}\tFILL_IN\tFILL_IN\t\n")
                else:
                    for r in range(1, 6):
                        fh.write(f"{i}\t{r}\tFILL_IN\tFILL_IN\tFILL_IN\tFILL_IN\t\n")

    print(f"prepared {len(kept)} Oblin-like ORFs")
    print(f"  {RESULTS}/validated_oblins.faa")
    print(f"  {outdir}/  (per-sequence FASTAs, batch.fasta, plddt.tsv, foldseek.tsv)")
    print()
    print("Next, manually:")
    print("  1. ColabFold AlphaFold2.ipynb — run each sequence, save PDBs into")
    print(f"     {outdir}/ and record mean pLDDT in plddt.tsv")
    print("  2. Foldseek — search each PDB against AFDB and PDB; record the top 5 in")
    print("     foldseek.tsv, or use scripts/15_foldseek.sh if the CLI is installed")
    print(f"  3. python3 scripts/14_colabfold_prep.py --validate")
    print()
    print("Both Foldseek outcomes are informative: no significant hits argues a genuinely")
    print("novel fold; hits to known Oblin-1 confirm a true Obelisk relative. Record which,")
    print("with E-values — an empty result you cannot substantiate is not a finding.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
