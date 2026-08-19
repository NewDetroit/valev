#!/usr/bin/env python3
"""Task 12 — the replication matrix. GATE 3, and the project's primary contamination control.

The argument that a candidate is a real biological entity rather than an artifact is that
it recurs in data that share no laboratory, no library prep, and no sequencing run. A
contaminant, an index-hop, or an assembly chimera does not do that. This script measures it.

Reads  : $RESULTS/candidates_novel.fna  (from scripts/10_exclude_known.sh)
         $WORK/vnom/<ACC>/*.fasta       (from scripts/08_run_vnom.sh)
         $WORK/run_metadata.tsv         (from scripts/04)
Writes : $RESULTS/replication_matrix.tsv
         $WORK/passed_ids.txt           plain ids, one per line — consumed by `seqkit grep -f`
         $RESULTS/validated.fna         written by the caller after this passes

Exit: 0 GATE 3 passed · 2 missing input · 3 GATE 3 not passed (a result, not a crash)

CORRECTION C-08 — the plan reads only the first VNom FASTA per accession:

    fas = [f for f in os.listdir(path) if f.endswith(".fasta")]
    if not fas: continue
    out = subprocess.run(["blastn","-query",f"{path}/{fas[0]}", ...

`fas[0]` from an unsorted os.listdir discards every other FASTA VNom produced for that
accession. Occurrences are undercounted, so n_bioprojects is undercounted, so genuine
candidates fail the `>= MIN_INDEPENDENT_BIOPROJECTS` test — nondeterministically, since
listdir order is not stable. On the step the plan calls "the core anti-contamination
control", that silently rejects real discoveries. Here every FASTA is read, in sorted order.
"""

from __future__ import annotations

import os
import subprocess
import sys
from collections import defaultdict


def die(msg: str, code: int = 1):
    sys.stderr.write(msg.rstrip("\n") + "\n")
    raise SystemExit(code)


def env(name: str) -> str:
    v = os.environ.get(name)
    if not v:
        die(f"[11_replication_matrix] {name} unset — run: source config/config.sh", 2)
    return v


def need(tool: str):
    from shutil import which
    if which(tool) is None:
        die(f"[11_replication_matrix] missing tool: {tool} — conda activate obelisk", 2)


def main() -> int:
    WORK, RESULTS = env("WORK"), env("RESULTS")
    min_bp = int(os.environ.get("MIN_INDEPENDENT_BIOPROJECTS", "3"))
    ident = os.environ.get("KNOWN_IDENTITY_PCT", "95")

    novel = os.path.join(RESULTS, "candidates_novel.fna")
    if not os.path.exists(novel):
        die(f"[11_replication_matrix] missing {novel} — run scripts/10_exclude_known.sh (exit 2)", 2)

    vnom_root = os.path.join(WORK, "vnom")
    if not os.path.isdir(vnom_root):
        die(f"[11_replication_matrix] missing {vnom_root} — run scripts/08_run_vnom.sh (exit 2)", 2)

    # GATE 3's negative branch is only meaningful if GATE 1 passed. Say so loudly.
    pc = os.path.join(RESULTS, "positive_control_evidence.tsv")
    pc_ok = os.path.exists(pc)
    if not pc_ok:
        sys.stderr.write(
            "\n*** WARNING: no results/positive_control_evidence.tsv ***\n"
            "  GATE 1 has not been recorded as passed. A negative result from this script\n"
            "  is only publishable because the pipeline is known to find Obelisks where\n"
            "  they exist. Without that, 'we found nothing' is uninterpretable.\n"
            "  Run scripts/03_positive_control.sh before relying on anything below.\n\n")

    need("makeblastdb")
    need("blastn")

    db = os.path.join(WORK, "novel_db")
    subprocess.run(["makeblastdb", "-in", novel, "-dbtype", "nucl", "-out", db],
                   check=True, stdout=subprocess.DEVNULL)

    # ---- occurrence scan (C-08) --------------------------------------------
    occ: dict[str, set[str]] = defaultdict(set)
    n_fastas = 0
    accs = sorted(d for d in os.listdir(vnom_root)
                  if os.path.isdir(os.path.join(vnom_root, d)))
    if not accs:
        die(f"[11_replication_matrix] no accession directories under {vnom_root} (exit 2)", 2)

    for acc in accs:
        path = os.path.join(vnom_root, acc)
        fastas = sorted(f for f in os.listdir(path)
                        if f.endswith((".fasta", ".fa", ".fna")))
        if not fastas:
            continue
        # C-08: every file, not fastas[0].
        for fa in fastas:
            n_fastas += 1
            r = subprocess.run(
                ["blastn", "-query", os.path.join(path, fa), "-db", db,
                 "-outfmt", "6 qseqid sseqid pident", "-evalue", "1e-20",
                 "-perc_identity", str(ident)],
                capture_output=True, text=True)
            if r.returncode != 0:
                sys.stderr.write(f"  blastn failed on {acc}/{fa}: {r.stderr.strip()[:200]}\n")
                continue
            for line in r.stdout.splitlines():
                if line.strip():
                    occ[line.split("\t")[1]].add(acc)

    print(f"accessions scanned : {len(accs):,}")
    print(f"VNom FASTAs read   : {n_fastas:,}   (the plan would have read {len(accs):,} — C-08)")
    print(f"candidates observed: {len(occ):,}")

    # ---- metadata join ------------------------------------------------------
    meta_path = os.path.join(WORK, "run_metadata.tsv")
    meta: dict[str, dict[str, str]] = {}
    if os.path.exists(meta_path):
        with open(meta_path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for ln in fh:
                f = ln.rstrip("\n").split("\t")
                if len(f) == len(hdr):
                    d = dict(zip(hdr, f))
                    meta[d.get("Run", "")] = d
    if not meta:
        die(f"[11_replication_matrix] missing or empty {meta_path}.\n"
            f"  Independence cannot be assessed without BioProject metadata, and\n"
            f"  independence is the entire contamination argument. Run scripts/04. (exit 2)", 2)

    rows = []
    for cand, cand_accs in occ.items():
        bps, cents, plats = set(), set(), set()
        unknown = 0
        for a in sorted(cand_accs):
            m = meta.get(a)
            if not m:
                unknown += 1
                continue
            for key, acc_set in (("BioProject", bps), ("CenterName", cents), ("Platform", plats)):
                v = (m.get(key) or "").strip()
                if v and v.upper() not in ("NA", "NAN", "NONE"):
                    acc_set.add(v)
        rows.append({
            "candidate": cand,
            "n_accessions": len(cand_accs),
            "n_bioprojects": len(bps),
            "n_centers": len(cents),
            "n_platforms": len(plats),
            "n_accessions_no_metadata": unknown,
            "bioprojects": ";".join(sorted(bps)) or "NA",
            "accessions": ";".join(sorted(cand_accs)),
        })

    rows.sort(key=lambda r: (-r["n_bioprojects"], -r["n_centers"], -r["n_accessions"]))
    out = os.path.join(RESULTS, "replication_matrix.tsv")
    hdr = ["candidate", "n_accessions", "n_bioprojects", "n_centers", "n_platforms",
           "n_accessions_no_metadata", "bioprojects", "accessions"]
    with open(out, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in hdr) + "\n")
    print(f"wrote {out}")

    passed = [r for r in rows if r["n_bioprojects"] >= min_bp]
    ids = os.path.join(WORK, "passed_ids.txt")
    # Plain one-id-per-line, no header, no quoting — `seqkit grep -f` consumes this.
    with open(ids, "w") as fh:
        for r in passed:
            fh.write(r["candidate"] + "\n")
    print(f"wrote {ids}")

    print(f"\ncandidates examined                    : {len(rows):,}")
    print(f"passing >= {min_bp} independent BioProjects : {len(passed):,}")
    if rows:
        print("\ntop by independence:")
        print(f"  {'candidate':<40} {'acc':>5} {'BP':>4} {'ctr':>4} {'plat':>5}")
        for r in rows[:20]:
            print(f"  {r['candidate'][:40]:<40} {r['n_accessions']:>5} "
                  f"{r['n_bioprojects']:>4} {r['n_centers']:>4} {r['n_platforms']:>5}")

    # ---- GATE 3 -------------------------------------------------------------
    print("\n" + "=" * 72)
    if passed:
        print(f"GATE 3 PASSED — {len(passed)} candidate(s) replicated across >= {min_bp} "
              f"independent BioProjects.")
        print("  Still required before claiming a discovery:")
        print("    scripts/12_junction_mapping.sh   physical evidence of circularity")
        print("    contaminant screen               no vector/PhiX/adapter hits")
        print("  Emit the validated set with:")
        print(f"    seqkit grep -f {ids} {novel} > {RESULTS}/validated.fna")
        print("=" * 72)
        return 0

    best = rows[0]["n_bioprojects"] if rows else 0
    if best >= 1:
        print(f"GATE 3 NOT PASSED (exit 3) — best candidate spans {best} BioProject(s), "
              f"need {min_bp}.")
        print("  Report these as TENTATIVE and flag them explicitly. An element seen in one")
        print("  or two projects cannot be distinguished from an artifact of those projects.")
        print("  Do not overclaim. Broadening the niche and re-sweeping is the honest next step.")
    else:
        print("GATE 3 NOT PASSED (exit 3) — no candidate replicated in any BioProject.")
        print()
        print("  This is a PUBLISHABLE NEGATIVE RESULT, stated as:")
        print('    "A systematic search of N accessions across M BioProjects in [niche],')
        print('     using a pipeline validated on a positive control, found no novel')
        print('     Obelisk-like elements."')
        print()
        if pc_ok:
            print("  GATE 1 passed, which is what makes that statement a finding rather than")
            print("  an admission of failure. This is a legitimate ISEF project — say so plainly.")
        else:
            print("  WARNING: GATE 1 is NOT recorded as passed, so you cannot make that claim")
            print("  yet. Run scripts/03_positive_control.sh first. Without it, a null result")
            print("  is indistinguishable from a broken pipeline.")
    print("=" * 72)
    return 3


if __name__ == "__main__":
    sys.exit(main())
