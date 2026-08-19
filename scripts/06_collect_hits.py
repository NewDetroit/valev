#!/usr/bin/env python3
"""Task 8 — collect triage hits, confirm with the profile HMM, evaluate GATE 2.

Reads  : $WORK/triage/<ACC>.tsv        per-accession DIAMOND output (from scripts/05)
         $WORK/run_metadata.tsv        Run -> BioProject/Center/Platform (from scripts/04)
Writes : $RESULTS/triage_hits.tsv      merged, with header
         $WORK/triage_best.faa         best hit per accession, for HMMER
         $WORK/triage_hmm.tbl          raw hmmsearch table
         $RESULTS/confirmed_hits.tsv   HMM-confirmed hits joined to run metadata
         $WORK/deep_accessions.txt     the deep path's input list

Exit: 0 GATE 2 passed · 2 missing input · 3 GATE 2 failed (a result, not a crash)

CORRECTIONS IMPLEMENTED HERE
  C-04  The plan wrote DIAMOND's `qseq` into a file named .faa and ran hmmsearch
        on it. In blastx `qseq` is nucleotide, so HMMER was being handed DNA to
        score against an amino-acid profile — the plan calls this step "a real
        second filter", and as written it filtered nothing. scripts/05 now emits
        `qseq_translated`; this script ASSERTS the collected sequences are protein
        before HMMER runs, and refuses to continue if they are not.
  C-06  The plan had 100 concurrent array tasks appending to one shared TSV.
        Merging happens here instead, single-threaded, after the sweep.
"""

from __future__ import annotations

import glob
import os
import re
import subprocess
import sys
from collections import defaultdict

COLS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qseq_translated"]
# Hit rate above which the plan says to stop. Far above biological plausibility.
MAX_HIT_RATE = 0.20
NT_ONLY = re.compile(r"^[ACGTUN]+$", re.IGNORECASE)
PROTEIN_OK = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYBXZJUO*\-]+$", re.IGNORECASE)


def die(msg: str, code: int = 1):
    """sys.exit(str) always exits 1; the exit-code convention needs explicit codes."""
    sys.stderr.write(msg.rstrip("\n") + "\n")
    raise SystemExit(code)


def env(name: str) -> str:
    v = os.environ.get(name)
    if not v:
        sys.exit(f"[06_collect_hits] {name} unset — run: source config/config.sh")
    return v


def main() -> int:
    WORK, RESULTS, REF = env("WORK"), env("RESULTS"), env("REF")
    hmm = env("REF_OBLIN_HMM")
    hmm_e = os.environ.get("HMM_EVALUE", "1e-10")
    hmm_e_relaxed = os.environ.get("HMM_EVALUE_RELAXED", "1e-5")
    min_bp = int(os.environ.get("MIN_INDEPENDENT_BIOPROJECTS", "3"))

    tdir = os.path.join(WORK, "triage")
    if not os.path.isdir(tdir):
        die(f"[06_collect_hits] no {tdir} — run scripts/05_logan_triage.sh first (exit 2)", 2)

    # ---- merge (C-06) -------------------------------------------------------
    done = sorted(glob.glob(os.path.join(tdir, "*.done")))
    if not done:
        die(f"[06_collect_hits] no completed accessions in {tdir} (exit 2)", 2)

    statuses = defaultdict(int)
    rows = []
    for d in done:
        acc = os.path.basename(d)[: -len(".done")]
        st_path = os.path.join(tdir, f"{acc}.status")
        st = open(st_path).read().strip() if os.path.exists(st_path) else "UNKNOWN"
        statuses[st] += 1
        tsv = os.path.join(tdir, f"{acc}.tsv")
        if st != "HIT" or not os.path.exists(tsv):
            continue
        for line in open(tsv):
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) != len(COLS):
                sys.stderr.write(
                    f"[06_collect_hits] malformed row in {tsv}: {len(f)} fields, "
                    f"expected {len(COLS)} — skipping\n")
                continue
            rows.append([acc] + f)

    n_done = len(done)
    n_searched = statuses["HIT"] + statuses["NONE"]
    hit_accs = sorted({r[0] for r in rows})

    out_hits = os.path.join(RESULTS, "triage_hits.tsv")
    with open(out_hits, "w") as fh:
        fh.write("accession\t" + "\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    print(f"accessions completed : {n_done:,}")
    for k in ("HIT", "NONE", "ABSENT", "FAIL", "UNKNOWN"):
        if statuses[k]:
            print(f"  {k:<8}           : {statuses[k]:,}")
    print(f"searchable (HIT+NONE): {n_searched:,}")
    print(f"hit rows             : {len(rows):,}")
    print(f"accessions with hits : {len(hit_accs):,}")
    print(f"wrote {out_hits}")

    if not hit_accs:
        print("\nGATE 2: zero triage hits.")
        print("  Either the niche genuinely lacks Oblin-like sequence — a real, reportable")
        print("  negative — or the threshold is too strict. Before concluding, re-run the")
        print(f"  sweep at DIAMOND_EVALUE={os.environ.get('DIAMOND_EVALUE_RELAXED','1e-3')}.")
        print("  Confirm GATE 1 passed first: results/positive_control_evidence.tsv")
        return 3

    # Hit rate sanity (plan Task 7 Step 7, implemented as a check rather than prose).
    if n_searched:
        rate = len(hit_accs) / n_searched
        print(f"hit rate             : {rate:.2%}")
        if rate > MAX_HIT_RATE:
            print(f"\nSTOP (exit 3): hit rate {rate:.1%} exceeds {MAX_HIT_RATE:.0%}.")
            print("  That is far above biological plausibility. Either DIAMOND_EVALUE is too")
            print("  permissive or the decoy calibration in scripts/02 did not actually run")
            print("  (see C-03). Check ref/calibration.tsv, tighten, and re-sweep.")
            return 3

    # ---- best hit per accession --------------------------------------------
    best: dict[str, list[str]] = {}
    for r in rows:
        acc, ev = r[0], float(r[5])
        if acc not in best or ev < float(best[acc][5]):
            best[acc] = r

    faa = os.path.join(WORK, "triage_best.faa")
    n_bad = 0
    with open(faa, "w") as fh:
        for acc, r in sorted(best.items()):
            seq = r[7].replace("-", "").strip().upper()
            # C-04: the assertion the plan needed and did not have.
            if not seq:
                n_bad += 1
                continue
            if NT_ONLY.match(seq) and len(seq) > 12:
                sys.stderr.write(
                    f"[06_collect_hits] FATAL: {acc} sequence looks like NUCLEOTIDE:\n"
                    f"  {seq[:60]}...\n"
                    f"  scripts/05 must emit `qseq_translated`, not `qseq`. In blastx,\n"
                    f"  DIAMOND's `qseq` is the aligned part of the nucleotide query.\n"
                    f"  Feeding it to hmmsearch scores DNA against a protein profile.\n"
                    f"  This is correction C-04. Re-run the sweep with the fixed worker.\n")
                return 1
            if not PROTEIN_OK.match(seq):
                n_bad += 1
                continue
            fh.write(f">{acc}\n{seq}\n")
    if n_bad:
        print(f"  skipped {n_bad} unusable sequences")
    n_faa = len(best) - n_bad
    print(f"wrote {faa} ({n_faa} sequences, verified amino-acid)")
    if n_faa == 0:
        die("[06_collect_hits] no usable protein sequences to confirm (exit 1)", 1)

    # ---- HMM confirmation ---------------------------------------------------
    if not os.path.exists(hmm):
        die(f"[06_collect_hits] missing {hmm} — run scripts/02_build_search_dbs.sh (exit 2)", 2)

    tbl = os.path.join(WORK, "triage_hmm.tbl")

    def hmmsearch(evalue: str) -> set[str]:
        subprocess.run(
            ["hmmsearch", "--tblout", tbl, "-E", evalue, hmm, faa],
            check=True, stdout=subprocess.DEVNULL)
        return {ln.split()[0] for ln in open(tbl) if ln.strip() and not ln.startswith("#")}

    confirmed = hmmsearch(hmm_e)
    print(f"HMM-confirmed (E<={hmm_e}) : {len(confirmed):,}")
    if not confirmed:
        confirmed = hmmsearch(hmm_e_relaxed)
        print(f"  retried at E<={hmm_e_relaxed}: {len(confirmed):,}")
        if not confirmed:
            print("\nGATE 2 FAILED (exit 3): no HMM-confirmed hits at either threshold.")
            print("  DIAMOND found candidates but the profile HMM rejects all of them, which")
            print("  usually means the DIAMOND hits are spurious. This is the filter doing")
            print("  its job. Treat as a negative result for this niche.")
            return 3

    # ---- join metadata, evaluate independence -------------------------------
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
    else:
        sys.stderr.write(
            f"[06_collect_hits] WARNING: no {meta_path}; BioProject independence cannot be\n"
            f"  assessed and GATE 2 cannot be evaluated properly. Run scripts/04.\n")

    extra = ["BioProject", "LibraryStrategy", "Platform", "CenterName", "ScientificName"]
    out_conf = os.path.join(RESULTS, "confirmed_hits.tsv")
    per_bp: dict[str, list[str]] = defaultdict(list)
    with open(out_conf, "w") as fh:
        fh.write("accession\t" + "\t".join(COLS) + "\t" + "\t".join(extra) + "\n")
        for acc in sorted(confirmed):
            r = best[acc]
            m = meta.get(acc, {})
            fh.write("\t".join(r) + "\t" + "\t".join(m.get(c, "NA") for c in extra) + "\n")
            per_bp[m.get("BioProject", "NA")].append(acc)
    print(f"wrote {out_conf}")

    deep = os.path.join(WORK, "deep_accessions.txt")
    with open(deep, "w") as fh:
        for acc in sorted(confirmed):
            fh.write(acc + "\n")
    print(f"wrote {deep} ({len(confirmed)} accessions for the deep path)")

    known = {bp: a for bp, a in per_bp.items() if bp not in ("NA", "")}
    print(f"\nindependent BioProjects with confirmed hits: {len(known)}")
    for bp, accs in sorted(known.items(), key=lambda x: -len(x[1]))[:10]:
        print(f"  {bp}: {len(accs)}")

    # ---- GATE 2 -------------------------------------------------------------
    print("\n" + "=" * 72)
    if len(known) >= min_bp:
        print(f"GATE 2 PASSED — hits in {len(known)} independent BioProjects (need {min_bp}).")
        print("Next: scripts/07_deep_assemble.sh via scripts/submit_array.sh")
        print("=" * 72)
        return 0
    if len(known) >= 1:
        print(f"GATE 2 WARNING (exit 3) — hits in only {len(known)} BioProject(s), need {min_bp}.")
        print("  Treat this as a contamination signal, not a discovery: an element seen in one")
        print("  project is indistinguishable from an artifact of that project. Broaden the")
        print("  niche (scripts/04) and re-sweep before going further.")
    else:
        print("GATE 2 FAILED (exit 3) — no BioProject metadata, so independence is unknown.")
        print("  Run scripts/04_select_accessions.sh to populate work/run_metadata.tsv.")
    print("=" * 72)
    return 3


if __name__ == "__main__":
    sys.exit(main())
