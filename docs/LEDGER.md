# EXECUTION LEDGER

Durable progress record for the subagent-driven build of the obelisk-hunt pipeline.
**This file, not conversation memory, is the source of truth for what is done.**
A task marked `complete` here must not be re-dispatched.

- **Branch:** `claude/skills-verification-i6l0ue`
- **Scope decided with user:** scaffold + live verification. Pipeline execution
  (SLURM, conda toolchain, NCBI) is out of scope in this container — see `docs/RUNBOOK.md`.
- **Governing decisions:** Logan **v1.2 primary, v1.0 pinnable**; reference set from the
  **published `paper/Obelisk/` artifacts with from-scratch fallback**; **correctness governs**
  over plan-verbatim, every deviation logged in `docs/CORRECTIONS.md`.

---

## Pre-flight

| Item | Status |
|---|---|
| Plan read end to end | complete |
| Conflict scan, batched to user | complete — 3 decisions returned |
| Ground truth re-verified against primary sources | complete — 13 corrections, 3 open items |
| `config/config.sh` written and sourced in all 3 modes | complete |
| `docs/CORRECTIONS.md` | complete |
| `docs/PLAN.original.md` preserved verbatim | complete |

---

## Task ledger

Status values: `pending` → `dispatched` → `implemented` → `reviewed` → `complete`.

| # | Task | Files | Status |
|---|---|---|---|
| T1 | Phase 0 — Logan verification | `scripts/00_verify_logan.sh` | complete |
| T2 | Phase 0 — environment | `env/environment.yml`, `docs/RUNBOOK.md` | complete |
| T3 | Phase 1 — reference acquisition | `scripts/01_fetch_references.sh`, `scripts/tsv_to_fasta.py` | complete |
| T4 | Phase 1 — search databases | `scripts/02_build_search_dbs.sh` | complete |
| T5 | Phase 2 — positive control (GATE 1) | `scripts/03_positive_control.sh` | complete |
| T6 | Phase 3 — niche selection | `config/niche_query.txt`, `scripts/04_select_accessions.sh`, `scripts/04b_select_accessions_offline.sh` | complete |
| T7 | Phase 3 — triage sweep | `scripts/05_logan_triage.sh`, `scripts/05_logan_triage.sbatch`, `scripts/submit_array.sh` | complete |
| T8 | Phase 3 — HMM confirmation (GATE 2) | `scripts/06_collect_hits.py` | complete |
| T9 | Phase 4 — deep assembly | `scripts/07_deep_assemble.sh`, `scripts/07_deep_assemble.sbatch` | complete |
| T10 | Phase 4 — VNom | `scripts/08_run_vnom.sh` | complete |
| T11 | Phase 5 — circularity + ORFs | `scripts/09_circularity_check.py` | complete |
| T12 | Phase 5 — replication matrix (GATE 3) | `scripts/11_replication_matrix.py` | complete |
| T13 | Phase 6 — RNA structure | `scripts/13_rnafold.sh` | complete |
| T14 | Phase 6 — protein structure | `scripts/14_colabfold_prep.py`, `scripts/15_foldseek.sh`, `scripts/16_cluster.sh` | complete |
| T15 | Phase 7 — wet-lab handoff | `scripts/17_wetlab_package.py`, `scripts/18_audit_claims.sh` | complete |
| T16 | Top-level docs | `PLAN.md`, `README.md`, `CLAUDE.md` | complete |

---

## Completion log

All 16 tasks complete. Eight subagents were dispatched per the chosen strategy; all eight
were terminated mid-flight by a session limit, leaving four partial files and no reports.
Execution continued inline from the controller's context rather than re-dispatching, since
a fresh wave would have hit the same limit and re-derived context already held.

Of the four surviving files, `estimate_sweep_cost.py` needed its data-access layer rewritten
(`ParquetFileFormat.make_fragment()` rejects a raw file object; a corrupted expression sat in
`seek()`), which is why unreviewed agent output was never committed.

Verified by execution, not inspection:
- `00_verify_logan.sh` — 14 checks, 0 failures, against live S3; failure path confirmed
- `estimate_sweep_cost.py` — 38,124,741 rows over HTTP Range, 339.9 MB not 1.8 GB; exits 2 on a bad column
- `09_circularity_check.py` — origin-spanning ORF found, homopolymer ends rejected (tests/test_circularity.sh)
- `05_logan_triage.sh` — ABSENT/skip/retry semantics all correct
- `06_collect_hits.py` — nucleotide input rejected, exit-code convention holds
- `submit_array.sh` — 250 accessions at chunk 100 gives 1-3%100
- `12_junction_mapping.sh` — CIGAR span parser correct incl. insertions
- `18_audit_claims.sh` — 9 passed, 0 failed, live
- Full syntax/compile sweep across all 25 scripts and tests

Two corrections were found by resolving the plan's own open unknowns: C-14 (every circuclust
flag wrong) and C-15 (the VNom invocation cannot work as written).
