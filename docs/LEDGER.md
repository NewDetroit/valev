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
| T1 | Phase 0 — Logan verification | `scripts/00_verify_logan.sh` | dispatched (A, sonnet) |
| T2 | Phase 0 — environment | `env/environment.yml`, `docs/RUNBOOK.md` | dispatched (A, sonnet) |
| T3 | Phase 1 — reference acquisition | `scripts/01_fetch_references.sh`, `scripts/tsv_to_fasta.py` | dispatched (B, opus) |
| T4 | Phase 1 — search databases | `scripts/02_build_search_dbs.sh` | dispatched (B, opus) |
| T5 | Phase 2 — positive control (GATE 1) | `scripts/03_positive_control.sh` | dispatched (C, sonnet) |
| T6 | Phase 3 — niche selection | `config/niche_query.txt`, `scripts/04_select_accessions.sh`, `scripts/04b_select_accessions_offline.sh` | dispatched (D, opus) |
| T7 | Phase 3 — triage sweep | `scripts/05_logan_triage.sh`, `scripts/05_logan_triage.sbatch`, `scripts/submit_array.sh` | dispatched (E, opus) |
| T8 | Phase 3 — HMM confirmation (GATE 2) | `scripts/06_collect_hits.py` | dispatched (E, opus) |
| T9 | Phase 4 — deep assembly | `scripts/07_deep_assemble.sh`, `scripts/07_deep_assemble.sbatch` | dispatched (F, sonnet) |
| T10 | Phase 4 — VNom | `scripts/08_run_vnom.sh` | dispatched (F, sonnet) |
| T11 | Phase 5 — circularity + ORFs | `scripts/09_circularity_check.py` | dispatched (G, opus) |
| T12 | Phase 5 — replication matrix (GATE 3) | `scripts/11_replication_matrix.py` | dispatched (G, opus) |
| T13 | Phase 6 — RNA structure | `scripts/12_rnafold.sh` | dispatched (H, sonnet) |
| T14 | Phase 6 — protein structure | `scripts/13_colabfold_prep.py`, `scripts/14_foldseek.sh`, `scripts/10_cluster.sh` | dispatched (H, sonnet) |
| T15 | Phase 7 — wet-lab handoff | `scripts/15_wetlab_package.py`, `scripts/15_audit_claims.sh` | dispatched (H, sonnet) |
| T16 | Top-level docs | `PLAN.md`, `README.md`, `CLAUDE.md` | in progress (controller) |

---

## Completion log

Appended as each task's review comes back clean. One line per task.
