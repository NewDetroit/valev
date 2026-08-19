# CLAUDE.md — working notes for agents in this repo

## What this is

`obelisk-hunt` — a computational pipeline that searches NCBI SRA data for novel
**Obelisk-like elements**: ~1 kb circular RNA agents with rod-like secondary structure that
encode a protein called Oblin-1. Target: ISEF 2027.

The deliverable is `results/validated.fna` plus the evidence that justifies calling its
contents real. A **negative** result is a legitimate deliverable here, but only because the
positive control at GATE 1 proves the pipeline can find what it is looking for.

## Architecture in one paragraph

Two paths joined at a triage point. The **fast path** streams pre-assembled contigs for
millions of accessions out of the public Logan S3 bucket, pipes them through DIAMOND without
touching disk, and keeps only accessions that plausibly contain an Oblin-like protein. The
**deep path** then runs the full published Zheludev et al. assembly pipeline
(`fasterq-dump → fastp → rnaSPAdes → VNom`) on just those survivors — hundreds of runs
instead of millions. Three gates sit between the stages and each one can stop the project.

## Before you change anything

1. **`config/config.sh` is the single source of truth.** Every path and threshold lives
   there. Never hardcode one. Every script starts by sourcing it.
2. **Read `docs/CORRECTIONS.md`.** The original plan (`docs/PLAN.original.md`, preserved
   verbatim) contains 13 defects, several of which silently produce wrong science. They are
   fixed in the code and each fix is numbered `C-01`..`C-13`. If code looks like it disagrees
   with the plan, it is probably deliberate — check there before "fixing" it back.
3. **`docs/LEDGER.md`** records what has been built and reviewed.

## Three facts that are easy to get wrong

- **Logan's current release is v1.2** (SRA cutoff 31 Dec 2025, 38,124,741 accessions), not the
  v1.0 December-2023 release the original plan assumed. `s3://logan-pub/c/` is v1.2;
  `s3://logan-pub/c1.0/` is the old one and is slated for deletion. Switch with
  `LOGAN_RELEASE=v1.0`. See `C-01`.
- **The Obelisk authors already ran a Logan-wide Oblin sweep** and published it at
  `s3://logan-pub/paper/Obelisk/`. We consume their reference artifacts rather than rebuilding
  from a journal supplement. Their sweep was on v1.0, so the ~10.9M accessions added in
  2024–2025 are the genuinely unscreened part. See `C-02`.
- **`diamond blastx` `qseq` is nucleotide, not protein.** Use `qseq_translated` anywhere the
  output feeds HMMER. See `C-04`.

## Scientific constants — do not touch

These are the published parameters. Changing them forfeits the replication claim.

```
fastp   --average_qual=30 --n_base_limit=0 --cut_front --cut_tail
VNom    -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1
Minia3  k=31                       (Logan's assembler; informational)
```

## Conventions

- Bash: `#!/usr/bin/env bash`, `set -euo pipefail`. Per-accession workers use `set -uo pipefail`
  so one failure does not kill a batch — always commented where used.
- Every script is idempotent and guards on file **existence** (`-e`) or an explicit sentinel,
  never on non-emptiness (`-s`). An empty output is a valid "ran, found nothing". See `C-05`.
- No two concurrent processes write the same file. Per-task outputs are merged by a later
  collection step. See `C-06`.
- Scripts that can produce an empty or partial result exit non-zero and explain, rather than
  printing advice and continuing.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | success |
| 1 | unexpected error |
| 2 | a required input is missing — run an earlier script first |
| 3 | a GATE failed its criterion. A scientific result, not a crash. Stop and read the output. |

## Environment notes

- `env/environment.yml` builds the conda environment. `diamond=2.1.9` is pinned because
  `qseq_translated` needs ≥2.0.8.
- Logan S3 needs no credentials: `aws s3 cp ... --no-sign-request`, or plain HTTPS against
  `https://s3.amazonaws.com/logan-pub/...`.
- Some environments block `eutils.ncbi.nlm.nih.gov`. `scripts/04b_select_accessions_offline.sh`
  is the NCBI-free path for accession enumeration. See `OPEN-A`.

## When claiming a number

`docs/CORRECTIONS.md` ends with an audited provenance table. Anything not in it — including
several figures quoted in the original plan — is **unverified**; `OPEN-C` lists exactly which.
Run `scripts/15_audit_claims.sh` to refresh what can be checked automatically. Do not put an
unverified figure on a poster.
