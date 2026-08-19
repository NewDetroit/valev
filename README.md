# obelisk-hunt

A computational pipeline for discovering novel **Obelisk-like elements** in an under-mined
sequencing niche.

Obelisks are ~1 kb circular RNA agents with rod-like secondary structure that encode a protein
family called Oblin-1. They were described in 2024 and are found across the microbiome, but
most of the public sequence archive has never been searched for them systematically.

**Project:** ISEF 2027 · **Status:** pipeline built and verified; not yet executed.

---

## The design decision

Two paths, joined at a triage point.

```
                            ┌──────────────────────────────────┐
                            │  Logan S3  ·  38,124,741 SRA     │
                            │  accessions, pre-assembled       │
                            └───────────────┬──────────────────┘
                                            │  stream, never store
                                            ▼
   FAST PATH        aws s3 cp - │ zstdcat │ diamond blastx  ──►  cheap, whole-archive
                                            │
                                            ▼
                                    ╔═══════════════╗
                                    ║    GATE 2     ║  Oblin-like protein present?
                                    ╚═══════╤═══════╝  survivors: tens–hundreds
                                            ▼
   DEEP PATH   fasterq-dump → fastp → rnaSPAdes → VNom   ──►  expensive, published method
                                            │
                                            ▼
                                    ╔═══════════════╗
                                    ║    GATE 3     ║  replicated across ≥3 BioProjects?
                                    ╚═══════╤═══════╝
                                            ▼
                                    results/validated.fna
```

Zheludev et al. ran the deep path on a curated sample set. Logan lets you pre-filter the whole
archive first, so the expensive step runs a few hundred times instead of a few million.

**One thing to be clear about.** Screening Logan for Oblins is not an unclaimed idea — the
Obelisk authors did it and published their results at `s3://logan-pub/paper/Obelisk/`. Their
sweep ran on Logan v1.0 (accessions up to December 2023). What is genuinely unscreened is the
**10,855,431 accessions added in v1.2** (2024–2025), plus a per-niche deep-assembly pass that
an archive-wide HMM scan does not do. See [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md) `C-02`.

---

## Quick start

```bash
git clone <this repo> obelisk-hunt && cd obelisk-hunt
conda env create -f env/environment.yml && conda activate obelisk
source config/config.sh

bash scripts/00_verify_logan.sh          # is Logan live, and is our ground truth still true?
bash scripts/01_fetch_references.sh      # reference set from the published Obelisk artifacts
bash scripts/02_build_search_dbs.sh      # DIAMOND db, profile HMM, calibrated decoy control
bash scripts/03_positive_control.sh      # GATE 1 — do not proceed until this passes
```

Full operator guide, including what must be run on a machine with NCBI access:
[`docs/RUNBOOK.md`](docs/RUNBOOK.md).

---

## The three gates

Each gate can stop the project, and each failure has a defined next action rather than a dead
end.

| Gate | Where | Criterion | If it fails |
|---|---|---|---|
| **1** | `scripts/03_positive_control.sh` | Recovers the known Obelisk-S.s from *S. sanguinis* SK36 RNA-seq | **Stop.** Every downstream result is meaningless until the pipeline can find a known answer. |
| **2** | `scripts/06_collect_hits.py` | HMM-confirmed hits in ≥3 independent BioProjects | 1–2 BioProjects is a contamination warning, not a discovery. Broaden the niche and re-sweep. |
| **3** | `scripts/11_replication_matrix.py` | ≥1 candidate replicated across ≥3 BioProjects, circular by junction mapping, Oblin-1 HMM positive, no contaminant hits | Report as tentative, or — if nothing passes — as a **publishable negative result**. |

The negative result is only meaningful *because* GATE 1 passed. That dependency is the most
important structural property of this design: "we searched N accessions with a pipeline proven
to find Obelisks, and found none" is a finding. "We found none" on its own is not.

---

## Layout

```
config/config.sh          single source of truth — all paths and thresholds
config/niche_query.txt    the Entrez query defining the target niche
env/environment.yml       conda spec
ref/                      built once by scripts/01–02, never edited by hand
scripts/                  00 → 18, run in order
work/                     scratch, safe to delete
results/                  the only directory that matters for the poster
docs/PLAN.original.md     the original plan, verbatim and unmodified
docs/CORRECTIONS.md       every deviation from it, numbered, with evidence
docs/LEDGER.md            what was built and reviewed
docs/RUNBOOK.md           operator guide
tests/                    unit tests for the pure logic
```

## Contamination control

The whole argument that a candidate is real rests on **independent replication**: the same
element recovered from unrelated BioProjects, sequencing centers, and platforms. A contaminant
or an index-hopping artifact does not do that. `results/replication_matrix.tsv` is the
evidence, and candidates are ranked for wet-lab confirmation by independence rather than
abundance — the most-replicated element is the least likely to be an artifact.

## Reproducibility

- `ref/PROVENANCE.tsv` records the source URI, size, SHA-256 and fetch time of every reference
  input.
- `LOGAN_RELEASE=v1.0` pins the pipeline to the December-2023 Logan release so the published
  sweep can be reproduced as a control.
- Published parameters from Zheludev et al. are reproduced exactly and are marked in the code
  as scientific constants.
- `scripts/18_audit_claims.sh` re-checks every automatable factual claim, so the record can be
  refreshed rather than re-derived.

## Before quoting any number

`docs/CORRECTIONS.md` ends with an audited provenance table. Figures not in it — including
several quoted in the original plan, such as the Oblin-1 mean pLDDT and the Foldseek E-value —
are **unverified**, and `OPEN-C` names each one. Do not put an unverified figure on a poster.

## Credits

Method: Zheludev et al., *Cell* 187:6521 (2024) · [VNom](https://github.com/Zheludev/VNom) ·
[circuclust](https://github.com/rcedgar/circuclust).
Data: [Logan](https://github.com/IndexThePlanet/Logan) (Chikhi, Raffestin, Korobeynikov, Edgar,
Babaian), via the AWS Registry of Open Data.
