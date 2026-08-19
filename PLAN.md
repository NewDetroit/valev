# Computational Pipeline Plan — Obelisk-like Element Discovery in an Under-Mined Niche

**Project:** ISEF 2027 (Los Angeles, May 8–14, 2027)
**Plan written:** 2026-08-17 · **Ground truth re-verified:** 2026-08-19
**Original, unmodified:** [`docs/PLAN.original.md`](docs/PLAN.original.md)
**Every deviation from it:** [`docs/CORRECTIONS.md`](docs/CORRECTIONS.md)

> This is the **working** plan. It differs from the original in 15 numbered places, each of
> which is a defect that would have produced a wrong result or wasted compute. Where a step
> below carries a `C-nn` tag, read that correction before changing the step back.

---

## Scope

One subsystem: the computational discovery and validation pipeline, ending at a ranked,
validated candidate set packaged for wet-lab confirmation. Excludes wet-lab protocols, ISEF
paperwork, and GenBank submission.

Two architectures joined at a triage point:

- **Fast path (Logan):** screen millions of pre-assembled SRA accessions cheaply to find which
  plausibly contain an Oblin-like protein.
- **Deep path (Zheludev replication):** run the full published assembly pipeline *only* on
  accessions that survive triage.

Zheludev et al. ran the deep path on a curated sample set. Logan lets you pre-filter the whole
archive first, so the expensive step runs a few hundred times instead of a few million.

---

## Verified ground truth

Checked against primary sources 2026-08-19. Do not modify without re-verifying; re-check with
`scripts/00_verify_logan.sh`, which fails loudly if Logan's release string changes.

| Fact | Value | Source |
|---|---|---|
| Logan S3 bucket | `s3://logan-pub/` | AWS Open Data Registry, `pasteur-logan` |
| Top-level prefixes | `c/` `c1.0/` `p/` `paper/` `stats/` `u/` | live bucket listing |
| Contig path | `s3://logan-pub/c/{ACC}/{ACC}.contigs.fa.zst` | `Contigs.md` |
| Unitig path | `s3://logan-pub/u/{ACC}/{ACC}.unitigs.fa.zst` | `Unitigs.md` |
| Auth | `--no-sign-request` (no AWS account) | `Accessions.md` |
| **Current release** | **v1.2, published 21 Apr 2026** | `Stats-v1.2.md` |
| **SRA cutoff** | **31 December 2025** | `Stats-v1.2.md` |
| **Accessions (unitigs)** | **38,124,741** (parquet rows) | parquet footer |
| Accessions (contigs) | 37,377,661 | `Stats-v1.2.md` |
| Raw bases assembled | 86.6 petabases | `Stats-v1.2.md` |
| Contigs compressed | 623 TB | `Stats-v1.2.md` |
| Unitigs compressed | 4.19 PB | `Stats-v1.2.md` |
| Superseded release | v1.0, Dec 2023 cutoff, 27,269,310 accessions, at `c1.0/` | `Stats-v1.md`, `Stats-v1.1.md` |
| Per-accession stats | `stats/logan-seqstats-contigs-v1.2.parquet` | `Stats-v1.2.md` |
| Size column | `contigs_after_compression_bytes` | parquet footer |
| Assembler | Minia3, k=31 | `Contigs.md` |

**C-01 — the original plan had this inverted.** It asserted "v1 = December 2023" and
*retracted a correct claim* of a December-2025 / 87-petabase release. It had read
`Stats-v1.md`, which the README files under "older release notes". Pinning to v1.0 would have
discarded **10,855,431 accessions** of 2024–2025 data — the least-screened part of the archive.

**C-02 — the incumbent is in the bucket, not the literature.** `s3://logan-pub/paper/Obelisk/`
holds the Obelisk authors' own Logan-wide Oblin sweep: `Obelisk_paper.hmm`,
`Round1_Oblin1_sig.fasta`, `Round1_Oblin2_sig.fasta`, `Obelisk_cen_cen_orf.dmnd`,
`novel_sequences.fasta`, `Circularity_presence.tsv`, `Domain_A_presence.tsv`, and three
executable notebooks. Babaian is co-lead of Logan *and* senior author of Zheludev et al.
Consequences both ways: the reference set is a download rather than a rebuild, **and** an
archive-wide Oblin screen is not an unclaimed idea. Their sweep ran on v1.0 — the v1.2 delta
and the per-niche deep pass are what remain novel. Say so in the write-up.

**Published Obelisk pipeline** (Zheludev et al., *Cell* 187:6521, 2024) — scientific constants,
reproduced exactly:

```
fasterq-dump → fastp (--average_qual=30 --n_base_limit=0 --cut_front --cut_tail)
             → rnaSPAdes (default settings)
             → VNom (-max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1)
```

**Positive control:** Obelisk-S.s in *Streptococcus sanguinis* SK36 RNA-seq — multiple
independent datasets, known answer.

---

## File structure

```
config/config.sh              single source of truth; LOGAN_RELEASE selects v1.2 or v1.0
config/niche_query.txt        Entrez query defining the target niche
env/environment.yml           conda spec (diamond=2.1.9 pinned — see C-04)
ref/                          built once by 01–02, never edited; PROVENANCE.tsv traces every file
scripts/00_verify_logan.sh              Phase 0
scripts/01_fetch_references.sh          Phase 1
scripts/02_build_search_dbs.sh          Phase 1
scripts/03_positive_control.sh          Phase 2 — GATE 1
scripts/04_select_accessions.sh         Phase 3 (needs NCBI)
scripts/04b_select_accessions_offline.sh  Phase 3 (no NCBI — see OPEN-A)
scripts/estimate_sweep_cost.py          Phase 3 — see C-07
scripts/05_logan_triage.sh / .sbatch    Phase 3
scripts/submit_array.sh                 SLURM bound computed at submit — see C-10
scripts/06_collect_hits.py              Phase 3 — GATE 2
scripts/07_deep_assemble.sh / .sbatch   Phase 4
scripts/08_run_vnom.sh                  Phase 4
scripts/09_circularity_check.py         Phase 5 — see C-12
scripts/10_exclude_known.sh             Phase 5
scripts/11_replication_matrix.py        Phase 5 — GATE 3, see C-08
scripts/12_junction_mapping.sh          Phase 5 — the strongest evidence in the project
scripts/13_rnafold.sh                   Phase 6
scripts/14_colabfold_prep.py            Phase 6
scripts/15_foldseek.sh                  Phase 6
scripts/16_cluster.sh                   Phase 6 — circuclust nomenclature
scripts/17_wetlab_package.py            Phase 7
scripts/18_audit_claims.sh              refreshes the OPEN-C claim table
work/ logs/                   disposable
results/                      the only directory that matters for the poster
tests/                        unit tests for the pure logic
```

`config/config.sh` is the single source of paths and thresholds — **no script hardcodes them.**
`ref/` is built once and never modified. `work/` is disposable. `results/` is the deliverable.

**Exit codes:** `0` success · `1` unexpected error · `2` missing input, run an earlier script ·
`3` a GATE failed — a scientific result, not a crash.

---

# PHASE 0 — Verification and environment

### Task 1 — Verify Logan is live and the ground truth still holds

`scripts/00_verify_logan.sh`

- [ ] Confirms the bucket is reachable anonymously and **reports every top-level prefix**.
- [ ] Fetches Logan's README and **fails (exit 3) if it no longer advertises v1.2** — this is
      the automated guard against C-01 recurring.
- [ ] Smoke-tests a known-present accession under `$LOGAN_C`.
- [ ] Confirms `$LOGAN_STATS_PARQUET` exists and reports its size. Does **not** download it by
      default; it is ~1.7 GB and needs an explicit `--download`.
- [ ] Confirms every `$LOGAN_OBELISK` artifact C-02 depends on.

**Decision point.** If a release newer than v1.2 exists, take it and update `config.sh` —
`LOGAN_RELEASE` is the only place that changes. If the freeze date limits your niche, state
that limitation explicitly in your paper. Disclosing it is a credibility gain, not a loss.

### Task 2 — Build the environment

`env/environment.yml`, `docs/RUNBOOK.md`

- [ ] `conda env create -f env/environment.yml && conda activate obelisk`
- [ ] Verify every tool resolves. **Do not proceed past any `MISSING`.**
- [ ] `diamond=2.1.9` is pinned deliberately: `qseq_translated` needs ≥2.0.8 (C-04).
- [ ] Clone VNom and circuclust; record the exact invocation each needs — neither CLI could be
      verified from the build environment (`OPEN-C`).

---

# PHASE 1 — Reference construction

### Task 3 — Obtain the Obelisk reference set

`scripts/01_fetch_references.sh`

**C-02: the primary path is a download from `$LOGAN_OBELISK`, not a journal supplement.** The
original plan called Table S1 "the backbone of the entire project" and required fetching it by
hand from *Cell*. Better data — already QC'd, already clustered, already named in the field's
`Obelisk_X_Y_Z` convention — is one anonymous `aws s3 cp` away. This dissolves the original
plan's Known Unknowns #2 and #3.

- [ ] Fetch the published HMM, Oblin-1 and Oblin-2 protein sets, nucleotide centroids, and the
      authors' circularity and domain tables.
- [ ] Write `ref/PROVENANCE.tsv` — source URI, bytes, SHA-256, fetch time for every file.
      This is what lets the write-up cite its inputs.
- [ ] `--from-scratch` retains the original mafft/hmmbuild path as a documented fallback.
- [ ] Sanity-check: `seqkit stats ref/obelisk_nt.fna` should show mean length near 1,000 nt.
      **If it is far off, the wrong object was selected.**

### Task 4 — Build search databases

`scripts/02_build_search_dbs.sh`

- [ ] DIAMOND protein database for translated triage.
- [ ] Profile HMM for sensitive confirmation — the published one, which also makes results
      directly comparable to the literature.
- [ ] **C-03 — the decoy control.** The original did
      `seqkit shuffle | seqkit mutate --any-point 100`. `seqkit shuffle` shuffles the *order of
      records*, not residues, and `--any-point` is not a real flag — so every "decoy" was a
      verbatim real Oblin hitting itself at E≈0. The original then instructed the operator to
      tighten `DIAMOND_EVALUE` until the decoy came back clean, which would have driven the
      threshold down until genuine Oblins stopped matching. **The step meant to control false
      positives instead guaranteed false negatives.** Now: residues are shuffled per record
      (seed 42, composition preserved), and a dirty decoy **exits 3** instead of emitting advice.

This calibration *is* the false-positive control. Do not skip it.

---

# PHASE 2 — Positive control (GATE 1)

> **The most important phase.** If the pipeline cannot re-find a known Obelisk where one is
> known to exist, every downstream result is meaningless — including a negative one.

### Task 5 — Recover Obelisk-S.s from *S. sanguinis* SK36

`scripts/03_positive_control.sh` — stages: `runinfo | triage | deep | confirm | all`

- [ ] Identify SK36 RNA-seq runs; select five spanning **≥2 BioProjects**, enforced as a real
      check. This deliberately exercises the same independence logic the contamination control
      relies on later.
- [ ] Run the Logan fast path on them.
- [ ] **GATE 1.** Hits → proceed. Zero hits → debug in order: (a) is the Oblin FASTA real
      protein, (b) do the S3 objects exist, (c) retry at `$DIAMOND_EVALUE_RELAXED`, (d) the runs
      may postdate the freeze — though under C-01 the freeze is now 2025-12-31, making (d) far
      less likely than the original plan assumed.
- [ ] Run the deep path on the best-scoring accession, with Zheludev's exact parameters.
- [ ] Confirm >95% identity to the known Obelisk-S.s; write
      `results/positive_control_evidence.tsv`.

**Do not continue past this gate on a failure.** A silent zero here becomes a false "we found
nothing novel" later. This output is a poster figure and the most persuasive thing you can show
a skeptical judge.

---

# PHASE 3 — Niche selection and triage

### Task 6 — Define and enumerate the target niche

`config/niche_query.txt`, `scripts/04_select_accessions.sh`,
`scripts/04b_select_accessions_offline.sh`, `scripts/estimate_sweep_cost.py`

- [ ] **Re-run the incumbent check before committing.** Search the literature *and* — per C-02
      — object storage. Record date and exact queries in your lab notebook. Repeat monthly.
      **If your niche is claimed, switch niches, not methods.** Everything downstream is reusable.
- [ ] Write the Entrez query. **RNA-Seq strategy is mandatory** — Obelisks are RNA elements.
- [ ] Enumerate accessions. Runinfo is parsed **by header name, never column position** (C-13).
      Sanity thresholds are enforced checks: under ~200 runs the niche is too small, over
      ~200,000 it is too large.
- [ ] **OPEN-A:** if NCBI is unreachable, `04b` enumerates from Logan's own
      `stats/logan_accessions_v1.2_SRA2025.csv.zst` instead.
- [ ] **C-07 — estimate cost before spending.** The original selected the size column by
      substring match, which matches *two* columns on the v1 parquet and takes `col[0]`
      arbitrarily, and produces *no output at all* if the schema shifts — on the one step whose
      job is to stop you overspending. Now the exact column comes from `$LOGAN_SIZE_COL` and a
      missing column raises. Null counts and matched fraction are reported alongside the total.

Read two things from the estimate: the match fraction tells you how much of your niche predates
the freeze, and the TB figure sizes the sweep — transfer volume, not disk, because you stream.

### Task 7 — Run the Logan triage sweep

`scripts/05_logan_triage.sh`, `scripts/05_logan_triage.sbatch`, `scripts/submit_array.sh`

- [ ] Stream contigs straight from S3 through `zstdcat` into DIAMOND. Nothing touches disk.
- [ ] **C-04** — request `qseq_translated`, not `qseq`. In blastx, `qseq` is nucleotide, and the
      original fed it to `hmmsearch` against an amino-acid profile.
- [ ] **C-05** — the idempotency guard tests the `.done` sentinel, not `-s`. The original's
      `-s` test meant every no-hit accession — the overwhelming majority — was re-downloaded
      and re-searched on any re-run.
- [ ] **C-06** — each task writes only its own file. The original had 100 concurrent array
      tasks appending to one shared TSV; with `qseq` in the payload, lines exceed `PIPE_BUF`
      and interleave.
- [ ] **C-10/C-11** — `submit_array.sh` computes the array bound at submit time instead of
      `sed`-patching the sbatch file, and the sbatch sources conda's shell hook before activating.
- [ ] Missing accessions are normal. The worker distinguishes *absent from Logan*, *download
      failed*, and *no hits*, so the coverage statistic is honest.
- [ ] **Test on ten accessions and time it** before scaling. Then a 5-task pilot. **Never
      submit the full array before the pilot is clean.**

### Task 8 — Confirm triage hits with the profile HMM (GATE 2)

`scripts/06_collect_hits.py`

- [ ] Collect per-accession outputs, take the best hit each, confirm with HMMER — a filter
      independent of DIAMOND, so it is a genuine second opinion.
- [ ] Assert the collected sequences are amino-acid before HMMER runs (C-04).
- [ ] **If the hit rate exceeds ~20%, stop (exit 3).** That is far above biological
      plausibility and means the threshold is loose or the C-03 decoy calibration failed.
- [ ] **GATE 2.** ≥3 independent BioProjects → Phase 4. 1–2 → contamination warning, broaden
      and re-sweep. Zero → either a real negative or too strict; retry at `$HMM_EVALUE_RELAXED`
      before concluding.

---

# PHASE 4 — Deep assembly on survivors

### Task 9 — Full Zheludev pipeline on confirmed accessions

`scripts/07_deep_assemble.sh`, `scripts/07_deep_assemble.sbatch`

- [ ] The deep list should be far shorter than the triage list — tens to a few hundred.
      **That reduction is the entire point of the hybrid architecture.**
- [ ] Exact published `fastp` parameters. Handles paired, single, and orphan reads.
- [ ] Raw FASTQs are deleted only **after** a confirmed successful assembly, so a retry does not
      re-download everything.
- [ ] **Time and measure peak memory on one accession first.** rnaSPAdes on a large
      metatranscriptome can exceed 64 GB; the `-m` value and `--mem` are derived from one source
      so they cannot drift apart.

### Task 10 — Run VNom to call circular candidates

`scripts/08_run_vnom.sh`

- [ ] Exact published VNom flags.
- [ ] Pool, filter to `$OBELISK_MIN_LEN`–`$OBELISK_MAX_LEN`, dedup with `seqkit rmdup -s` so one
      element recovered from twenty libraries counts once — while **retaining the
      candidate→accession map**, which a bare dedup destroys and which Task 12 requires.

---

# PHASE 5 — Validation (GATE 3)

> Everything before this produces *candidates*. This produces *evidence*. Judges and reviewers
> will attack precisely here.

### Task 11 — Circularity and ORF confirmation

`scripts/09_circularity_check.py`, `scripts/10_exclude_known.sh`

- [ ] Terminal-repeat detection reports the **length** of the maximal repeat, not a boolean, and
      rejects low-complexity repeats that would otherwise false-positive on ~1 kb sequences.
      This is a heuristic — true circularity is settled in Task 12.
- [ ] **C-12** — ORFs are found on the **doubled** sequence. Obelisk genomes are circular and
      Oblin-1 ORFs frequently span the point where the assembler linearised the molecule; a
      linear-only search truncates exactly the ORFs you most want. Origin-spanning ORFs are
      reported in their own column.
- [ ] Confirm ORFs are Oblin-like against the profile HMM.
- [ ] **Exclude known Obelisks — this defines what is actually new.** Re-discovering known
      Obelisks is a good sign: it independently validates the pipeline. Report both numbers;
      the ratio is itself a credibility statistic.

### Task 12 — The replication matrix (GATE 3)

`scripts/11_replication_matrix.py`, `scripts/12_junction_mapping.sh`

The core anti-contamination control: which candidates recur across unrelated BioProjects,
sequencing centers, and platforms.

- [ ] **C-08** — every VNom FASTA per accession is read, in sorted order. The original used
      `fastas[0]` from an unsorted `listdir`, silently undercounting occurrences and therefore
      failing genuine candidates at the `≥3 BioProjects` test.
- [ ] Emit `results/replication_matrix.tsv` and `results/validated.fna`. **This file is the
      project.**
- [ ] **Confirm circularity by junction read-mapping.** Build a doubled reference and count
      reads whose alignment **spans position L**. The original printed mean depth over the
      doubled reference, which does not demonstrate circularity at all. This is the single most
      defensible validation step in the plan — make it a poster figure.
- [ ] Screen against contaminants. Anything hitting a cloning vector or PhiX is a contaminant:
      remove it and **say so in your paper.** Documenting removals is a credibility gain.
- [ ] **GATE 3.** ≥1 candidate replicated, circular by junction mapping, HMM-positive, no
      contaminant hits → Phase 6, you have a discovery. Only 1–2 BioProjects → report as
      *tentative*, explicitly flagged; do not overclaim. Nothing passes → a **publishable
      negative result**: *"a systematic search of N accessions across M BioProjects in [niche],
      using a pipeline validated on a positive control, found no novel Obelisk-like elements."*
      The GATE 1 result is what makes that statement meaningful rather than an admission of
      failure. This is a legitimate ISEF project — say so plainly.

---

# PHASE 6 — Structural characterization

### Task 13 — RNA secondary structure

`scripts/13_rnafold.sh`

- [ ] Fold every validated candidate; generate plots for the top candidate.
- [ ] **Look for the rod-like fold** — part of the Obelisk definition, and a striking visual.
      Rod-likeness is computed as a number, not eyeballed.
- [ ] Search Rfam for ribozyme motifs. Some Obelisks carry hammerheads; absence is not
      disqualifying.

### Task 14 — Protein structure, and the phylogeny

`scripts/14_colabfold_prep.py`, `scripts/15_foldseek.sh`, `scripts/16_cluster.sh`

- [ ] Predict structures with ColabFold. **Record pLDDT for every model**; treat <70 as low
      confidence and say so. (The original quotes a mean pLDDT of 83.8 for Oblin-1 — that figure
      is **unverified**, see `OPEN-C`. Do not put it on a poster until you have read it in the
      paper.)
- [ ] Foldseek against AFDB and PDB. **Both outcomes are good:** no significant hits argues a
      genuinely novel fold; hits to known Oblin-1 confirm a true Obelisk relative.
- [ ] Build the phylogeny. **A well-supported clade separate from all published Oblins is the
      phylogenetic version of "this is new."**
- [ ] Cluster with circuclust for nomenclature. The field's convention is `Obelisk_X_Y_Z`
      (X = 80% nt cluster, Y = 95%, Z = strain). Following it makes results directly comparable
      and signals you know the literature — and `$LOGAN_OBELISK` contains real examples to
      validate your naming against.

---

# PHASE 7 — Package for wet-lab handoff

### Task 15 — Build the confirmation request

`scripts/17_wetlab_package.py`, `scripts/18_audit_claims.sh`

- [ ] Rank by **independence**, not abundance: `n_bioprojects*3 + n_centers*2 + n_platforms`.
      The most-replicated element is the least likely to be an artifact.
- [ ] Design **outward-facing (divergent) primers** — they amplify only across the circular
      junction and give no product from a linear template. That is the standard RT-PCR test for
      circularity. Target 100–300 bp.
- [ ] Identify obtainable sample sources from each candidate's BioProjects.
- [ ] Write `results/WETLAB_BRIEF.md`: what an Obelisk is in three sentences, what you found and
      the evidence tier for each candidate, the exact RT-PCR requested, the primers, the sample
      type, and what a positive result would establish.

**Bring this, not an idea.** A one-page request with validated targets and designed primers is a
far easier yes than a conversation about a possibility.

---

## Known unknowns

Flagged rather than papered over. Full detail in `docs/CORRECTIONS.md`.

1. **`OPEN-A`** — NCBI E-utilities were unreachable from the build environment, so
   `04_select_accessions.sh` and the runinfo stages of `03_positive_control.sh` are written to
   spec but **unvalidated against a live NCBI response**. `04b` is the NCBI-free alternative.
2. **`OPEN-B`** — the skills-registry question is still open. Re-attempted across two channels;
   the control query fails identically, so this reproduces the original non-result rather than
   resolving it. `skills.sh` is not reachable from the build environment.
3. **`OPEN-C`** — several literature figures the original plan quotes could not be checked from
   the build environment, including the 1,744 stringent Obelisks, the 40 marine Obelisks, the
   Oblin-1 mean pLDDT of 83.8, and the Foldseek E-value of 0.31. They are carried forward
   **unaudited and marked as such**. `scripts/18_audit_claims.sh` refreshes what can be
   automated.
4. **VNom and circuclust CLIs** — the *parameters* come from published methods and are reliable;
   the *invocation* is isolated behind one variable per tool so there is a single place to fix.
5. **Cluster specifics** — partition names, account strings, and module systems vary; marked
   `ADJUST`.

## What changed, and why it matters

The original plan's self-review concluded "No gaps." It was wrong in thirteen places, and the
two that mattered most were not coding errors but **verification** failures: a single-source
read of a superseded documentation page (C-01), and an incumbent check aimed only at the
literature when the decisive artifact was in object storage (C-02). Both went undetected
because the plan checked *more*, not *wider*.

The pipeline is stronger for it. It now searches 38.1M accessions instead of 27.3M, starts from
the field's own reference data instead of a hand-downloaded supplement, and has a false-positive
control that actually controls false positives.
