# CORRECTIONS — deviations from the original PLAN.md

**Purpose.** This file exists so nobody has to re-derive why the code differs from the plan.
Every deviation is numbered, quotes the plan text it overrides, states the evidence, and
names the files that implement it. Revert any single entry by reading it top to bottom.

**Verification date:** 2026-08-19. **Verifier:** automated agent session, branch
`claude/skills-verification-i6l0ue`.
**Original plan:** `docs/PLAN.original.md` (verbatim, unmodified — the corrected working
plan is `PLAN.md` at the repo root).

**Re-verify before ISEF.** Everything below is a snapshot of live network resources.
Re-run `scripts/00_verify_logan.sh` and re-read this file before you present.

---

## Severity key

| | Meaning |
|---|---|
| **CRITICAL** | Silently produces a wrong scientific conclusion. Fixed. |
| **MAJOR** | Wastes substantial compute or forfeits data. Fixed. |
| **MINOR** | Fails loudly, or is cosmetic. Fixed. |
| **OPEN** | Not resolved. Named here so it is not mistaken for settled. |

---

## C-01 — CRITICAL — The Logan freeze date is wrong, and the plan's own "correction" caused it

**Plan text (Verified ground truth table, and the note beneath it):**

> | Logan v1 freeze | **December 2023**, 50 petabases, ~27M accessions | Logan README |
>
> **Correction to an earlier claim:** I previously said the downloadable Logan covers a
> December 2025 freeze at ~87 petabases. The repository states **v1 = December 2023**.
> I could not confirm a newer public release.

**What is actually true.** The retracted claim was correct. The "correction" introduced the
error. From `https://raw.githubusercontent.com/IndexThePlanet/Logan/main/README.md`, fetched
2026-08-19:

> "Logan is a dataset of DNA and RNA sequences. It has been constructed by performing genome
> assembly over a **end-of-2025 freeze** of the entire NCBI Sequence Read Archive, which at
> the time contained **87 petabases** of public raw data."
>
> "## v1.2 release — Nearly all sequencing experiments from a **December 2025 freeze** of the
> SRA have been reconstructed and made available as unitigs and contigs as the **v1.2 release
> of Logan**."
>
> "For older releases notes, see the original release (December 2023 freeze) in [Stats v1]..."

**Root cause.** The plan quoted `Stats-v1.md`, which the README explicitly files under
*"older release notes"*, and treated it as the current release. One page, one channel, no
cross-check against the release index. This is the failure the `verification` discipline
is meant to catch: a single-source read that happened to land on a superseded page.

**Consequence had it shipped.** The pipeline would have been pinned to a deprecated prefix
and would have ignored **10,855,431 accessions** (38,124,741 − 27,269,310) deposited during
2024–2025 — the freshest and least-screened data in the archive, and the part no published
Obelisk survey has ever seen.

**Fix.** `config/config.sh` defaults to `LOGAN_RELEASE=v1.2` and derives the prefix, stats
parquet, size column and accession count from it. `LOGAN_RELEASE=v1.0` pins the old release
for replication. `scripts/00_verify_logan.sh` re-checks the release live and fails loudly
if the README no longer says v1.2.

---

## C-02 — CRITICAL — The incumbent analysis is published inside the bucket the plan streams from

**Plan text (Task 6, Step 1):**

> Search bioRxiv, PubMed, and Google Scholar for: `obelisk rumen`, `obelisk poultry` ...
> **If your niche is claimed, switch niches — not methods.**

The plan directs the incumbent check at the literature only. The decisive artifact is not in
the literature — it is object storage, in the same bucket the pipeline already reads.

**What exists.** `s3://logan-pub/paper/Obelisk/` — the Obelisk authors' own Logan-wide Oblin
sweep, anonymously readable. Artem Babaian is co-lead of Logan *and* senior author of
Zheludev et al. Confirmed present 2026-08-19 (sizes from HTTP `HEAD`):

| Object | Size | Supersedes |
|---|---|---|
| `01diamond_hmm/Obelisk_paper.hmm` | 177 KB | Task 4's `mafft`+`hmmbuild` |
| `01diamond_hmm/Round1_Oblin1_sig.fasta` | 9.7 MB | Task 3 Step 5 (Oblin-1 proteins) |
| `01diamond_hmm/Round1_Oblin2_sig.fasta` | 37 KB | — (Oblin-2; plan ignores Oblin-2 entirely) |
| `02_Building_Obelisk_DB/Obelisk_sig_nt.fasta` | 25.7 MB | Task 3 Step 3 (nucleotide set) |
| `03_Obelisk_DB_QC/Obelisk_db_centroids_nt.fasta` | 7.7 MB | Task 3 Step 3 (centroids) |
| `03_Obelisk_DB_QC/Obelisk_cen_cen_orf.dmnd` | — | Task 4's `diamond makedb` |
| `02_Building_Obelisk_DB/case1/novel_sequences.fasta` | 14.3 MB | — (their novel calls) |
| `03_Obelisk_DB_QC/Circularity_presence.tsv` | — | cross-check for Task 11 |
| `03_Obelisk_DB_QC/Domain_A_presence.tsv` | — | cross-check for Task 11 |

Also present: `01Downloading_+_HMM_Obelisk_Logan.ipynb`, `02_Building_Obelisk_DB.ipynb`,
`03_Obelisk_DB_QC.ipynb` — their full method, executable.

**Two consequences, opposite in sign.**

1. *Good.* Task 3, called "the backbone of the entire project" and requiring a manual
   download of a *Cell* supplementary table, is unnecessary. Better data — already
   QC'd, already clustered, already named in the field's `Obelisk_X_Y_Z` convention — is
   one anonymous `aws s3 cp` away. This dissolves plan Known Unknowns #2 and #3.
2. *Bad, and load-bearing for the project's framing.* "Use Logan to pre-screen the archive
   for Oblins" is **not** an unclaimed idea. The people who built Logan and the people who
   discovered Obelisks are overlapping sets, and they did exactly this.

**What survives as defensible novelty.** Their sweep ran on **v1.0 contigs (accessions up to
December 2023)** — stated explicitly in `Proteins.md` and consistent with the `c1.0` naming
throughout `paper/Obelisk/`. Two things they could not have done:

- **The v1.2 delta.** ~10.9M accessions from 2024–2025, never screened by any published
  Obelisk survey.
- **Niche depth.** A per-niche deep-assembly pass (rnaSPAdes + VNom) on triage survivors,
  rather than an archive-wide HMM pass over pre-computed contigs.

State both plainly in the write-up. Claiming the archive-wide screen as novel would not
survive one informed question at a poster session.

**Fix.** `scripts/01_fetch_references.sh` pulls the published artifacts as the primary path
and keeps the plan's from-scratch build as a documented fallback. `PLAN.md` Task 6 Step 1 now
includes the object-storage channel in the incumbent check.

---

## C-03 — CRITICAL — The decoy control is inverted; `seqkit shuffle` does not shuffle residues

**Plan text (Task 4, Step 1):**

> ```bash
> seqkit shuffle -s 42 "$REF/oblin1_centroids.faa" \
>   | seqkit mutate --any-point 100 -s 42 \
>   > "$REF/decoy_shuffled.faa" || \
>   python - <<'EOF'
> ```

**Two defects, compounding.**

1. `seqkit shuffle` shuffles the **order of records in a file**. It does not permute residues
   within a sequence. Its output is the identical protein set in a different order.
2. `seqkit mutate` has no `--any-point` flag.

**Why this is worse than a crash.** The plan's Step 4 then runs the decoy against the real
database and instructs:

> Expected: **zero or very few lines.** ... **If the decoy hits heavily, your E-value
> threshold is too permissive.** Tighten `DIAMOND_EVALUE` in `config.sh` until the decoy is
> clean ... **This calibration *is* your false-positive control — do not skip it.**

Every "decoy" is a verbatim real Oblin, so every one hits itself at E≈0 — a maximal,
unfixable hit rate. Following the instruction drives `DIAMOND_EVALUE` toward zero until
genuine Oblins stop matching, destroying the sensitivity of the entire triage sweep. The
step that exists to calibrate false positives instead silently guarantees false negatives.

The Python fallback in the plan is correct — it shuffles residues per record — but it is
gated behind `||`, so it runs only if the `seqkit` pipeline *fails*. Because `seqkit shuffle`
succeeds, the correct code never executes.

**Fix.** `scripts/02_build_search_dbs.sh` generates the decoy with the residue-shuffling
Python implementation unconditionally (mono-residue composition preserved, motif order
destroyed, fixed seed 42), then asserts the decoy is clean and **exits non-zero** if it is
not, rather than emitting advice that makes things worse.

---

## C-04 — CRITICAL — `diamond blastx` emits nucleotides; Task 8 feeds them to `hmmsearch`

**Plan text (Task 7 Step 1 writes, Task 8 Step 1 consumes):**

> ```
> --outfmt 6 qseqid sseqid pident length evalue bitscore qseq
> ```
> ```python
> fh.write(f">{r.accession}\n{r.qseq}\n")   # written to triage_best.faa
> subprocess.run(["hmmsearch", ..., f"{ROOT}/ref/oblin1.hmm", faa], check=True)
> ```

In `blastx` mode the query is nucleotide, and DIAMOND's `qseq` field is the aligned segment
of the **query** — i.e. nucleotide sequence. The file is named `.faa` and handed to
`hmmsearch` against an amino-acid profile. HMMER will either reject the alphabet or, worse,
interpret `ACGT` as amino acids and score garbage. The plan calls this step "a real second
filter"; as written it is not a filter at all.

**Fix.** Triage requests `qseq_translated` (DIAMOND ≥ 2.0.8, satisfied by the pinned 2.1.9)
and `scripts/06_collect_hits.py` asserts the collected sequences are amino-acid before
invoking HMMER — failing loudly on any residue outside the protein alphabet.

---

## C-05 — MAJOR — The idempotency guard skips almost nothing

**Plan text (Task 7, Step 1):**

> ```bash
> [ -s "$OUT" ] && { echo "$ACC already done"; exit 0; }   # idempotent
> ...
> else
>   touch "$OUT"; echo "$ACC none"
> ```

`-s` tests "exists **and is non-empty**". No-hit accessions are recorded with `touch`, so
they are empty. On any re-run every no-hit accession fails the guard and is re-downloaded
and re-searched. Since the plan expects hits to be "a small fraction", the guard skips a
small fraction and re-does the overwhelming majority — the opposite of the stated intent,
on the step explicitly designed to survive job failures.

**Fix.** `scripts/05_logan_triage.sh` guards on `[ -e "$OUT" ]` and records completion with
an explicit sentinel so "done, no hits" is distinguishable from "never ran".

---

## C-06 — MAJOR — 100 concurrent SLURM tasks append to one shared file

**Plan text (Task 7 Step 1, with `#SBATCH --array=1-1000%100`):**

> ```bash
> awk -v a="$ACC" '{print a"\t"$0}' "$OUT" >> "$RESULTS/triage_hits.tsv"
> ```

`O_APPEND` writes are atomic only up to `PIPE_BUF` (4096 bytes on Linux) and only per
`write()` call. With `qseq` in the output, hit blocks routinely exceed that, and 100
concurrent writers on a shared filesystem will interleave partial lines. The corrupted
rows then flow into `pd.read_csv(..., names=cols)` in Task 8, which will either raise on
ragged rows or silently mis-assign fields.

**Fix.** Each array task writes only its own `$WORK/triage/$ACC.tsv`. `scripts/06_collect_hits.py`
concatenates them at collection time. No shared-file writes anywhere in the sweep.

---

## C-07 — MAJOR — The cost estimator picks its column by ambiguous substring match

**Plan text (Task 6, Step 5):**

> ```python
> col = [c for c in sub.columns if "contig" in c and "after" in c]
> ...
> if col:
>     tb = sub[col[0]].sum()/1e12
> ```

Verified by reading both parquet footers directly over HTTP range requests on 2026-08-19:

| Parquet | Rows | Columns matching the selector |
|---|---|---|
| `logan-seqstats.parquet` (v1) | 27,269,310 | **2** — `contigs_after_compression`, `size_contigs_after_compression` |
| `logan-seqstats-contigs-v1.2.parquet` (v1.2) | 38,124,741 | **1** — `contigs_after_compression_bytes` |

On v1 the selector matches two columns and `col[0]` takes whichever pandas happens to order
first; only one of them is documented as bytes. The resulting TB figure — the number that
decides whether to spend the compute — is arbitrary. The `if col:` guard also means a schema
change produces **no output at all** rather than an error, on the one step whose entire job
is to stop you overspending.

Both parquets additionally carry nulls for accessions without contigs (v1.2 covers 38.1M
accessions but only 37,377,661 have contigs), which `.sum()` silently drops.

**Fix.** `config/config.sh` exports `LOGAN_SIZE_COL` as an exact name per release.
`scripts/04_select_accessions.sh` raises `KeyError` if the column is absent, and reports
null counts and the matched fraction alongside the total.

---

## C-08 — MAJOR — The contamination control reads only the first VNom FASTA

**Plan text (Task 12, Step 1):**

> ```python
> fas = [f for f in os.listdir(path) if f.endswith(".fasta")]
> if not fas: continue
> out = subprocess.run(["blastn","-query",f"{path}/{fas[0]}", ...
> ```

`fas[0]` — from an unsorted `os.listdir` — discards every other FASTA VNom produced for that
accession. Occurrences are undercounted, so `n_bioprojects` is undercounted, so genuine
candidates fail the `>= MIN_INDEPENDENT_BIOPROJECTS` test at GATE 3. The plan calls this
"the core anti-contamination control" and "the single most defensible validation step";
a nondeterministic undercount here rejects real discoveries silently.

**Fix.** `scripts/11_replication_matrix.py` iterates every `.fasta` under each accession
directory in sorted order and unions the occurrences.

---

## C-09 — MINOR — `hmmpress` is for `hmmscan`, not `hmmsearch`

**Plan text (Task 4, Step 1):** `hmmpress "$REF/oblin1.hmm"`

`hmmsearch` — the only HMMER search the pipeline uses (Tasks 8, 11) — reads the plain `.hmm`.
`hmmpress` builds binary indices for `hmmscan`. Harmless, except it exits non-zero if the
`.h3*` files already exist, which breaks the `set -euo pipefail` script on any re-run.

**Fix.** `hmmpress -f` (idempotent), retained so `hmmscan` remains available, with a comment
explaining it is not required by this pipeline.

---

## C-10 — MINOR — `sed -i` array-bound patch fails silently after a manual edit

**Plan text (Task 9, Step 5):**

> ```bash
> sed -i "s/#SBATCH --array=1-100%20/#SBATCH --array=1-${N}%20/" ...
> ```

The plan's own Step 4 comment tells you to `ADJUST` that line by hand. Once you have, the
literal pattern no longer matches, `sed` reports success, and you submit a 100-task array
over a longer list — silently truncating the deep-assembly stage.

**Fix.** The sbatch scripts read their bound from the accession file at submit time via a
wrapper (`scripts/submit_array.sh`) instead of rewriting themselves.

---

## C-11 — MINOR — `conda activate` in a non-interactive SLURM shell

**Plan text (Tasks 7 and 9 sbatch):** `source activate obelisk 2>/dev/null || conda activate obelisk`

Neither form works in a non-interactive batch shell until conda's shell hook is sourced;
`conda activate` fails with "shell not properly configured".

**Fix.** Both sbatch scripts source `$(conda info --base)/etc/profile.d/conda.sh` first, and
abort with a clear message if conda cannot be located.

---

## C-12 — MINOR — ORFs that cross the circular origin are missed

**Plan text (Task 11, Step 1):** `longest_orf` translates the linear sequence in six frames.

Obelisk genomes are circular; Oblin-1 ORFs frequently span the point where the assembler
linearised the molecule. Searching only the linear form truncates exactly those ORFs, biasing
the `ORF>=100aa` count downward on true positives.

**Fix.** `scripts/09_circularity_check.py` searches the doubled sequence and caps ORF length
at the molecule length, reporting origin-spanning ORFs in a separate column.

---

## C-13 — MINOR — `cut -d, -f22` assumes a fixed runinfo column order

**Plan text (Task 6, Step 3):** `cut -d, -f22 "$WORK/niche_runinfo.csv"   # BioProject column`

SRA `runinfo` column order is not contractual. Position 22 silently yields the wrong field
if NCBI reorders.

**Fix.** All runinfo parsing selects by header name.

---

## C-14 — MAJOR — Every circuclust flag in the plan is wrong

**Plan text (Task 14, Step 5):**

> ```bash
> "$CIRCUCLUST/bin/circuclust" --input "$WORK/all_obelisks.fna" --id 0.80 \
>   --output "$RESULTS/clusters80.tsv"
> ```
> (Confirm exact flags against circuclust's README.)

The plan flagged this as unverified (its Known Unknown #5). It is now verified, and all
three flags are wrong. From circuclust's README, fetched 2026-08-19 from
`raw.githubusercontent.com/rcedgar/circuclust/master/README.md`:

> ```
> circuclust -cluster seqs.fa -id 0.9 -fastaout centroids.fa -tsvout hits.tsv
> ```

Single dashes, not double. The input file is the **value of `-cluster`**, not of `--input`.
Output is split between `-fastaout` (centroids) and `-tsvout` (assignments); there is no
`--output` at all. The plan's command would fail immediately.

**Fix.** `scripts/16_cluster.sh` uses the documented form and runs both the 80% and 95%
levels, which the `Obelisk_X_Y_Z` convention requires.

---

## C-15 — CRITICAL — The plan's VNom invocation cannot work

**Plan text (Task 10 Step 1, and identically in Task 5 Step 6):**

> ```bash
> python "$VNOM/VNom.py" -i "$T" \
>   -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1 \
>   -o "$WORK/vnom/$ACC"
> ```
> where `$T` is `$WORK/deep/<ACC>/transcripts.fasta`.

The plan flagged the invocation as unverified (its Known Unknown #4). It is now verified
against VNom's README, fetched 2026-08-19, whose own worked example is:

> ```
> python ../VNom.py -i peach_subset -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1
> ```

Four separate incompatibilities, any one of which breaks the run:

1. **`-i` takes a basename, not a path, and without the extension.** The README is explicit:
   *"you must specify this single underscore name without the file ending for VNom"*. The
   plan passes a full path ending in `.fasta`.
2. **There is no `-o` flag.** Output goes to a `4_final_clusters` directory relative to the
   working directory: *"outputs are stored to `4_final_clusters` (so if this dir wasn't
   written, VNom failed to nominate viroid-like contigs)"*. Isolating accessions therefore
   requires running each in its own directory.
3. **The filename must contain exactly one underscore** and end in `.fasta` — *"X_Y.fasta is
   good, but XY.fasta is bad"*. `transcripts.fasta` contains none.
4. **seqIDs must keep the default rnaSPAdes layout.** *"adding more underscores will cause
   VNom to crash"*. The README's example substitutes the accession for the literal `NODE`.

Two further constraints the plan never mentions, both from the same README:

- **Input must come from stranded RNA-seq.** VNom's third filter keeps only clusters
  containing both polarities, which an unstranded library cannot satisfy. The author notes
  this is where VNom most often stops.
- **VNom depends on circuclust, USEARCH and MARS** installed under `VNom/dependencies/`.
  The plan's Task 2 installs none of them.

**Consequence.** Since the plan's downstream globs look for `*.fasta` directly under the
accession directory, and VNom writes to `4_final_clusters`, even a VNom run that succeeded
would have produced zero pooled candidates — a silent empty result at Task 10, propagating
to an empty GATE 3 and a false negative finding.

**Fix.** `scripts/08_run_vnom.sh` renames input to `<ACC>_contigs.fasta`, substitutes the
accession for `NODE`, strips N-containing contigs as the README's example does, runs VNom
from inside a per-accession directory with the bare stem, and treats a missing
`4_final_clusters` as the documented "nothing nominated" outcome rather than an error.
`scripts/03_positive_control.sh` and `scripts/11_replication_matrix.py` read from
`4_final_clusters`.

---

## OPEN-A — NCBI E-utilities are unreachable from the build environment

`esearch`/`efetch` (Task 5 Step 1, Task 6 Step 3) could not be exercised here:
`eutils.ncbi.nlm.nih.gov:443` is refused by this container's egress policy
(`connect_rejected`, gateway 403). The scripts are written as specified and marked
cluster-only; they are **unvalidated against a live NCBI response**.

**Resolved for niche selection.** `scripts/04b_select_accessions_offline.sh` enumerates
accessions from Logan's own `stats/logan_accessions_v1.2_SRA2025.csv.zst` (11.0 GiB) joined
to `stats/sra_taxid.csv.zst` (424 MiB). Decompressing the headers on 2026-08-19 confirmed
the accession table carries `acc, assay_type, center_name, librarysource, organism,
platform, bioproject` — including **BioProject and CenterName**. That was not a given, and
it matters: GATE 3 rests entirely on BioProject independence, so an offline path lacking it
could have found candidates but never validated them.

Still NCBI-only: the SK36 run enumeration in `scripts/03_positive_control.sh`, which accepts
`--runs FILE` as a manual substitute. `docs/RUNBOOK.md` says what to run where.

**A caveat on niche size.** Filtering a 120 MB slice of the accession table found 11 rumen
keyword hits, all AMPLICON or WGS against METAGENOMIC/GENOMIC source — i.e. DNA — and none
passing the RNA-Seq filter. Do **not** extrapolate a niche size from that: the table is
sorted by `center_name`, so a prefix slice is a biased sample of sequencing centers rather
than a random sample of the archive. A full pass is required before concluding whether the
rumen metatranscriptome clears the plan's own 200-run floor. The vocabulary itself is
confirmed correct: `RNA-Seq` appears 75,366 times and `METATRANSCRIPTOMIC` 2,101 times in
that slice.

---

## OPEN-B — The skills-registry question from plan Known Unknown #7 is still open

**Plan text:**

> `npx skills find` returned nothing for six bioinformatics queries **and also nothing for the
> control query `react`** ... **no conclusion can be drawn about what bioinformatics skills
> exist.** Run `npx skills find bioinformatics` on your own machine.

Re-attempted 2026-08-19 in this environment, across two channels:

| Channel | Result |
|---|---|
| `npx -y skills find bioinformatics` | `No skills found` |
| `npx -y skills find react` (control) | `No skills found` |
| `curl https://skills.sh` | exit 56, HTTP `000` — connection refused |
| `curl https://api.skills.sh` | exit 56, HTTP `000` — connection refused |

The control query fails identically here, so this reproduces the plan's non-result rather
than resolving it. `skills.sh` is not reachable from this container's egress policy.
**No conclusion can be drawn about what bioinformatics skills exist.** Still to be run on a
machine with unrestricted network access.

---

## OPEN-C — Claims that could not be checked from here

Stated so they are not mistaken for verified. These are the plan's numbers, carried forward
unaudited.

| Claim | Plan's source | Status |
|---|---|---|
| 1,744 stringent Obelisks; 788 with ≥2 ORFs | Zheludev et al. Table S1 | Not checked — table is behind the *Cell* site; `paper/Obelisk/` artifacts may supersede it. Verify with `seqkit stats` after `scripts/01`. |
| ~1,700 Oblin-1 centroid proteins | Zheludev et al. | Not checked; `Round1_Oblin1_sig.fasta` is 9.7 MB, consistent in magnitude but not counted. |
| 40 marine Obelisks | López-Simón et al., *ISME J* 19:wraf033 | Not checked — publisher not reachable from here. |
| Mean pLDDT 83.8 for Oblin-1 | Zheludev et al. | Not checked. Do not quote on a poster until you have read it in the paper. |
| Foldseek best E-value 0.31 | Urayama et al., *Nat Commun* 17:3041 | Not checked. |
| 17 SK36 runs, four labs | *J Mol Evol* 93:370 (2025) | Not checked. |
| VNom CLI entry point and flags | VNom README | Not checked — `api.github.com` returns 403 here. `raw.githubusercontent.com` is reachable, so `scripts/01` fetches the README directly for you to read. |
| circuclust flags | circuclust README | Same as above. |

`scripts/18_audit_claims.sh` re-runs every check in this table that can be automated, so
this section can be refreshed rather than re-derived.

---

## Numbers this project may quote, with provenance

Audited 2026-08-19. Anything not in this table should not go on a poster without a source.

| Value | Figure | Source, verified |
|---|---|---|
| Logan v1.2 release date | 21 Apr 2026 | `Stats-v1.2.md` |
| v1.2 SRA cutoff | 31 Dec 2025 | `Stats-v1.2.md` |
| v1.2 accessions, unitigs | 38,100,000 (README) / **38,124,741** (parquet rows) | `Stats-v1.2.md`; parquet footer |
| v1.2 accessions, contigs | 37,377,661 | `Stats-v1.2.md` |
| v1.2 raw bases assembled | 86.6 Pbp | `Stats-v1.2.md` |
| v1.2 contigs compressed | 623 TB | `Stats-v1.2.md` |
| v1.2 unitigs compressed | 4.19 PB | `Stats-v1.2.md` |
| v1.0 cutoff | 10 Dec 2023 | `Stats-v1.md` |
| v1.0 accessions, unitigs | 27,300,000 (README) / **27,269,310** (parquet rows) | `Stats-v1.md`; parquet footer |
| v1.0 accessions, contigs | 26,788,835 | `Stats-v1.1.md` |
| v1.0 raw bases | 48.2 Pbp | `Stats-v1.md` |
| **Unscreened delta, v1.0 → v1.2** | **10,855,431 accessions** | computed: 38,124,741 − 27,269,310 |
| SRA at Dec 2025 | 87 Pbp, 38M accessions | Logan README |
| Assembler | Minia3, k=31 | `Contigs.md` |
| Auth | `--no-sign-request` | `Contigs.md`, `Accessions.md` |

**Two caveats on the above.**

1. **`Contigs.md` is stale.** It still reports "315 terabytes" and "26.7 million accessions"
   — v1.1-era figures — while `Stats-v1.2.md` reports 623 TB and 37.4M for the same `c/`
   prefix. The plan cites `Contigs.md`, so it quoted the page accurately but the page is
   superseded. Prefer the release-notes figures.
2. **README round numbers differ from parquet row counts** (38.1M vs 38,124,741; 27.3M vs
   27,269,310). Both are correct at their stated precision. Quote the exact figure only when
   citing the parquet, and never present a rounded README figure as exact.

The v1-era library breakdown the plan cites — metagenomes 4,791,129 accessions / 43.4 TB —
is **v1.0 only**; no v1.2 equivalent is published. Note also that the same page reports
**metatranscriptomes: 129,974 accessions, 2.6 TB**, which is the library type this project
actually targets, and which the plan never cites.
