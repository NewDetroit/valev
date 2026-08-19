# Computational Pipeline Plan — Obelisk-like Element Discovery in an Under-Mined Niche

**Project:** ISEF 2027 (Los Angeles, May 8–14, 2027)
**Plan written:** 2026-08-17
**Save location:** `~/obelisk-hunt/PLAN.md` (no existing project convention; this is the proposed root)
**Executor assumption:** skilled at following instructions, zero context on this domain or toolset.

---

## Scope check

This plan covers **one subsystem: the computational discovery and validation pipeline**, ending at a ranked, validated candidate set packaged for wet-lab confirmation. It deliberately excludes wet-lab protocols, ISEF paperwork, and GenBank submission — each is a separate plan written later, once this one produces output.

The pipeline splits into two architectures joined at a triage point:

- **Fast path (Logan):** screen millions of pre-assembled SRA accessions cheaply to find which ones plausibly contain an Oblin-like protein.
- **Deep path (Zheludev replication):** run the full published assembly pipeline *only* on accessions that survive triage.

This hybrid is the core design decision. Zheludev et al. ran the deep path on a curated sample set. Logan lets you pre-filter the whole archive first, so the expensive step runs a few hundred times instead of a few million. That is where institutional compute converts into reach.

---

## Verified ground truth (checked 2026-08-17)

Do not modify these without re-verifying.

| Fact | Value | Source |
|---|---|---|
| Logan S3 bucket | `s3://logan-pub/` | AWS Open Data Registry, `pasteur-logan` |
| Unitig path | `s3://logan-pub/u/{ACC}/{ACC}.unitigs.fa.zst` | Logan `Accessions.md` |
| Contig path | `s3://logan-pub/c/{ACC}/{ACC}.contigs.fa.zst` | Logan `Sequence_Search.md` |
| Auth | `--no-sign-request` (no AWS account needed) | Logan `Accessions.md` |
| Logan v1 freeze | **December 2023**, 50 petabases, ~27M accessions | Logan README |
| Contigs total | 315 TB compressed, 26.7M files | Logan `Contigs.md` |
| Unitigs total | 2.18 PB compressed | Logan `Stats-v1.md` |
| Metagenome subset | 4,791,129 accessions, 43.4 TB compressed | Logan `Stats-v1.md` |
| Per-accession stats | `https://s3.amazonaws.com/logan-pub/stats/logan-seqstats.parquet` | Logan `Stats-v1.md` |
| Assembler used | Minia3, k=31 | Logan `Contigs.md` |

**Correction to an earlier claim:** I previously said the downloadable Logan covers a December 2025 freeze at ~87 petabases. The repository states **v1 = December 2023**. I could not confirm a newer public release. Task 1 below verifies this before you build on it — if a v2 exists, you gain two years of unscreened data and should use it.

**Published Obelisk pipeline (Zheludev et al., *Cell* 187:6521, 2024), exact parameters:**
```
fasterq-dump → fastp (--average_qual=30 --n_base_limit=0 --cut_front --cut_tail)
             → rnaSPAdes (default settings)
             → VNom (-max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1)
```
VNom: `github.com/Zheludev/VNom`. Circular clustering: `circuclust`, `github.com/rcedgar/circuclust`.

**Reference material to obtain:** `Supplementary_table_1_stringent_Obelisk_clustering_011724.tsv` (Zheludev et al.) — 1,744 stringent Obelisks, 788 with ≥2 ORFs; ~1,700 Oblin-1 centroid proteins. Marine set: 40 Obelisks (López-Simón et al., *ISME J* 19:wraf033, 2025).

**Positive control:** Obelisk-S.s in *Streptococcus sanguinis* SK36 RNA-seq. Independently reproduced in *J Mol Evol* 93:370 (2025) across 17 SK36 monoculture runs from four independent labs — 13 of which Zheludev et al. did **not** analyze. This is an unusually clean control: multiple independent datasets, known answer.

---

## File structure

```
~/obelisk-hunt/
├── PLAN.md                          # this file
├── config/
│   ├── config.sh                    # all paths and thresholds, sourced by every script
│   └── niche_query.txt              # Entrez query defining the target niche
├── env/
│   └── environment.yml              # conda spec
├── ref/
│   ├── obelisk_stringent.tsv        # raw supplementary table
│   ├── oblin1_centroids.faa         # Oblin-1 proteins → DIAMOND/HMM input
│   ├── oblin1.hmm                   # profile HMM
│   ├── oblin1_dmnd.dmnd             # DIAMOND database
│   ├── obelisk_nt.fna               # known Obelisk nucleotide seqs
│   ├── positive_control.fna         # Obelisk-S.s
│   └── decoy_shuffled.faa           # negative control query
├── scripts/
│   ├── 00_verify_logan.sh
│   ├── 01_fetch_references.sh
│   ├── 02_build_search_dbs.sh
│   ├── 03_positive_control.sh       # GATE 1
│   ├── 04_select_accessions.sh
│   ├── 05_logan_triage.sh
│   ├── 05_logan_triage.sbatch
│   ├── 06_collect_hits.py
│   ├── 07_deep_assemble.sh
│   ├── 07_deep_assemble.sbatch
│   ├── 08_run_vnom.sh
│   ├── 09_circularity_check.py
│   ├── 10_cluster.sh
│   ├── 11_replication_matrix.py     # GATE 3 — contamination control
│   ├── 12_rnafold.sh
│   ├── 13_colabfold_prep.py
│   └── 14_foldseek.sh
├── work/                            # scratch, safe to delete
├── logs/
└── results/
    ├── triage_hits.tsv
    ├── candidates.fna
    ├── validated.fna
    └── figures/
```

**Responsibility boundaries:** `config/config.sh` is the single source of paths and thresholds — no script hardcodes them. `ref/` is built once and never modified. `work/` is disposable. `results/` is the only directory that matters for the poster.

---

# PHASE 0 — Verification and environment

### Task 1: Verify Logan is live and confirm the freeze date

**Files:**
- Create: `~/obelisk-hunt/scripts/00_verify_logan.sh`

- [ ] **Step 1: Make the project root**
      ```bash
      mkdir -p ~/obelisk-hunt/{config,env,ref,scripts,work,logs,results/figures}
      cd ~/obelisk-hunt
      ```

- [ ] **Step 2: Confirm AWS CLI exists**
      Run: `aws --version`
      Expected: `aws-cli/2.x.x ...`
      If missing: `pip install --user awscli` or load your cluster's module.

- [ ] **Step 3: Pull one known accession as a smoke test**
      Run:
      ```bash
      cd ~/obelisk-hunt/work
      aws s3 cp s3://logan-pub/c/DRR000016/DRR000016.contigs.fa.zst . --no-sign-request
      ```
      Expected: a `.zst` file downloads with no credential prompt. `DRR000016` is the accession used in Logan's own tutorial, so it is guaranteed present.
      If this fails with a credentials error, you omitted `--no-sign-request`.

- [ ] **Step 4: Confirm the file decompresses and looks like FASTA**
      Run: `zstdcat DRR000016.contigs.fa.zst | head -4`
      Expected: a `>` header line containing `ka:f:` (k-mer abundance annotation), followed by sequence.
      If `zstdcat` is missing: `conda install -c conda-forge zstd`.

- [ ] **Step 5: Check whether a release newer than v1 (Dec 2023) exists**
      Run:
      ```bash
      aws s3 ls s3://logan-pub/ --no-sign-request
      ```
      Expected: prefixes including `u/`, `c/`, `stats/`. **Record every prefix you see.**
      Then open `https://github.com/IndexThePlanet/Logan` and read the release section.
      **Decision:** if a v2 or newer freeze exists as downloadable data, use it and note the freeze date in `config.sh`. If only v1 exists, proceed with Dec 2023 and state that limitation explicitly in your paper — it is a real constraint, and disclosing it is a credibility gain, not a loss.

- [ ] **Step 6: Download the per-accession stats table**
      Run:
      ```bash
      cd ~/obelisk-hunt/ref
      wget https://s3.amazonaws.com/logan-pub/stats/logan-seqstats.parquet
      ```
      Expected: a parquet file, several hundred MB.

- [ ] **Step 7: Confirm the stats table loads and inspect its columns**
      Run:
      ```bash
      python -c "
      import pandas as pd
      df = pd.read_parquet('$HOME/obelisk-hunt/ref/logan-seqstats.parquet')
      print(df.shape)
      print(df.columns.tolist())
      print(df.head())
      "
      ```
      Expected: roughly 27 million rows; columns including `accession` and contig/unitig size fields.
      This table tells you how large each accession's contigs are **before** you download anything — essential for cost estimation in Task 9.

---

### Task 2: Build the conda environment

**Files:**
- Create: `~/obelisk-hunt/env/environment.yml`

- [ ] **Step 1: Write the environment file**
      ```yaml
      name: obelisk
      channels:
        - conda-forge
        - bioconda
      dependencies:
        - python=3.11
        - pandas
        - pyarrow
        - biopython
        - numpy
        - matplotlib
        - seaborn
        - zstd
        - diamond=2.1.9
        - hmmer=3.4
        - mafft
        - iqtree
        - seqkit
        - samtools
        - minimap2
        - bowtie2
        - spades
        - fastp
        - sra-tools
        - viennarna
        - infernal
        - cd-hit
        - blast
        - entrez-direct
        - pip
        - pip:
          - rcgrep
          - logomaker
      ```

- [ ] **Step 2: Create the environment**
      Run: `conda env create -f ~/obelisk-hunt/env/environment.yml`
      Expected: solves and installs. Takes 10–25 minutes.
      If solving hangs, install `mamba` and use `mamba env create -f ...` instead.

- [ ] **Step 3: Activate and verify every tool resolves**
      Run:
      ```bash
      conda activate obelisk
      for t in diamond hmmbuild hmmsearch mafft iqtree seqkit samtools minimap2 \
               rnaspades.py fastp fasterq-dump RNAfold cmsearch cd-hit blastn esearch; do
        printf '%-16s ' "$t"; command -v "$t" >/dev/null && echo OK || echo MISSING
      done
      ```
      Expected: every line reads `OK`.
      **Do not proceed past any `MISSING`.** Install it individually with `conda install -c bioconda <tool>` before continuing.

- [ ] **Step 4: Install VNom (not on conda)**
      Run:
      ```bash
      cd ~/obelisk-hunt
      git clone https://github.com/Zheludev/VNom.git
      ```
      Expected: clones successfully.
      Read `VNom/README.md` and follow its install instructions. **Record the exact invocation path** — you need it in Task 14.

- [ ] **Step 5: Install circuclust**
      Run:
      ```bash
      cd ~/obelisk-hunt
      git clone https://github.com/rcedgar/circuclust.git
      ```
      Expected: clones successfully. This is a precompiled binary distribution (`v1.0.i86linux64` is the version cited in the literature); confirm the binary runs with `./circuclust/bin/circuclust --version` or equivalent per its README.

- [ ] **Step 6: Write the config file**
      Create `~/obelisk-hunt/config/config.sh`:
      ```bash
      #!/usr/bin/env bash
      # Single source of truth for all paths and thresholds.
      export ROOT="$HOME/obelisk-hunt"
      export REF="$ROOT/ref"
      export WORK="$ROOT/work"
      export LOGS="$ROOT/logs"
      export RESULTS="$ROOT/results"
      export SCRIPTS="$ROOT/scripts"
      export VNOM="$ROOT/VNom"
      export CIRCUCLUST="$ROOT/circuclust"

      # Logan
      export LOGAN_FREEZE="2023-12"        # UPDATE if Task 1 Step 5 found a newer release
      export LOGAN_C="s3://logan-pub/c"
      export LOGAN_U="s3://logan-pub/u"

      # Search thresholds
      export DIAMOND_EVALUE="1e-5"
      export DIAMOND_SENS="--very-sensitive"
      export HMM_EVALUE="1e-10"

      # Validation thresholds
      export MIN_INDEPENDENT_BIOPROJECTS=3
      export OBELISK_MIN_LEN=800
      export OBELISK_MAX_LEN=1200

      # Compute
      export THREADS=16
      ```

- [ ] **Step 7: Verify config sources cleanly**
      Run: `source ~/obelisk-hunt/config/config.sh && echo "$ROOT $LOGAN_C $MIN_INDEPENDENT_BIOPROJECTS"`
      Expected: `/home/<you>/obelisk-hunt s3://logan-pub/c 3`

---

# PHASE 1 — Reference construction

### Task 3: Obtain the published Obelisk reference set

**Files:**
- Create: `~/obelisk-hunt/scripts/01_fetch_references.sh`
- Output: `~/obelisk-hunt/ref/obelisk_stringent.tsv`

- [ ] **Step 1: Download the Zheludev supplementary table**
      Open the *Cell* paper (doi:10.1016/j.cell.2024.09.033) and the bioRxiv preprint (doi:10.1101/2024.01.20.576352). Download **Table S1** / `Supplementary_table_1_stringent_Obelisk_clustering_011724.tsv`.
      Save to `~/obelisk-hunt/ref/obelisk_stringent.tsv`.
      **Verification:** this file is the backbone of the entire project. Confirm before proceeding.

- [ ] **Step 2: Confirm the table's shape and contents**
      Run:
      ```bash
      cd ~/obelisk-hunt/ref
      head -1 obelisk_stringent.tsv | tr '\t' '\n' | cat -n
      wc -l obelisk_stringent.tsv
      ```
      Expected: a column list including sequence and cluster identifiers; roughly 1,744 stringent Obelisk rows (the literature reports 1,744 stringent Obelisks in 1,744 clusters at 80% nt identity, drawn from 29,959 total at 90% identity).
      **Record the actual column names** — the next step depends on them.

- [ ] **Step 3: Extract nucleotide sequences to FASTA**
      Adapt the column names from Step 2 into this script and run it:
      ```python
      # ~/obelisk-hunt/scripts/tsv_to_fasta.py
      import pandas as pd, sys
      tsv, seqcol, idcol, out = sys.argv[1:5]
      df = pd.read_csv(tsv, sep='\t')
      with open(out, 'w') as fh:
          n = 0
          for _, r in df.iterrows():
              s = str(r[seqcol]).strip().upper()
              if s and s != 'NAN':
                  fh.write(f">{r[idcol]}\n{s}\n"); n += 1
      print(f"wrote {n} sequences to {out}")
      ```
      Run:
      ```bash
      python ~/obelisk-hunt/scripts/tsv_to_fasta.py \
        ~/obelisk-hunt/ref/obelisk_stringent.tsv \
        <SEQUENCE_COLUMN_NAME> <ID_COLUMN_NAME> \
        ~/obelisk-hunt/ref/obelisk_nt.fna
      ```
      Expected: `wrote 1744 sequences to .../obelisk_nt.fna` (or close to it).

- [ ] **Step 4: Sanity-check the length distribution**
      Run: `seqkit stats ~/obelisk-hunt/ref/obelisk_nt.fna`
      Expected: average length near 1,000 nt. Obelisks are ~1 kb.
      **If the mean is far from ~1,000, you extracted the wrong column.** Go back to Step 3.

- [ ] **Step 5: Obtain the Oblin-1 protein set**
      The literature references ~1,700 Oblin-1 centroid proteins from Zheludev et al. If the supplementary provides protein sequences directly, extract them the same way into `~/obelisk-hunt/ref/oblin1_centroids.faa`.
      If it provides only nucleotide sequences plus ORF coordinates, translate:
      ```bash
      seqkit subseq --bed <orf_coords.bed> ~/obelisk-hunt/ref/obelisk_nt.fna \
        | seqkit translate --frame 1 --trim \
        > ~/obelisk-hunt/ref/oblin1_centroids.faa
      ```
      Expected: on the order of 1,700 protein sequences.

- [ ] **Step 6: Verify the protein set**
      Run: `seqkit stats ~/obelisk-hunt/ref/oblin1_centroids.faa`
      Expected: hundreds to ~1,700 sequences, mean length roughly 200–300 aa.

- [ ] **Step 7: Add the marine and hot-spring Obelisks**
      Retrieve the 40 marine Obelisks (López-Simón et al., *ISME J* 19:wraf033, 2025) and the hot-spring HsOblin-1 sequences (Urayama et al., *Nat Commun* 17:3041, 2026) from their supplementary data. Append to your reference FASTAs.
      **Why this matters:** including them makes your search *more* sensitive to divergent Oblins, and it guarantees you can recognize and exclude a re-discovery of their findings rather than claiming it as new.

---

### Task 4: Build search databases

**Files:**
- Create: `~/obelisk-hunt/scripts/02_build_search_dbs.sh`

- [ ] **Step 1: Write the script**
      ```bash
      #!/usr/bin/env bash
      set -euo pipefail
      source "$HOME/obelisk-hunt/config/config.sh"

      # DIAMOND protein database for fast translated triage
      diamond makedb --in "$REF/oblin1_centroids.faa" --db "$REF/oblin1_dmnd"

      # Multiple alignment → profile HMM for sensitive confirmation
      mafft --auto --thread "$THREADS" "$REF/oblin1_centroids.faa" > "$REF/oblin1_aln.faa"
      hmmbuild --amino "$REF/oblin1.hmm" "$REF/oblin1_aln.faa"
      hmmpress "$REF/oblin1.hmm"

      # Negative-control decoy: same composition, destroyed motifs
      seqkit shuffle -s 42 "$REF/oblin1_centroids.faa" \
        | seqkit mutate --any-point 100 -s 42 \
        > "$REF/decoy_shuffled.faa" || \
        python - <<'EOF'
      import random
      from Bio import SeqIO
      import os
      random.seed(42)
      ref = os.path.expanduser("~/obelisk-hunt/ref")
      out = []
      for rec in SeqIO.parse(f"{ref}/oblin1_centroids.faa", "fasta"):
          s = list(str(rec.seq)); random.shuffle(s)
          rec.seq = type(rec.seq)("".join(s)); rec.id = "decoy_" + rec.id
          rec.description = ""; out.append(rec)
      SeqIO.write(out, f"{ref}/decoy_shuffled.faa", "fasta")
      print(f"wrote {len(out)} decoys")
      EOF

      diamond makedb --in "$REF/decoy_shuffled.faa" --db "$REF/decoy_dmnd"
      echo "Databases built."
      ```

- [ ] **Step 2: Run it**
      Run: `bash ~/obelisk-hunt/scripts/02_build_search_dbs.sh`
      Expected: ends with `Databases built.`; creates `oblin1_dmnd.dmnd`, `oblin1.hmm`, `decoy_dmnd.dmnd`.

- [ ] **Step 3: Verify the HMM built correctly**
      Run: `head -20 ~/obelisk-hunt/ref/oblin1.hmm`
      Expected: a header block with `HMMER3/f`, a `LENG` line, and `ALPH amino`.

- [ ] **Step 4: Confirm the decoy is genuinely non-matching**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      diamond blastp -q "$REF/decoy_shuffled.faa" -d "$REF/oblin1_dmnd" \
        -e "$DIAMOND_EVALUE" --very-sensitive -o "$WORK/decoy_selfcheck.tsv" -p "$THREADS"
      wc -l "$WORK/decoy_selfcheck.tsv"
      ```
      Expected: **zero or very few lines.** The shuffled decoys must not hit real Oblins at your chosen E-value.
      **If the decoy hits heavily, your E-value threshold is too permissive.** Tighten `DIAMOND_EVALUE` in `config.sh` until the decoy is clean, then re-run. This calibration *is* your false-positive control — do not skip it.

---

# PHASE 2 — Positive control (GATE 1)

> **This is the most important phase in the plan.** If the pipeline cannot re-find a known Obelisk in data where one is known to exist, every downstream result is meaningless. Do not proceed until this passes.

### Task 5: Recover Obelisk-S.s from S. sanguinis SK36

**Files:**
- Create: `~/obelisk-hunt/scripts/03_positive_control.sh`

- [ ] **Step 1: Identify SK36 RNA-seq runs**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      esearch -db sra -query '"Streptococcus sanguinis"[Organism] AND "RNA-Seq"[Strategy]' \
        | efetch -format runinfo > "$WORK/sk36_runinfo.csv"
      wc -l "$WORK/sk36_runinfo.csv"
      ```
      Expected: tens of rows. The *J Mol Evol* 93:370 (2025) paper used 17 SK36 monoculture runs from four independent labs.

- [ ] **Step 2: Pick five runs from at least two different BioProjects**
      Run:
      ```bash
      python - <<'EOF'
      import pandas as pd, os
      w = os.path.expanduser("~/obelisk-hunt/work")
      df = pd.read_csv(f"{w}/sk36_runinfo.csv")
      sel = df.groupby("BioProject").head(3).head(5)
      sel[["Run","BioProject"]].to_csv(f"{w}/pc_runs.tsv", sep="\t", index=False)
      print(sel[["Run","BioProject","ScientificName"]])
      EOF
      ```
      Expected: five run accessions spanning ≥2 BioProjects.
      **Choosing runs from multiple BioProjects is deliberate** — it exercises the same independence logic you rely on later for contamination control.

- [ ] **Step 3: Run the Logan fast path on those accessions**
      ```bash
      #!/usr/bin/env bash
      set -euo pipefail
      source "$HOME/obelisk-hunt/config/config.sh"
      mkdir -p "$WORK/pc"
      tail -n +2 "$WORK/pc_runs.tsv" | cut -f1 | while read -r ACC; do
        echo "=== $ACC ==="
        aws s3 cp "$LOGAN_C/$ACC/$ACC.contigs.fa.zst" - --no-sign-request 2>/dev/null \
          | zstdcat \
          | diamond blastx --db "$REF/oblin1_dmnd" --query - \
              --evalue "$DIAMOND_EVALUE" $DIAMOND_SENS --threads "$THREADS" \
              --outfmt 6 qseqid sseqid pident length evalue bitscore \
              --out "$WORK/pc/$ACC.oblin.tsv" || echo "  no Logan contigs for $ACC"
        [ -f "$WORK/pc/$ACC.oblin.tsv" ] && echo "  hits: $(wc -l < "$WORK/pc/$ACC.oblin.tsv")"
      done
      ```
      Run: `bash ~/obelisk-hunt/scripts/03_positive_control.sh`
      Expected: **at least one accession reports a non-zero hit count.**
      Note the streaming pattern — `aws s3 cp ... -` pipes straight into `zstdcat` and then DIAMOND, so nothing is written to disk. This is the pattern from Logan's own `Chickens.md` tutorial and it is what makes a large sweep affordable.

- [ ] **Step 4: GATE 1 — decide**
      - **Hits found →** the fast path works. Proceed to Step 5.
      - **Zero hits everywhere →** stop and debug in this order: (a) is `oblin1_centroids.faa` real protein sequence? `head` it; (b) does `aws s3 cp` succeed for these accessions, or are they absent from Logan? (c) loosen `DIAMOND_EVALUE` to `1e-3` and retry; (d) if still nothing, the accessions may not be in the Dec 2023 freeze — pick older runs.
      **Do not continue past this gate on a failure.** A silent zero here becomes a false "we found nothing novel" later.

- [ ] **Step 5: Run the deep path on one positive accession**
      ```bash
      source ~/obelisk-hunt/config/config.sh
      ACC=<accession that produced hits>
      cd "$WORK/pc"
      fasterq-dump --split-3 -e "$THREADS" "$ACC"
      fastp -i "${ACC}_1.fastq" -I "${ACC}_2.fastq" \
            -o "${ACC}_1.clean.fastq" -O "${ACC}_2.clean.fastq" \
            --average_qual=30 --n_base_limit=0 --cut_front --cut_tail \
            -j "${ACC}.fastp.json" -h "${ACC}.fastp.html"
      rnaspades.py -1 "${ACC}_1.clean.fastq" -2 "${ACC}_2.clean.fastq" \
                   -o "${ACC}_rnaspades" -t "$THREADS"
      ```
      Expected: `${ACC}_rnaspades/transcripts.fasta` exists and is non-empty.
      These fastp flags are Zheludev et al.'s exact published parameters — do not change them, because matching the reference pipeline is what lets you claim replication.

- [ ] **Step 6: Run VNom on the assembly**
      Run:
      ```bash
      python "$VNOM/VNom.py" -i "$WORK/pc/${ACC}_rnaspades/transcripts.fasta" \
        -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1 \
        -o "$WORK/pc/${ACC}_vnom"
      ```
      (Adjust the entry-point path to match VNom's README from Task 2 Step 4.)
      Expected: an output set of viroid-like circular candidate sequences.

- [ ] **Step 7: Confirm you recovered Obelisk-S.s**
      Run:
      ```bash
      makeblastdb -in "$REF/obelisk_nt.fna" -dbtype nucl -out "$WORK/obelisk_db"
      blastn -query "$WORK/pc/${ACC}_vnom/<vnom_output>.fasta" -db "$WORK/obelisk_db" \
             -outfmt "6 qseqid sseqid pident length evalue" -evalue 1e-20 \
        | sort -k5,5g | head
      ```
      Expected: a high-identity hit (>95%) to the known Obelisk-S.s sequence.
      **This is the moment the pipeline is proven.** Save this output — it is a poster figure and the single most persuasive thing you can show a skeptical judge.

---

# PHASE 3 — Niche selection and triage

### Task 6: Define and enumerate the target niche

**Files:**
- Create: `~/obelisk-hunt/config/niche_query.txt`
- Create: `~/obelisk-hunt/scripts/04_select_accessions.sh`

- [ ] **Step 1: Re-run the incumbent check before committing**
      Search bioRxiv, PubMed, and Google Scholar for: `obelisk rumen`, `obelisk poultry`, `obelisk chicken gut`, `Oblin ruminant`, `viroid-like rumen metatranscriptome`.
      **Record the date and the exact queries in your lab notebook.** Repeat monthly.
      **If your niche is claimed, switch niches — not methods.** Everything downstream is reusable.

- [ ] **Step 2: Write the Entrez query**
      Create `~/obelisk-hunt/config/niche_query.txt`. Rumen example:
      ```
      (rumen[All Fields] OR ruminal[All Fields]) AND "metatranscriptomic"[Source] AND "RNA-Seq"[Strategy]
      ```
      Poultry example:
      ```
      (chicken[Organism] OR poultry[All Fields]) AND (cecum[All Fields] OR caecal[All Fields] OR gut[All Fields]) AND "RNA-Seq"[Strategy]
      ```
      **RNA-Seq strategy is mandatory.** Obelisks are RNA elements; DNA libraries cannot contain them.

- [ ] **Step 3: Enumerate the accessions**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      esearch -db sra -query "$(cat $ROOT/config/niche_query.txt)" \
        | efetch -format runinfo > "$WORK/niche_runinfo.csv"
      wc -l "$WORK/niche_runinfo.csv"
      cut -d, -f22 "$WORK/niche_runinfo.csv" | sort -u | head   # BioProject column
      ```
      Expected: hundreds to tens of thousands of runs across many BioProjects.
      **Sanity thresholds:** under ~200 runs, the niche is too small — broaden it. Over ~200,000, narrow it or you will not finish the sweep.

- [ ] **Step 4: Extract the run list and the run→BioProject map**
      Run:
      ```bash
      python - <<'EOF'
      import pandas as pd, os
      w = os.path.expanduser("~/obelisk-hunt/work")
      df = pd.read_csv(f"{w}/niche_runinfo.csv")
      df = df[df["Run"].notna()]
      df["Run"].to_csv(f"{w}/accessions.txt", index=False, header=False)
      df[["Run","BioProject","LibraryStrategy","Platform","ScientificName"]].to_csv(
          f"{w}/run_metadata.tsv", sep="\t", index=False)
      print("runs:", len(df), "bioprojects:", df["BioProject"].nunique())
      EOF
      ```
      Expected: prints both counts. **You need many BioProjects** — replication across independent projects is the whole contamination defense.

- [ ] **Step 5: Estimate download volume before spending anything**
      Run:
      ```bash
      python - <<'EOF'
      import pandas as pd, os
      r = os.path.expanduser("~/obelisk-hunt/ref"); w = os.path.expanduser("~/obelisk-hunt/work")
      stats = pd.read_parquet(f"{r}/logan-seqstats.parquet")
      accs = [l.strip() for l in open(f"{w}/accessions.txt") if l.strip()]
      sub = stats[stats["accession"].isin(accs)]
      col = [c for c in sub.columns if "contig" in c and "after" in c]
      print("matched in Logan:", len(sub), "of", len(accs))
      if col:
          tb = sub[col[0]].sum()/1e12
          print(f"total compressed contigs: {tb:.2f} TB")
      EOF
      ```
      Expected: a coverage fraction and a TB estimate.
      **Two things to read here.** The match fraction tells you how much of your niche predates the Dec 2023 freeze — anything missing is invisible to the fast path. The TB figure sizes the sweep; because you stream rather than store, this is transfer volume, not disk requirement.

---

### Task 7: Run the Logan triage sweep

**Files:**
- Create: `~/obelisk-hunt/scripts/05_logan_triage.sh`
- Create: `~/obelisk-hunt/scripts/05_logan_triage.sbatch`

- [ ] **Step 1: Write the per-accession worker**
      ```bash
      #!/usr/bin/env bash
      # Usage: 05_logan_triage.sh <ACCESSION>
      set -uo pipefail
      source "$HOME/obelisk-hunt/config/config.sh"
      ACC="$1"
      OUT="$WORK/triage/$ACC.tsv"
      mkdir -p "$WORK/triage"
      [ -s "$OUT" ] && { echo "$ACC already done"; exit 0; }   # idempotent

      aws s3 cp "$LOGAN_C/$ACC/$ACC.contigs.fa.zst" - --no-sign-request 2>/dev/null \
        | zstdcat \
        | diamond blastx --db "$REF/oblin1_dmnd" --query - \
            --evalue "$DIAMOND_EVALUE" $DIAMOND_SENS --threads 2 \
            --outfmt 6 qseqid sseqid pident length evalue bitscore qseq \
            --out "$OUT" 2>>"$LOGS/triage.err"

      if [ -s "$OUT" ]; then
        awk -v a="$ACC" '{print a"\t"$0}' "$OUT" >> "$RESULTS/triage_hits.tsv"
        echo "$ACC HIT $(wc -l < "$OUT")"
      else
        touch "$OUT"; echo "$ACC none"
      fi
      ```
      Note `--outfmt` includes `qseq` — you keep the matching sequence itself, so a hit is immediately usable without re-downloading.
      The `[ -s "$OUT" ] && exit 0` guard makes the script idempotent: re-running after a job failure skips completed accessions instead of redoing them.

- [ ] **Step 2: Test on ten accessions before scaling**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      head -10 "$WORK/accessions.txt" | while read -r A; do
        bash "$SCRIPTS/05_logan_triage.sh" "$A"
      done
      ```
      Expected: ten lines, each `<ACC> HIT <n>` or `<ACC> none`.
      **Time this step.** Multiply by your accession count to project total runtime before submitting anything large.

- [ ] **Step 3: Write the SLURM array script**
      ```bash
      #!/usr/bin/env bash
      #SBATCH --job-name=logan-triage
      #SBATCH --array=1-1000%100
      #SBATCH --cpus-per-task=2
      #SBATCH --mem=8G
      #SBATCH --time=02:00:00
      #SBATCH --output=%x-%A_%a.out
      # ADJUST: --partition and --account for your cluster.
      # ADJUST: --array upper bound to your accession count (see Step 4).
      # %100 caps concurrent tasks — raise or lower to respect S3 egress and local policy.

      source "$HOME/obelisk-hunt/config/config.sh"
      source activate obelisk 2>/dev/null || conda activate obelisk

      CHUNK=100
      START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
      END=$(( SLURM_ARRAY_TASK_ID * CHUNK ))
      sed -n "${START},${END}p" "$WORK/accessions.txt" | while read -r A; do
        [ -n "$A" ] && bash "$SCRIPTS/05_logan_triage.sh" "$A"
      done
      ```

- [ ] **Step 4: Compute the correct array bound**
      Run:
      ```bash
      N=$(wc -l < ~/obelisk-hunt/work/accessions.txt)
      echo "accessions: $N ; array tasks needed: $(( (N + 99) / 100 ))"
      ```
      Edit `--array=1-<that number>%100` in the sbatch file.

- [ ] **Step 5: Submit a 5-task pilot**
      Run: `sbatch --array=1-5%5 ~/obelisk-hunt/scripts/05_logan_triage.sbatch`
      Expected: a job ID. Watch with `squeue -u $USER`.
      Confirm the `.out` files show per-accession lines and no repeated S3 or conda errors.
      **Never submit the full array before the pilot is clean.**

- [ ] **Step 6: Submit the full sweep**
      Run: `sbatch ~/obelisk-hunt/scripts/05_logan_triage.sbatch`

- [ ] **Step 7: Monitor**
      Run:
      ```bash
      echo "done: $(ls ~/obelisk-hunt/work/triage | wc -l)"
      echo "hits: $(cut -f1 ~/obelisk-hunt/results/triage_hits.tsv 2>/dev/null | sort -u | wc -l)"
      tail -20 ~/obelisk-hunt/logs/triage.err
      ```
      Expected: `done` climbs steadily; `hits` is a small fraction of it.
      **If the hit rate exceeds ~20%, stop.** That is far above biological plausibility and means your threshold is too loose or the decoy calibration in Task 4 Step 4 failed. Re-tighten and re-run.

---

### Task 8: Confirm triage hits with the profile HMM

**Files:**
- Create: `~/obelisk-hunt/scripts/06_collect_hits.py`

- [ ] **Step 1: Write the collector**
      ```python
      #!/usr/bin/env python3
      """Collapse DIAMOND triage hits to per-accession, HMM-confirm, emit ranked list."""
      import os, subprocess, pandas as pd
      from collections import defaultdict

      ROOT = os.path.expanduser("~/obelisk-hunt")
      hits_f = f"{ROOT}/results/triage_hits.tsv"
      cols = ["accession","qseqid","sseqid","pident","length","evalue","bitscore","qseq"]
      df = pd.read_csv(hits_f, sep="\t", names=cols)
      print(f"raw hit rows: {len(df)}  accessions: {df.accession.nunique()}")

      # Write best hit per accession as protein for HMM confirmation
      best = df.sort_values("evalue").groupby("accession").first().reset_index()
      faa = f"{ROOT}/work/triage_best.faa"
      with open(faa, "w") as fh:
          for _, r in best.iterrows():
              fh.write(f">{r.accession}\n{r.qseq}\n")
      print(f"wrote {len(best)} best-hit sequences")

      # HMM confirmation — independent of DIAMOND, so it is a real second filter
      tbl = f"{ROOT}/work/triage_hmm.tbl"
      subprocess.run(["hmmsearch","--tblout",tbl,"-E",os.environ.get("HMM_EVALUE","1e-10"),
                      f"{ROOT}/ref/oblin1.hmm", faa], check=True,
                     stdout=subprocess.DEVNULL)

      confirmed = set()
      for line in open(tbl):
          if not line.startswith("#"):
              confirmed.add(line.split()[0])
      print(f"HMM-confirmed accessions: {len(confirmed)}")

      # Attach BioProject so independence can be assessed immediately
      meta = pd.read_csv(f"{ROOT}/work/run_metadata.tsv", sep="\t")
      out = best[best.accession.isin(confirmed)].merge(
          meta, left_on="accession", right_on="Run", how="left")
      out.to_csv(f"{ROOT}/results/confirmed_hits.tsv", sep="\t", index=False)

      per_bp = defaultdict(list)
      for _, r in out.iterrows():
          per_bp[r.get("BioProject","NA")].append(r.accession)
      print(f"independent BioProjects with hits: {len(per_bp)}")
      for bp, accs in sorted(per_bp.items(), key=lambda x: -len(x[1]))[:10]:
          print(f"  {bp}: {len(accs)}")
      ```

- [ ] **Step 2: Run it**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      python ~/obelisk-hunt/scripts/06_collect_hits.py
      ```
      Expected: prints raw rows, best-hit count, HMM-confirmed count, and BioProject spread.

- [ ] **Step 3: GATE 2 — assess the BioProject spread**
      - **Hits in ≥3 independent BioProjects →** proceed to Phase 4.
      - **Hits in only 1–2 →** treat as a contamination warning, not a discovery. Broaden the niche (Task 6) and re-sweep.
      - **Zero confirmed hits →** either the niche genuinely lacks Obelisks (a real, reportable negative result) or the threshold is too strict. Loosen `HMM_EVALUE` to `1e-5` and re-run Step 2 before concluding.

---

# PHASE 4 — Deep assembly on survivors

### Task 9: Full Zheludev pipeline on confirmed accessions

**Files:**
- Create: `~/obelisk-hunt/scripts/07_deep_assemble.sh`
- Create: `~/obelisk-hunt/scripts/07_deep_assemble.sbatch`

- [ ] **Step 1: Write the worker**
      ```bash
      #!/usr/bin/env bash
      # Usage: 07_deep_assemble.sh <ACCESSION>
      set -uo pipefail
      source "$HOME/obelisk-hunt/config/config.sh"
      ACC="$1"; D="$WORK/deep/$ACC"
      [ -s "$D/transcripts.fasta" ] && { echo "$ACC already assembled"; exit 0; }
      mkdir -p "$D"; cd "$D"

      fasterq-dump --split-3 -e "$THREADS" -O . "$ACC" || { echo "$ACC fasterq FAIL"; exit 1; }

      if [ -f "${ACC}_1.fastq" ] && [ -f "${ACC}_2.fastq" ]; then
        fastp -i "${ACC}_1.fastq" -I "${ACC}_2.fastq" \
              -o r1.clean.fastq -O r2.clean.fastq \
              --average_qual=30 --n_base_limit=0 --cut_front --cut_tail \
              -j fastp.json -h fastp.html
        rnaspades.py -1 r1.clean.fastq -2 r2.clean.fastq -o . -t "$THREADS" -m 64
      else
        fastp -i "${ACC}.fastq" -o r.clean.fastq \
              --average_qual=30 --n_base_limit=0 --cut_front --cut_tail \
              -j fastp.json -h fastp.html
        rnaspades.py -s r.clean.fastq -o . -t "$THREADS" -m 64
      fi

      rm -f ./*.fastq                          # raw reads are large; assemblies are what matter
      [ -s transcripts.fasta ] && echo "$ACC OK $(grep -c '^>' transcripts.fasta) transcripts" \
                               || echo "$ACC ASSEMBLY FAIL"
      ```

- [ ] **Step 2: Build the deep-path accession list**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      cut -f1 "$RESULTS/confirmed_hits.tsv" | tail -n +2 | sort -u > "$WORK/deep_accessions.txt"
      wc -l "$WORK/deep_accessions.txt"
      ```
      Expected: far fewer than the triage list — typically tens to a few hundred. **That reduction is the entire point of the hybrid architecture.**

- [ ] **Step 3: Test on one accession**
      Run: `bash ~/obelisk-hunt/scripts/07_deep_assemble.sh $(head -1 ~/obelisk-hunt/work/deep_accessions.txt)`
      Expected: `<ACC> OK <n> transcripts`.
      **Time and measure peak memory.** rnaSPAdes on a large metatranscriptome can exceed 64 GB; raise `-m` and `--mem` together if it fails.

- [ ] **Step 4: Write the SLURM array**
      ```bash
      #!/usr/bin/env bash
      #SBATCH --job-name=deep-assemble
      #SBATCH --array=1-100%20
      #SBATCH --cpus-per-task=16
      #SBATCH --mem=64G
      #SBATCH --time=12:00:00
      #SBATCH --output=%x-%A_%a.out
      # ADJUST: --partition, --account, --array bound (= line count of deep_accessions.txt)

      source "$HOME/obelisk-hunt/config/config.sh"
      source activate obelisk 2>/dev/null || conda activate obelisk
      ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$WORK/deep_accessions.txt")
      [ -n "$ACC" ] && bash "$SCRIPTS/07_deep_assemble.sh" "$ACC"
      ```

- [ ] **Step 5: Set the array bound and submit**
      Run:
      ```bash
      N=$(wc -l < ~/obelisk-hunt/work/deep_accessions.txt)
      sed -i "s/#SBATCH --array=1-100%20/#SBATCH --array=1-${N}%20/" \
        ~/obelisk-hunt/scripts/07_deep_assemble.sbatch
      sbatch ~/obelisk-hunt/scripts/07_deep_assemble.sbatch
      ```

- [ ] **Step 6: Verify completion**
      Run:
      ```bash
      ls ~/obelisk-hunt/work/deep/*/transcripts.fasta 2>/dev/null | wc -l
      grep -l "ASSEMBLY FAIL" ~/obelisk-hunt/deep-assemble-*.out 2>/dev/null | wc -l
      ```
      Expected: assembly count close to the accession count; few failures.

---

### Task 10: Run VNom to call circular candidates

**Files:**
- Create: `~/obelisk-hunt/scripts/08_run_vnom.sh`

- [ ] **Step 1: Write the script**
      ```bash
      #!/usr/bin/env bash
      set -euo pipefail
      source "$HOME/obelisk-hunt/config/config.sh"
      mkdir -p "$WORK/vnom"
      for T in "$WORK"/deep/*/transcripts.fasta; do
        ACC=$(basename "$(dirname "$T")")
        [ -d "$WORK/vnom/$ACC" ] && continue
        echo "=== VNom $ACC ==="
        python "$VNOM/VNom.py" -i "$T" \
          -max 2000 -CF_k 10 -CF_simple 0 -CF_tandem 1 -USG_vs_all 1 \
          -o "$WORK/vnom/$ACC" 2>>"$LOGS/vnom.err" || echo "  VNom failed on $ACC"
      done
      echo "VNom complete: $(ls -d "$WORK"/vnom/*/ 2>/dev/null | wc -l) accessions"
      ```
      The flags are Zheludev et al.'s published values verbatim. Adjust only the entry-point path to match VNom's README.

- [ ] **Step 2: Run it**
      Run: `bash ~/obelisk-hunt/scripts/08_run_vnom.sh`
      Expected: ends with a count of processed accessions.

- [ ] **Step 3: Pool candidates and filter to Obelisk size range**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      cat "$WORK"/vnom/*/*.fasta 2>/dev/null \
        | seqkit seq -m "$OBELISK_MIN_LEN" -M "$OBELISK_MAX_LEN" \
        | seqkit rmdup -s \
        > "$RESULTS/candidates.fna"
      seqkit stats "$RESULTS/candidates.fna"
      ```
      Expected: a candidate set with mean length near 1,000 nt.
      `rmdup -s` collapses identical sequences so one element recovered from twenty libraries counts once.

---

# PHASE 5 — Validation (GATE 3)

> Everything up to here produces *candidates*. This phase produces *evidence*. Judges and reviewers will attack precisely here.

### Task 11: Circularity and ORF confirmation

**Files:**
- Create: `~/obelisk-hunt/scripts/09_circularity_check.py`

- [ ] **Step 1: Write the checker**
      ```python
      #!/usr/bin/env python3
      """Confirm head-to-tail circularity and find Oblin-like ORFs."""
      import os, subprocess
      from Bio import SeqIO
      from Bio.Seq import Seq

      ROOT = os.path.expanduser("~/obelisk-hunt")
      MIN_OVERLAP = 20

      def is_circular(seq, k=MIN_OVERLAP):
          """A circular element assembled linearly repeats its start at its end."""
          s = str(seq).upper()
          return len(s) > 2*k and s[:k] == s[-k:]

      def longest_orf(seq):
          best = ""
          for strand, nt in [(1, seq), (-1, seq.reverse_complement())]:
              for frame in range(3):
                  prot = str(nt[frame:].translate(to_stop=False))
                  for piece in prot.split("*"):
                      if "M" in piece:
                          cand = piece[piece.index("M"):]
                          if len(cand) > len(best): best = cand
          return best

      rows, faa = [], []
      for rec in SeqIO.parse(f"{ROOT}/results/candidates.fna", "fasta"):
          circ = is_circular(rec.seq)
          orf = longest_orf(rec.seq)
          rows.append((rec.id, len(rec.seq), circ, len(orf)))
          if len(orf) >= 100:
              faa.append(f">{rec.id}\n{orf}\n")

      with open(f"{ROOT}/work/candidate_orfs.faa","w") as fh: fh.writelines(faa)
      with open(f"{ROOT}/results/circularity.tsv","w") as fh:
          fh.write("id\tlength\tcircular\torf_aa\n")
          for r in rows: fh.write("\t".join(map(str,r))+"\n")

      print(f"candidates: {len(rows)}")
      print(f"circular:   {sum(1 for r in rows if r[2])}")
      print(f"ORF>=100aa: {len(faa)}")
      ```
      Note: `is_circular` uses terminal-repeat detection, which is a *heuristic*. Report it as supporting evidence, and confirm true circularity by read-mapping across the junction (Task 12 Step 4).

- [ ] **Step 2: Run it**
      Run: `python ~/obelisk-hunt/scripts/09_circularity_check.py`
      Expected: three counts printed.

- [ ] **Step 3: Confirm ORFs are Oblin-like, not something else**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      hmmsearch --tblout "$WORK/candidate_oblin.tbl" -E 1e-5 \
        "$REF/oblin1.hmm" "$WORK/candidate_orfs.faa" > /dev/null
      grep -vc '^#' "$WORK/candidate_oblin.tbl"
      ```
      Expected: a count of ORFs matching the Oblin-1 profile.

- [ ] **Step 4: Exclude known Obelisks — this defines what is actually new**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      makeblastdb -in "$REF/obelisk_nt.fna" -dbtype nucl -out "$WORK/known_db"
      blastn -query "$RESULTS/candidates.fna" -db "$WORK/known_db" \
             -outfmt "6 qseqid sseqid pident length evalue" -evalue 1e-10 \
             -out "$WORK/vs_known.tsv"
      awk '$3 >= 95' "$WORK/vs_known.tsv" | cut -f1 | sort -u > "$WORK/already_known.txt"
      wc -l "$WORK/already_known.txt"
      seqkit grep -v -f "$WORK/already_known.txt" "$RESULTS/candidates.fna" \
        > "$RESULTS/candidates_novel.fna"
      seqkit stats "$RESULTS/candidates_novel.fna"
      ```
      Expected: a count of re-discoveries, and a novel-only FASTA.
      **Re-discovering known Obelisks is a good sign, not a bad one** — it independently validates the pipeline. Report both numbers; the ratio is itself a credibility statistic.

---

### Task 12: The replication matrix — primary contamination control

**Files:**
- Create: `~/obelisk-hunt/scripts/11_replication_matrix.py`

- [ ] **Step 1: Write the analysis**
      ```python
      #!/usr/bin/env python3
      """Independence evidence: which candidates appear across unrelated BioProjects,
      sequencing centers, and platforms. This is the core anti-contamination control."""
      import os, subprocess, pandas as pd
      from collections import defaultdict

      ROOT = os.path.expanduser("~/obelisk-hunt")
      MIN_BP = int(os.environ.get("MIN_INDEPENDENT_BIOPROJECTS", 3))

      # Map each novel candidate back to every accession whose VNom output contained it
      subprocess.run(["makeblastdb","-in",f"{ROOT}/results/candidates_novel.fna",
                      "-dbtype","nucl","-out",f"{ROOT}/work/novel_db"],
                     check=True, stdout=subprocess.DEVNULL)

      occ = defaultdict(set)
      for d in os.listdir(f"{ROOT}/work/vnom"):
          path = f"{ROOT}/work/vnom/{d}"
          if not os.path.isdir(path): continue
          fas = [f for f in os.listdir(path) if f.endswith(".fasta")]
          if not fas: continue
          out = subprocess.run(
              ["blastn","-query",f"{path}/{fas[0]}","-db",f"{ROOT}/work/novel_db",
               "-outfmt","6 qseqid sseqid pident","-evalue","1e-20","-perc_identity","95"],
              capture_output=True, text=True).stdout
          for line in out.strip().split("\n"):
              if line: occ[line.split("\t")[1]].add(d)

      meta = pd.read_csv(f"{ROOT}/work/run_metadata.tsv", sep="\t").set_index("Run")
      rows = []
      for cand, accs in occ.items():
          bps, cents, plats = set(), set(), set()
          for a in accs:
              if a in meta.index:
                  r = meta.loc[a]
                  bps.add(r.get("BioProject","NA"))
                  cents.add(r.get("CenterName","NA"))
                  plats.add(r.get("Platform","NA"))
          rows.append({"candidate":cand,"n_accessions":len(accs),
                       "n_bioprojects":len(bps),"n_centers":len(cents),
                       "n_platforms":len(plats),
                       "bioprojects":";".join(sorted(map(str,bps)))})

      df = pd.DataFrame(rows).sort_values("n_bioprojects", ascending=False)
      df.to_csv(f"{ROOT}/results/replication_matrix.tsv", sep="\t", index=False)

      passed = df[df.n_bioprojects >= MIN_BP]
      print(f"candidates examined: {len(df)}")
      print(f"passing >={MIN_BP} independent BioProjects: {len(passed)}")
      print(passed.head(20).to_string(index=False))

      passed[["candidate"]].to_csv(f"{ROOT}/work/passed_ids.txt", index=False, header=False)
      ```

- [ ] **Step 2: Run it**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      python ~/obelisk-hunt/scripts/11_replication_matrix.py
      ```
      Expected: a table ranked by independent BioProject count.

- [ ] **Step 3: Emit the validated set**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      seqkit grep -f "$WORK/passed_ids.txt" "$RESULTS/candidates_novel.fna" \
        > "$RESULTS/validated.fna"
      seqkit stats "$RESULTS/validated.fna"
      ```
      Expected: your validated discovery set. **This file is the project.**

- [ ] **Step 4: Confirm circularity by junction read-mapping**
      For each of your top candidates:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      CAND=<candidate_id>; ACC=<an accession containing it>
      cd "$WORK"
      seqkit grep -p "$CAND" "$RESULTS/validated.fna" > cand.fna
      # Build a doubled reference so reads spanning the junction can align
      python -c "
      from Bio import SeqIO
      r = next(SeqIO.parse('cand.fna','fasta'))
      print('>'+r.id+'_doubled'); print(str(r.seq)*2)
      " > cand_doubled.fna
      bowtie2-build cand_doubled.fna cand_idx
      fasterq-dump --split-3 -e "$THREADS" -O . "$ACC"
      bowtie2 -x cand_idx -1 "${ACC}_1.fastq" -2 "${ACC}_2.fastq" -p "$THREADS" \
        | samtools sort -o cand.bam -
      samtools index cand.bam
      samtools depth -a cand.bam | awk '{s+=$3; n++} END {print "mean depth:", s/n}'
      ```
      Expected: reads mapping **across the junction point** (position ≈ L, the original length) in the doubled reference. That is direct physical evidence of circularity, not a heuristic.
      **This is the single most defensible validation step in the plan.** Make it a poster figure.

- [ ] **Step 5: Screen against known contaminants**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      blastn -query "$RESULTS/validated.fna" -db nt -remote \
             -outfmt "6 qseqid sseqid pident length evalue stitle" \
             -evalue 1e-10 -max_target_seqs 5 > "$WORK/validated_vs_nt.tsv"
      cut -f6 "$WORK/validated_vs_nt.tsv" | sort | uniq -c | sort -rn | head -20
      ```
      Expected: no strong hits to vectors, PhiX, adapters, or common lab strains.
      **Any candidate hitting a cloning vector or PhiX is a contaminant. Remove it and say so in your paper.** Documenting removals is a credibility gain.

- [ ] **Step 6: GATE 3 — decide**
      - **≥1 candidate passing ≥3 BioProjects, circular by junction mapping, Oblin-1 HMM positive, no contaminant hits →** proceed to Phase 6. You have a discovery.
      - **Candidates only in 1–2 BioProjects →** report as *tentative*, explicitly flagged. Do not overclaim.
      - **Nothing passes →** this is a **publishable negative result**: "systematic search of N accessions across M BioProjects in [niche] using a pipeline validated on a positive control found no novel Obelisk-like elements." The positive control from Task 5 is what makes that statement meaningful rather than an admission of failure. This is a legitimate ISEF project — say so plainly.

---

# PHASE 6 — Structural characterization

### Task 13: RNA secondary structure

- [ ] **Step 1: Fold every validated candidate**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      mkdir -p "$RESULTS/figures/rnafold" && cd "$RESULTS/figures/rnafold"
      RNAfold --noPS < "$RESULTS/validated.fna" > validated_fold.txt
      grep -A2 '^>' validated_fold.txt | head -30
      ```
      Expected: dot-bracket structures with minimum free energy values.

- [ ] **Step 2: Generate structure plots for the top candidate**
      Run:
      ```bash
      seqkit grep -p <top_candidate_id> "$RESULTS/validated.fna" | RNAfold
      ```
      Expected: an `.ps` file per sequence. Convert with `ps2pdf`.
      **Look for the rod-like fold.** Obelisks are defined partly by this morphology — a genuine rod-like structure is strong corroborating evidence and a striking visual.

- [ ] **Step 3: Search for ribozyme motifs**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      cd "$WORK"
      wget -q https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
      gunzip -f Rfam.cm.gz && cmpress Rfam.cm
      cmscan --tblout validated_rfam.tbl --cut_ga Rfam.cm "$RESULTS/validated.fna" > /dev/null
      grep -v '^#' validated_rfam.tbl | awk '{print $2, $3, $16}' | sort -u | head
      ```
      Expected: possible hits to hammerhead ribozyme families. Some Obelisks carry them; absence is not disqualifying.

---

### Task 14: Protein structure — the differentiator

- [ ] **Step 1: Extract Oblin-like ORFs from validated candidates**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      cut -f1 "$WORK/passed_ids.txt" > "$WORK/keep.txt"
      seqkit grep -f "$WORK/keep.txt" "$WORK/candidate_orfs.faa" \
        > "$RESULTS/validated_oblins.faa"
      seqkit stats "$RESULTS/validated_oblins.faa"
      ```

- [ ] **Step 2: Predict structures with ColabFold**
      Open ColabFold in Google Colab (`AlphaFold2.ipynb` from the ColabFold repository). Paste each Oblin-like protein sequence. Run with default settings.
      Expected: a predicted PDB and a pLDDT score.
      **Record pLDDT for every model.** Zheludev et al. reported a mean pLDDT of 83.8 for Oblin-1; treat anything below ~70 as low confidence and say so.
      Save models to `~/obelisk-hunt/results/figures/colabfold/`.

- [ ] **Step 3: Structural homology search with Foldseek**
      Go to the Foldseek web server, upload each predicted PDB, and search against AFDB and PDB.
      Expected, and both outcomes are good:
      - **No significant hits** → a genuinely novel fold. The hot-spring paper reported a best E-value of 0.31, i.e. nothing significant, and used exactly that to argue novelty.
      - **Hits to known Oblin-1** → confirms your element is a true Obelisk relative.
      **Record the top five hits and E-values for each.** This table goes on the poster.

- [ ] **Step 4: Build the phylogeny**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      cat "$RESULTS/validated_oblins.faa" "$REF/oblin1_centroids.faa" > "$WORK/tree_in.faa"
      mafft --auto --thread "$THREADS" "$WORK/tree_in.faa" > "$WORK/tree_aln.faa"
      iqtree -s "$WORK/tree_aln.faa" -m MFP -bb 1000 -nt "$THREADS" \
             -pre "$RESULTS/figures/oblin_tree"
      ```
      Expected: a treefile with bootstrap support.
      **Where your sequences fall matters.** A well-supported clade separate from all published Oblins is the phylogenetic version of "this is new."

- [ ] **Step 5: Cluster with circuclust for nomenclature**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      cat "$RESULTS/validated.fna" "$REF/obelisk_nt.fna" > "$WORK/all_obelisks.fna"
      "$CIRCUCLUST/bin/circuclust" --input "$WORK/all_obelisks.fna" --id 0.80 \
        --output "$RESULTS/clusters80.tsv"
      ```
      (Confirm exact flags against circuclust's README.)
      The field's convention is `Obelisk_X_Y_Z` — X = cluster at 80% nt identity, Y = at 95%, Z = strain. Following it makes your results directly comparable to published work and signals you know the literature.

---

# PHASE 7 — Package for wet-lab handoff

### Task 15: Build the confirmation request

- [ ] **Step 1: Rank candidates for RT-PCR**
      Run:
      ```bash
      source ~/obelisk-hunt/config/config.sh
      python - <<'EOF'
      import pandas as pd, os
      R = os.path.expanduser("~/obelisk-hunt/results")
      df = pd.read_csv(f"{R}/replication_matrix.tsv", sep="\t")
      df["score"] = df.n_bioprojects*3 + df.n_centers*2 + df.n_platforms
      top = df.sort_values("score", ascending=False).head(5)
      top.to_csv(f"{R}/wetlab_priority.tsv", sep="\t", index=False)
      print(top.to_string(index=False))
      EOF
      ```
      Expected: five ranked candidates. Prioritizing by *independence* rather than abundance is deliberate — the most-replicated element is the least likely to be an artifact.

- [ ] **Step 2: Design primers**
      For each of the top three, design outward-facing (divergent) primers — these amplify only across the circular junction and produce no product from a linear template. That is the standard RT-PCR test for circularity. Use Primer3 or Benchling; target 100–300 bp products.
      Save to `~/obelisk-hunt/results/primers.tsv`.

- [ ] **Step 3: Identify obtainable sample sources**
      From `results/wetlab_priority.tsv`, look up each candidate's source BioProjects and record the sample type and origin. Note which are plausibly obtainable locally (a rumen sample from a veterinary school, a poultry sample from an agricultural program).

- [ ] **Step 4: Write the one-page brief**
      Create `~/obelisk-hunt/results/WETLAB_BRIEF.md` containing: what an Obelisk is (three sentences); what you found and the evidence tier for each candidate; the exact RT-PCR you are requesting; the primers; the sample type needed; and what a positive result would establish.
      **Bring this, not an idea.** A one-page request with validated targets and designed primers is a far easier yes than a conversation about a possibility.

---

## Self-review

Run against the spec.

**1. Spec coverage.**

| Requirement | Task |
|---|---|
| Verify Logan works | Task 1 |
| Environment / tooling | Task 2 |
| Reference data | Task 3 |
| Search databases | Task 4 |
| Positive control | Task 5 (GATE 1) |
| Niche definition | Task 6 |
| Large-scale triage on HPC | Task 7 |
| Confirmation filter | Task 8 (GATE 2) |
| Deep assembly | Task 9 |
| Candidate calling | Task 10 |
| Circularity + ORFs | Task 11 |
| Contamination control | Task 12 (GATE 3) |
| RNA structure | Task 13 |
| Protein structure / Foldseek | Task 14 |
| Wet-lab packaging | Task 15 |
| Scoop monitoring | Task 6 Step 1 (monthly) |

No gaps.

**2. Placeholder scan.** Four items require values only you can supply, each with an explicit resolution step rather than a silent blank: column names in Task 3 Step 3 (resolved by Step 2); VNom's entry path in Task 10 (resolved by Task 2 Step 4); SLURM partition/account in Tasks 7 and 9 (marked `ADJUST`); circuclust flags in Task 14 Step 5 (verify against README). No "TBD," no "add error handling," no undefined references.

**3. Naming consistency.** Verified across tasks: `$REF`/`$WORK`/`$RESULTS`/`$SCRIPTS`/`$VNOM`/`$CIRCUCLUST` defined once in `config.sh`; `oblin1_centroids.faa` → `oblin1_dmnd` → `oblin1.hmm` chain consistent; `candidates.fna` → `candidates_novel.fna` → `validated.fna` chain consistent; `passed_ids.txt` written in Task 12 and consumed in Tasks 12 Step 3 and 14 Step 1.

**4. Gate logic.** Three gates, each with an explicit failure branch, and each failure branch leads to a defined action rather than a dead end. GATE 3 failure yields a publishable negative result *because* GATE 1 passed — that dependency is the plan's most important structural property.

---

## Known unknowns

Flagged rather than papered over:

1. **Logan freeze date.** Repo says v1 = Dec 2023. My earlier claim of a Dec 2025 / 87-petabase release is **unverified**. Task 1 Step 5 resolves this. If only v1 exists, anything deposited in SRA after Dec 2023 is invisible to the fast path — a real limitation to state in your paper.
2. **Supplementary table structure.** Column names in Zheludev's Table S1 are unverified; Task 3 Step 2 discovers them before use.
3. **Whether Oblin-1 protein sequences are directly provided** or must be translated from ORF coordinates. Task 3 Step 5 handles both.
4. **VNom's exact CLI.** The *parameters* are from the published methods and are reliable; the *invocation path* must come from its README.
5. **circuclust flags.** Version `v1.0.i86linux64` is cited in the literature; verify flags locally.
6. **Cluster specifics.** Partition names, account strings, and module systems vary; marked `ADJUST`.
7. **Skills registry.** `npx skills find` returned nothing for six bioinformatics queries **and also nothing for the control query `react`**, which per the skill's own documentation should have 100K+ install results. The registry is unreachable from this sandbox (`skills.sh` is outside the allowed-domain list), so **no conclusion can be drawn about what bioinformatics skills exist.** Run `npx skills find bioinformatics` on your own machine.

---

## Execution handoff

Two options:

1. **Subagent-driven (recommended)** — a fresh subagent per task, with review between tasks. Best fit here: tasks are independent, and the gates give natural review points.
2. **Inline** — execute in one session with checkpoints at each gate.

Either way, **do not skip Task 5.** Every downstream claim depends on the positive control passing first.
