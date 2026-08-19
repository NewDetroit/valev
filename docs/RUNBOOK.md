# RUNBOOK

Operator guide. Read `docs/CORRECTIONS.md` first if you are about to change anything —
several scripts deliberately disagree with `docs/PLAN.original.md`.

---

## 0. Before you start

**On a machine with unrestricted network**, do these three things once. Each resolves
something the build environment could not.

1. **Check for a newer Logan release.**
   ```bash
   bash scripts/00_verify_logan.sh
   ```
   Exits non-zero if Logan's README stops advertising v1.2. If a v1.3 appears, add a case to
   `config/config.sh` — `LOGAN_RELEASE` is the only place that needs to change.

2. **Resolve the skills-registry question** (`OPEN-B`) — it failed from the build sandbox
   *and* its control query failed, so nothing can be concluded from that attempt:
   ```bash
   npx skills find bioinformatics
   npx skills find react          # control: if this is also empty, the registry is unreachable
   ```

3. **Read the two literature sources that carry unverified numbers** (`OPEN-C`). Do not put
   the pLDDT 83.8 figure or the Foldseek E-value 0.31 on a poster until you have seen them
   yourself.

---

## 1. Environment

```bash
conda env create -f env/environment.yml
conda activate obelisk
source config/config.sh
```

`diamond=2.1.9` is pinned deliberately — `qseq_translated` needs ≥2.0.8 (`C-04`).

**VNom and circuclust are not on conda.** VNom additionally needs USEARCH and MARS in its
`dependencies/` directory; its README has the exact steps, and `scripts/08_run_vnom.sh`
depends on that layout.

```bash
git clone https://github.com/Zheludev/VNom.git
git clone https://github.com/rcedgar/circuclust.git   # or a release binary
```

---

## 2. Pipeline order

Run in this order. Every script is idempotent — re-running skips completed work.

| Step | Command | Notes |
|---|---|---|
| 0 | `bash scripts/00_verify_logan.sh` | Fails loudly if ground truth moved |
| 1 | `bash scripts/01_fetch_references.sh` | Published artifacts (`C-02`); `--from-scratch` for the fallback |
| 2 | `bash scripts/02_build_search_dbs.sh` | **Exits 3 if the decoy control is dirty** (`C-03`) |
| 3 | `bash scripts/03_positive_control.sh all` | **GATE 1 — do not proceed on failure** |
| 4 | `bash scripts/04_select_accessions.sh` | Needs NCBI. Otherwise use `04b` |
| 4b | `bash scripts/04b_select_accessions_offline.sh` | No NCBI. Streams 11 GiB once, then caches |
| — | `python3 scripts/estimate_sweep_cost.py` | **Size the sweep before running it** |
| 5 | `bash scripts/submit_array.sh scripts/05_logan_triage.sbatch work/accessions.txt 100 100` | Pilot first — see below |
| 6 | `python3 scripts/06_collect_hits.py` | **GATE 2** |
| 7 | `bash scripts/submit_array.sh scripts/07_deep_assemble.sbatch work/deep_accessions.txt 20 1` | |
| 8 | `bash scripts/08_run_vnom.sh` | |
| 9 | `python3 scripts/09_circularity_check.py` | |
| 10 | `bash scripts/10_exclude_known.sh` | |
| 11 | `python3 scripts/11_replication_matrix.py` | **GATE 3** |
| — | `seqkit grep -f work/passed_ids.txt results/candidates_novel.fna > results/validated.fna` | |
| 12 | `bash scripts/12_junction_mapping.sh --top 3` | The strongest evidence in the project |
| 13 | `bash scripts/13_rnafold.sh [--rfam]` | |
| 14 | `python3 scripts/14_colabfold_prep.py` | Then ColabFold by hand |
| 15 | `bash scripts/15_foldseek.sh all` | |
| 16 | `bash scripts/16_cluster.sh` | Nomenclature |
| 17 | `python3 scripts/17_wetlab_package.py` | |
| 18 | `bash scripts/18_audit_claims.sh` | Re-verify before presenting |

### Always pilot the sweep first

```bash
head -10 work/accessions.txt | while read -r A; do
  time bash scripts/05_logan_triage.sh "$A"
done
```
Multiply by your accession count to project total runtime. Then a small array. **Never
submit the full array before a pilot comes back clean.**

---

## 3. The three gates

Exit code `3` from any script means a gate failed. That is a **scientific result**, not a
crash — read the output.

| Gate | Script | Passing | Failing |
|---|---|---|---|
| 1 | `03_positive_control.sh` | Recovers Obelisk-S.s at ≥95% identity | **Stop.** Nothing downstream is interpretable, including a negative. |
| 2 | `06_collect_hits.py` | HMM-confirmed hits in ≥3 BioProjects | 1–2 → contamination warning, broaden and re-sweep. 0 → retry relaxed, then report as negative. |
| 3 | `11_replication_matrix.py` | ≥1 candidate across ≥3 BioProjects | 1–2 → tentative, flag explicitly. 0 → publishable negative, **but only because GATE 1 passed**. |

`06` also exits 3 if the hit rate exceeds 20%, which means the threshold is loose or the
`C-03` decoy calibration did not really run.

---

## 4. Switching Logan release

```bash
LOGAN_RELEASE=v1.0 source config/config.sh     # pin the Dec-2023 release
```

Use v1.0 only to reproduce the published sweep as a control. It is deprecated — Logan's own
`Stats-v1.1.md` said `c1.0/` would be deleted about a year after January 2025, and that
window has passed. Default v1.2 covers 38,124,741 accessions against v1.0's 27,269,310.

---

## 5. When something fails

| Symptom | Likely cause |
|---|---|
| `exit 2` | A required input is missing. Run the earlier script the message names. |
| Decoy calibration exits 3 | `C-03`. **Do not tighten the E-value** — that is the plan's advice and it inverts the control. Check the decoys really are shuffled. |
| Collector reports NUCLEOTIDE | `C-04`. The triage worker emitted `qseq` instead of `qseq_translated`. Re-run the sweep. |
| VNom produces no `4_final_clusters` | Normal — nothing nominated, usually at the dual-polarity filter. Confirm the library is **stranded** RNA-seq; VNom cannot work otherwise (`C-15`). |
| `conda activate` fails in a SLURM job | `C-11`. The sbatch scripts source conda's hook; if you wrote your own, do the same. |
| Array submitted with the wrong bound | Use `submit_array.sh`. Never hand-edit the `#SBATCH --array` line (`C-10`). |
| Triage re-downloads everything on re-run | You are not using the `.done` sentinel (`C-05`). |

---

## 6. Before you present

```bash
bash scripts/18_audit_claims.sh
```

Re-checks the Logan release string, the published-artifact URLs, the parquet schema, both
third-party READMEs, and the decoy calibration. It **cannot** check the six literature
figures in `OPEN-C` — those need a human with journal access, and it prints them as a
reminder rather than passing them silently.
