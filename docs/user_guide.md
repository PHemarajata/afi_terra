# AFI Rickettsiales Pipeline — User Guide

**Workflow version:** 0.4.1  
**Platform:** Terra.bio (Cromwell / Google Batch)  
**WDL standard:** 1.0  
**Last updated:** 2026-08

---

## Table of Contents

1. [Overview](#1-overview)
2. [Pipeline Architecture and Data Flow](#2-pipeline-architecture-and-data-flow)
3. [Interpretation Logic and Thresholds](#3-interpretation-logic-and-thresholds)
4. [Prerequisites](#4-prerequisites)
5. [Building the Docker Image](#5-building-the-docker-image)
6. [Setting Up Terra Data Tables](#6-setting-up-terra-data-tables)
7. [Running the Batch Workflow](#7-running-the-batch-workflow)
8. [Running a Single Sample](#8-running-a-single-sample)
9. [Understanding the Outputs](#9-understanding-the-outputs)
10. [Adjusting Detection Thresholds](#10-adjusting-detection-thresholds)
11. [Troubleshooting](#11-troubleshooting)
12. [Reference: Sample Types and Modes](#12-reference-sample-types-and-modes)
13. [Reference: Docker Images](#13-reference-docker-images)

---

## 1. Overview

The `afi_terra` pipeline detects *Orientia*, *Rickettsia*, and a broad panel of other bacterial genera from paired-end 16S metagenomic sequencing reads. It is designed to run on [Terra.bio](https://app.terra.bio) and produces laboratory-grade outputs for routine clinical use alongside validation controls.

### Key design principles

- **Single-pass batch.** NTC (no-template control) background is computed automatically from the NTC samples included in the same submission. No manual two-pass workflow is required.
- **Dual-module interpretation.** *Orientia* and *Rickettsia* are called using 16S minimap2 alignment evidence (higher specificity). All other genera are called from Centrifuge classification counts.
- **Per-run PC8 validity.** Each run is evaluated against its PC_MIX8 positive control. A `run_pc8_valid` flag in the run summary tells you whether the run passed QC (≥ 6 of 8 expected organisms detected).
- **Terra data-model native.** Sample metadata lives in Terra data tables (sample entity + sample_set entity). No JSON file editing is needed once the tables are loaded.

### Two workflows

| WDL file | Use case |
|---|---|
| `wdl/AFI_Rickettsiales_Batch.wdl` | Full sequencing run — multiple samples in one Terra submission. Automatic NTC background computation. Recommended for routine clinical use. |
| `wdl/AFI_Rickettsiales_Main.wdl` | Single sample — standalone processing. Requires a pre-computed NTC background file. |

---

## 2. Pipeline Architecture and Data Flow

Every sample passes through the same processing steps. In the batch workflow, Phase 1 and Phase 2 run as parallel scatters; a single gather task between them computes the NTC background automatically.

```
Phase 1 — per sample, in parallel
───────────────────────────────────────────────────────────────────────
  [optional] NCBI SRA Human Scrubber    remove human reads
       │
       ▼
  fastp                                 adapter trimming, quality filtering
       │
       ├────────────────────────────────┐
       ▼                                ▼
  Centrifuge + kreport parse       minimap2 → samtools sort/index
  genus counts TSV                 align to 16S Rickettsiales panel
       │                                │
       │                                ▼
       │                      ExtractMetrics (Python)
       │                      per-genus reads + breadth TSV
       │                                │
       └──── NTC/NC samples expose their output files ────┘

Gather — once per batch
───────────────────────────────────────────────────────────────────────
  BuildNTCBackground            aggregate NTC align + centrifuge metrics
                                → ntc_background.tsv

Phase 2 — per sample, in parallel
───────────────────────────────────────────────────────────────────────
  InterpretCalls (call_taxa.py)     NTC-aware dual-module calling
                                    → calls.tsv
       │
       ├── validation mode: CompareExpectedConcordance → validation_summary.tsv
       └── routine mode:    SummarizeRoutineTaxa       → routine_summary.tsv

Final gather — once per batch
───────────────────────────────────────────────────────────────────────
  BuildRunSummary               merge all summaries → run_summary.tsv
                                annotate with run_id and run_pc8_valid
```

### Step-by-step details

| Step | Tool | Docker image | Output |
|---|---|---|---|
| Human scrubbing (optional) | NCBI SRA Scrubber | `sra-human-scrubber:2.2.1` | Dehosted R1/R2 fastq |
| QC / adapter trimming | fastp | `staphb/fastp:0.23.4` | Clean R1/R2 fastq |
| Taxonomic classification | Centrifuge | `phemarajata614/centrifuger:1.1.0` | kreport TSV |
| Kreport parsing | Python script | `afi-terra:0.4.1` | `genus_counts.tsv` |
| 16S alignment | minimap2 + samtools | `afi-terra:0.4.1` | `align.bam`, `align.bam.bai` |
| Alignment metrics | Python script | `afi-terra:0.4.1` | `align_metrics.tsv` |
| NTC background (batch only) | Python script | `afi-terra:0.4.1` | `ntc_background.tsv` |
| Interpretation | Python script | `afi-terra:0.4.1` | `calls.tsv` |
| Validation summary | Inline Python | `afi-terra:0.4.1` | `validation_summary.tsv` |
| Routine summary | Inline Python | `afi-terra:0.4.1` | `routine_summary.tsv` |
| Run summary (batch only) | Inline Python | `afi-terra:0.4.1` | `run_summary.tsv` |

### How NTC background is computed automatically

The batch workflow uses a WDL 1.0-compatible pattern to identify NTC/NC samples without requiring a separate submission:

1. Inside the Phase 1 scatter, each sample checks its own `sample_type`.
2. If `sample_type` is `NTC` or `NC`, the sample's `align_metrics.tsv` and `genus_counts.tsv` are exposed as optional (`File?`) output variables.
3. After the scatter, `select_all()` collects only the non-null values — i.e., only the NTC/NC files.
4. `BuildNTCBackground` receives these pre-filtered file lists directly and computes per-genus maximum read counts.

This is why a batch with no NTC/NC sample will produce an empty background file (all thresholds default to zero). **Always include at least one NTC sample per batch.**

---

## 3. Interpretation Logic and Thresholds

### Module 3 — Alignment (Orientia and Rickettsia only)

*Orientia* and *Rickettsia* are evaluated exclusively from 16S minimap2 alignment evidence. Four call tiers are evaluated in order; the first matching tier is assigned:

| Call | Criteria |
|---|---|
| **Confirmed** | `mapped_reads ≥ 100` AND `max_breadth ≥ 0.25` AND `reads ≥ 5× NTC alignment reads` |
| **Probable** | `reads ≥ 50` AND `max_breadth ≥ 0.20` AND `reads > NTC alignment reads` |
| **Not_Confirmed** | `reads ≥ 50` AND `reads ≤ NTC alignment reads` |
| **Negative** | `reads < 50` |

Centrifuge counts for *Orientia* and *Rickettsia* are collected but **not used for calling** — alignment is the authoritative source for these two genera.

### Module 1 — Centrifuge (all other genera)

Every genus other than *Orientia* and *Rickettsia* is called from Centrifuge kreport genus-level clade read counts. Two tiers:

| Call | Criteria |
|---|---|
| **Detected** | `reads ≥ 500` AND `reads ≥ 5× NTC centrifuge reads` |
| **Not_Detected** | otherwise |

Clade reads are used (column 1 of the kreport, rank `G`) so that reads assigned to species within a genus are counted toward the genus total.

### Positive call summary

The following calls are considered "positive" and drive concordance checks and taxa lists:

| Source | Positive calls |
|---|---|
| Alignment (Orientia/Rickettsia) | Confirmed, Probable |
| Centrifuge (all other genera) | Detected |

`Not_Confirmed` is **not** a positive call. It indicates signal above the minimum read threshold but indistinguishable from NTC noise.

### NTC background

For each genus, the NTC background value is the **maximum** read count across all NTC/NC samples in the batch, separately for alignment (`align_ntc_reads`) and Centrifuge (`cfr_ntc_reads`) sources. Taking the maximum is conservative — it uses the worst-case background level observed in the run.

Genera absent from the NTC background default to 0, which effectively disables the fold-change gate for those genera.

### PC8 validity

For any `PC_MIX8` sample processed in validation mode:

| Value | Condition |
|---|---|
| `pc8_pass = true` | ≥ 6 of the expected taxa received a positive call |
| `pc8_pass = false` | < 6 of the expected taxa received a positive call |
| `pc8_pass = not_applicable` | Sample is not PC_MIX8, or no expected taxa were provided |

The `run_pc8_valid` column in `run_summary.tsv` carries this value across every row in the run, making it easy to flag an entire run as QC-failed at a glance.

**Recommended `expected_taxa` for Zymo PC_MIX8:**
```
Orientia;Rickettsia;Escherichia;Salmonella;Staphylococcus;Pseudomonas;Listeria;Enterococcus
```

---

## 4. Prerequisites

### Terra workspace

- An active Terra workspace linked to a Google Cloud project.
- The workspace must have the Google Batch API enabled (the default for Terra workspaces).

### GCS data

Upload these files to your workspace bucket before running:

| Item | Notes |
|---|---|
| R1 and R2 FASTQ files | Gzipped FASTQ (`.fastq.gz`) recommended. One pair per sample. |
| 16S Rickettsiales panel | FASTA file. Stored in the repo under `align_rickettsiales_16S/refs/rickettsiales_panel_16S.clean.fa`. Pre-built `.mmi` index also available. |
| Centrifuge database | One or more `.tar.gz` archives containing the combined bacteria + archaea + Rickettsiales Centrifuge index. |

### Docker image

You must build and push `phemarajata614/afi-terra:0.4.1` (or your own tagged version) before running. The other three images used by the pipeline (`staphb/fastp`, `phemarajata614/centrifuger`, `sra-human-scrubber`) are public and do not require a build step.

### Local tools (for Docker build only)

- Docker Engine or Docker Desktop with linux/amd64 build support (use `--platform linux/amd64` on Apple Silicon).
- Docker Hub account with push access to your image repository.
- A local clone of this repository.

---

## 5. Building the Docker Image

The `afi-terra` image contains the Python analysis scripts, samtools, minimap2, and pandas. Rebuild and push a new image any time you change a script in `scripts/` or the `Dockerfile`.

```bash
# From the repo root on your local machine
docker login
bash scripts/build_push_afi_core_image.sh phemarajata614 0.4.1 linux/amd64
```

The script:
1. Builds the image for `linux/amd64` (required — Terra runs on GCP VMs).
2. Runs a smoke test verifying samtools, minimap2, Python/pandas, and all scripts respond to `--help`.
3. Pushes `phemarajata614/afi-terra:0.4.1` to Docker Hub.

> **Tag discipline:** Terra caches Docker images. Pushing a content change to the same tag may leave Terra using the stale cached version. Always increment the tag (e.g., `:0.4.2`) when pushing a change, then update `afi_core_docker` in your Terra input config to match.

### Contents of the afi-terra image

| Component | Purpose |
|---|---|
| `python=3.11` + `pandas` | All Python scripts |
| `samtools` | BAM metrics extraction (`samtools idxstats`, `samtools depth`) |
| `minimap2` | 16S alignment (MinimapRick16S task) |
| Scripts under `/opt/afi/scripts/` | Pipeline logic |

Scripts bundled in the image:

| Script | Invoked by task |
|---|---|
| `parse_centrifuge_kreport.py` | ParseCentrifugerKreport |
| `extract_rick16s_metrics.py` | ExtractMetrics |
| `build_ntc_background.py` | BuildNTCBackground |
| `call_taxa.py` | InterpretCalls |

### Equivalent manual build commands

```bash
docker build --platform linux/amd64 -t phemarajata614/afi-terra:0.4.1 .
docker push phemarajata614/afi-terra:0.4.1
```

---

## 6. Setting Up Terra Data Tables

The batch workflow reads sample metadata from Terra's entity tables. Three table uploads are required before launching a run. Template files are in `wdl/inputs/`.

### 6.1 Sample entity table — `terra_sample_table.tsv`

**Upload:** Data tab → Import Data → Upload TSV

Each row is one sample. Column descriptions:

| Column | Required | Description |
|---|---|---|
| `entity:sample_id` | Yes | Unique sample identifier. Becomes the entity name in Terra. |
| `r1_fastq` | Yes | `gs://` path to R1 FASTQ file. |
| `r2_fastq` | Yes | `gs://` path to R2 FASTQ file. |
| `sample_type` | Yes | See Section 12 for allowed values. |
| `mode` | Yes | `routine` or `validation`. |
| `expected_taxa` | Required for validation samples; blank otherwise | Semicolon-delimited genera expected to be detected. Leave the cell empty for NTC and routine clinical samples. |

**Example rows:**

```
entity:sample_id  r1_fastq                               r2_fastq                               sample_type  mode        expected_taxa
NTC_S1            gs://fc-.../NTC_S1_R1.fastq.gz         gs://fc-.../NTC_S1_R2.fastq.gz         NTC          routine
PC_S12            gs://fc-.../PC_S12_R1.fastq.gz         gs://fc-.../PC_S12_R2.fastq.gz         PC_MIX8      validation  Orientia;Rickettsia;Escherichia;Salmonella;Staphylococcus;Pseudomonas;Listeria;Enterococcus
CLINICAL_S2       gs://fc-.../CLINICAL_S2_R1.fastq.gz    gs://fc-.../CLINICAL_S2_R2.fastq.gz    clinical     routine
```

### 6.2 Sample_set entity table — `terra_sample_set_table.tsv`

**Upload:** Data tab → Import Data → Upload TSV

One row per sequencing run. Stores the `run_id` attribute that is written into every row of `run_summary.tsv`.

```
entity:sample_set_id   run_id
RUN_2026_04_03_001     RUN_2026_04_03_001
```

`run_id` can be any string — use your laboratory run identifier or a date-based ID.

> **One submission = one run.** The batch workflow computes a single NTC background from all NTC/NC samples in the submitted sample_set. If you need to process samples from multiple distinct sequencing runs, either submit them separately or use a naming convention that distinguishes the runs within your QC documentation.

### 6.3 Membership table — `terra_sample_set_membership.tsv`

**Upload:** Data tab → Import Data → Upload TSV

Links each sample to its run. One row per sample.

```
membership:sample_set_id   sample
RUN_2026_04_03_001         NTC_S1
RUN_2026_04_03_001         PC_S12
RUN_2026_04_03_001         CLINICAL_S2
RUN_2026_04_03_001         CLINICAL_S3
```

### 6.4 Workspace attribute — `terra_workspace_attributes.tsv`

**Upload:** Workspace tab → Edit → Import

Stores the `use_human_scrub` flag that applies to every sample in every run from this workspace.

```
attribute          value
use_human_scrub    true
```

Set to `false` if samples have already been dehosted upstream or if you are processing non-clinical data.

---

## 7. Running the Batch Workflow

### 7.1 Add the workflow to your Terra workspace

1. Navigate to **Workflows → Find a Workflow**.
2. Import `AFI_Rickettsiales_Batch` from your method repository (GitHub/Dockstore) or upload the WDL directly.
3. When uploading via ZIP, the archive must preserve the relative directory structure including all task WDLs and the NCBI scrub WDL:

```
wdl/AFI_Rickettsiales_Batch.wdl
wdl/AFI_Rickettsiales_Main.wdl
wdl/tasks/align.wdl
wdl/tasks/classify.wdl
wdl/tasks/interpret.wdl
wdl/tasks/metrics.wdl
wdl/tasks/preprocess.wdl
wdl/tasks/validate.wdl
NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl
```

### 7.2 Configure workflow inputs

Open the workflow → **Inputs** tab. Map each input to a Terra data table expression or a literal value:

| Workflow input | Terra mapping | Notes |
|---|---|---|
| `run_id` | `this.run_id` | From the sample_set entity |
| `sample_ids` | `this.samples.sample_id` | |
| `r1_fastqs` | `this.samples.r1_fastq` | |
| `r2_fastqs` | `this.samples.r2_fastq` | |
| `sample_types` | `this.samples.sample_type` | |
| `modes` | `this.samples.mode` | |
| `expected_taxa` | `this.samples.expected_taxa` | Empty string for non-validation samples |
| `use_human_scrub` | `workspace.use_human_scrub` | Or hardcode `true` / `false` |
| `rickettsiales_panel` | Literal GCS path | e.g., `"gs://fc-.../rickettsiales_panel_16S.clean.fa"` |
| `centrifuger_db` | Literal string | Index prefix name (no path or extension) |
| `centrifuger_db_archives` | Literal JSON array | `["gs://fc-.../centrifuger_index.tar.gz"]` |
| `afi_core_docker` | Literal | `"phemarajata614/afi-terra:0.4.1"` |

The full input config template is at `wdl/inputs/terra_workflow_input_config.json`.

### 7.3 Set compute resources for Centrifuge (important)

The Centrifuge classification step loads a large database (~67 GB compressed, ~90–110 GB in memory). If Terra allocates the default small VM, the task will segfault silently or produce an empty kreport.

Set these in your workflow inputs:

| Input | Recommended value | Reason |
|---|---|---|
| `centrifuger_memory` | `"128G"` | Combined bacteria + archaea + Rickettsiales database |
| `centrifuger_disks` | `"local-disk 500 HDD"` | Extracted index is much larger than the compressed archive |
| `classify_threads` | `16` | Matches the `cpu` request in the task runtime block |

### 7.4 Launch

1. Navigate to **Workflows** and select `AFI_Rickettsiales_Batch`.
2. Choose **Run workflow(s) with inputs defined by data table**.
3. Select **sample_set** as the entity type.
4. Select the row for your run (e.g., `RUN_2026_04_03_001`).
5. Click **Run Analysis**.

### 7.5 Required sample composition per batch

Every batch submission must include:

| Control | Minimum | Reason |
|---|---|---|
| NTC or NC sample | 1 | Without it, NTC background is all zeros — fold-change thresholds are disabled |
| PC_MIX8 in validation mode | Strongly recommended | Required to compute `run_pc8_valid` |

Without a PC_MIX8 sample, `run_pc8_valid` will read `no_pc8_in_run` in the run summary.

---

## 8. Running a Single Sample

Use `AFI_Rickettsiales_Main` when you need to process one sample in isolation — for example, to re-interpret a sample with a different NTC background, or for development and debugging.

> The single-sample workflow requires a pre-computed `ntc_background.tsv` file. Use the output from a previous batch run, or supply a placeholder file (all zeros) for a first-pass look without NTC correction.

### Creating a placeholder NTC background

Create a file named `ntc_background.placeholder.tsv` with this content:

```
genus	align_ntc_reads	cfr_ntc_reads
Orientia	0	0
Rickettsia	0	0
```

Upload it to GCS. Use its `gs://` path as the `ntc_background` workflow input.

### Key inputs

| Input | Type | Default | Description |
|---|---|---|---|
| `sample_id` | String | — | Unique sample name |
| `sample_type` | String | `"clinical"` | See Section 12 |
| `mode` | String | `"routine"` | `routine` or `validation` |
| `use_human_scrub` | Boolean | `true` | Enable NCBI SRA human read scrubbing |
| `r1_fastq` | File | — | GCS path to R1 FASTQ |
| `r2_fastq` | File | — | GCS path to R2 FASTQ |
| `rickettsiales_panel` | File | — | GCS path to 16S reference FASTA |
| `ntc_background` | File | — | Pre-computed or placeholder NTC background TSV |
| `centrifuger_db` | String | `""` | Centrifuge index prefix |
| `centrifuger_db_archives` | Array[File] | `[]` | GCS paths to Centrifuge archive files |
| `expected_taxon` | String | `""` | Semicolon-delimited expected genera (validation mode) |
| `centrifuger_memory` | String | `"128G"` | Memory for Centrifuge task |
| `centrifuger_disks` | String | `"local-disk 500 HDD"` | Disk for Centrifuge task |
| `classify_threads` | Int | `16` | CPU threads for Centrifuge |
| `afi_core_docker` | String | `"phemarajata614/afi-terra:0.4.1"` | Core image |
| `fastp_docker` | String | `"staphb/fastp:0.23.4"` | QC trimming image |
| `centrifuger_docker` | String | `"phemarajata614/centrifuger:1.1.0"` | Classifier image |

### Example input JSON (routine clinical sample)

```json
{
  "AFI_Rickettsiales_Main.sample_id": "CLINICAL_001",
  "AFI_Rickettsiales_Main.sample_type": "clinical",
  "AFI_Rickettsiales_Main.mode": "routine",
  "AFI_Rickettsiales_Main.use_human_scrub": true,
  "AFI_Rickettsiales_Main.r1_fastq": "gs://YOUR_BUCKET/fastq/CLINICAL_001_R1.fastq.gz",
  "AFI_Rickettsiales_Main.r2_fastq": "gs://YOUR_BUCKET/fastq/CLINICAL_001_R2.fastq.gz",
  "AFI_Rickettsiales_Main.rickettsiales_panel": "gs://YOUR_BUCKET/ref/rickettsiales_panel_16S.clean.fa",
  "AFI_Rickettsiales_Main.ntc_background": "gs://YOUR_BUCKET/ref/ntc_background.placeholder.tsv",
  "AFI_Rickettsiales_Main.centrifuger_db": "centrifuger_bact_arch_plus_rickettsiales",
  "AFI_Rickettsiales_Main.centrifuger_db_archives": [
    "gs://YOUR_BUCKET/db/centrifuger_index.tar.gz"
  ],
  "AFI_Rickettsiales_Main.afi_core_docker": "phemarajata614/afi-terra:0.4.1",
  "AFI_Rickettsiales_Main.centrifuger_docker": "phemarajata614/centrifuger:1.1.0"
}
```

### Example input JSON (validation — PC_MIX8)

```json
{
  "AFI_Rickettsiales_Main.sample_id": "PC_S12",
  "AFI_Rickettsiales_Main.sample_type": "PC_MIX8",
  "AFI_Rickettsiales_Main.mode": "validation",
  "AFI_Rickettsiales_Main.expected_taxon": "Orientia;Rickettsia;Escherichia;Salmonella;Staphylococcus;Pseudomonas;Listeria;Enterococcus",
  "AFI_Rickettsiales_Main.use_human_scrub": true,
  "AFI_Rickettsiales_Main.r1_fastq": "gs://YOUR_BUCKET/fastq/PC_S12_R1.fastq.gz",
  "AFI_Rickettsiales_Main.r2_fastq": "gs://YOUR_BUCKET/fastq/PC_S12_R2.fastq.gz",
  "AFI_Rickettsiales_Main.rickettsiales_panel": "gs://YOUR_BUCKET/ref/rickettsiales_panel_16S.clean.fa",
  "AFI_Rickettsiales_Main.ntc_background": "gs://YOUR_BUCKET/ref/run_ntc_background.tsv",
  "AFI_Rickettsiales_Main.centrifuger_db": "centrifuger_bact_arch_plus_rickettsiales",
  "AFI_Rickettsiales_Main.centrifuger_db_archives": [
    "gs://YOUR_BUCKET/db/centrifuger_index.tar.gz"
  ],
  "AFI_Rickettsiales_Main.afi_core_docker": "phemarajata614/afi-terra:0.4.1",
  "AFI_Rickettsiales_Main.centrifuger_docker": "phemarajata614/centrifuger:1.1.0"
}
```

---

## 9. Understanding the Outputs

### 9.1 Per-sample outputs

#### `calls.tsv` — taxa calls (primary result per sample)

One row per genus evaluated. Both *Orientia* and *Rickettsia* always appear (even if zero reads), from the alignment source. All other genera appear only if detected by Centrifuge.

| Column | Description |
|---|---|
| `sample` | Sample identifier |
| `genus` | Genus name |
| `source` | `alignment` (Orientia/Rickettsia) or `centrifuge` (all others) |
| `reads` | Mapped reads (alignment) or clade reads (centrifuge) |
| `breadth` | Fraction of reference covered (alignment only; blank for centrifuge) |
| `ntc_reads` | NTC background reads for this genus and source |
| `call` | **Confirmed / Probable / Not_Confirmed / Negative** (alignment) or **Detected / Not_Detected** (centrifuge) |

#### `align_metrics.tsv` — raw alignment metrics

| Column | Description |
|---|---|
| `genus` | `Orientia` or `Rickettsia` |
| `mapped_reads` | Reads mapping to this genus's 16S references |
| `max_breadth` | Maximum breadth of coverage across this genus's references (0.0–1.0) |

#### `genus_counts.tsv` — centrifuge genus-level counts

| Column | Description |
|---|---|
| `genus` | Genus name |
| `reads` | Clade read count from kreport (includes reads to child species) |

#### `validation_summary.tsv` — concordance check (validation mode only)

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `sample_type` | e.g., `PC_MIX8` |
| `expected_taxa` | Comma-separated expected genera |
| `detected_expected_taxa` | Expected genera that were detected |
| `missing_expected_taxa` | Expected genera that were not detected |
| `unexpected_detected_taxa` | Detected genera not in the expected list |
| `validation_result` | `Concordant` (all expected detected) or `Discordant` |
| `pc8_pass` | `true`, `false`, or `not_applicable` |

#### `routine_summary.tsv` — detected taxa list (routine mode only)

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `sample_type` | e.g., `clinical` |
| `taxa_present` | Comma-separated genera with a positive call |
| `n_taxa_present` | Count of detected genera |
| `routine_positive_control` | `true` if sample_type is PC_MIX8 or PC |

### 9.2 Run-level outputs (batch workflow only)

#### `ntc_background.tsv` — auto-computed NTC thresholds

| Column | Description |
|---|---|
| `genus` | Genus name |
| `align_ntc_reads` | Maximum alignment reads across all NTC/NC samples in this run |
| `cfr_ntc_reads` | Maximum centrifuge reads across all NTC/NC samples in this run |

*Orientia* and *Rickettsia* always appear in this file even if absent from the NTC kreport (with zero values).

This file can be downloaded from Terra and used as the `ntc_background` input for future single-sample reruns from the same run context.

#### `run_summary.tsv` — final run report

One row per sample. This is the primary deliverable for a run.

| Column | Description |
|---|---|
| `run_id` | Run identifier from the sample_set entity |
| `sample_id` | Sample identifier |
| `sample_type` | e.g., `clinical`, `NTC`, `PC_MIX8` |
| `detected_taxa` | Comma-separated genera with a positive call |
| `n_detected` | Count of detected genera |
| `validation_result` | `Concordant`, `Discordant`, or blank for routine samples |
| `pc8_pass` | `true`, `false`, `not_applicable`, or blank |
| `run_pc8_valid` | Run-level PC8 flag: `true`, `false`, or `no_pc8_in_run` |

A `run_pc8_valid = false` means the PC_MIX8 control failed for this run. All clinical results should be interpreted with caution.

### 9.3 Intermediate outputs retained in Terra

Terra stores all task outputs in your workspace GCS bucket under `submissions/<submission-id>/`. The following intermediate files are exposed as workflow outputs for audit purposes:

- Dehosted FASTQ files (if `use_human_scrub = true`)
- QC-trimmed FASTQ files (post-fastp)
- Centrifuge kreport TSVs
- minimap2 BAM + BAI files

---

## 10. Adjusting Detection Thresholds

All thresholds are workflow-level inputs with sensible defaults. Override them in the Terra inputs panel or your JSON config without touching any code.

### Alignment thresholds (Orientia / Rickettsia)

| Input | Default | Effect |
|---|---|---|
| `align_confirm_reads` | `100` | Minimum mapped reads for a **Confirmed** call |
| `align_confirm_breadth` | `0.25` | Minimum breadth of coverage for a **Confirmed** call |
| `align_fold` | `5.0` | Minimum fold above NTC alignment reads for a **Confirmed** call |

Probable threshold (not configurable): reads ≥ 50 AND breadth ≥ 0.20 AND reads > NTC reads.

### Centrifuge thresholds (all other genera)

| Input | Default | Effect |
|---|---|---|
| `cfr_floor` | `500` | Minimum clade reads for a **Detected** call |
| `cfr_fold` | `5.0` | Minimum fold above NTC centrifuge reads for a **Detected** call |

### Guidance for threshold adjustment

| Situation | Suggested adjustment |
|---|---|
| High NTC background for a genus | The NTC mechanism handles this automatically — no adjustment needed |
| Very low-abundance clinical samples | Lower `align_confirm_reads` to 50 and `cfr_floor` to 200–250 |
| Higher specificity required | Raise `align_confirm_reads` to 200 or `align_confirm_breadth` to 0.30 |
| Frequent Probable calls believed to be real | First check NTC kreport for that genus; if NTC is clean, lower `align_fold` to 3.0 |

---

## 11. Troubleshooting

### Centrifuge segfaults or produces an empty kreport

**Cause:** Terra allocated a VM with insufficient memory for the database.

**Fix:** Set `centrifuger_memory = "128G"` and `centrifuger_disks = "local-disk 500 HDD"` in your workflow inputs. Without these, Terra may launch the task on a 1 CPU / 2 GB VM, which is far too small.

---

### `build_ntc_background.py: error: the following arguments are required: --sample-types-file`

**Cause:** Terra is running an old Docker image (`afi-terra:0.3.x`) that has a different script signature. The `--sample-types-file` argument was removed in `0.4.0`.

**Fix:** Build and push the new image (Section 5), then set `afi_core_docker = "phemarajata614/afi-terra:0.4.1"` in your workflow inputs. Always increment the tag to force Terra to pull the new image rather than using a cached version.

---

### `fastp: command not found` in the FastpClean task

**Cause:** The `fastp_docker` input is pointing to an image that does not contain fastp. As of `0.4.0`, fastp was removed from the `afi-terra` image and delegated to `staphb/fastp`.

**Fix:** Ensure `fastp_docker = "staphb/fastp:0.23.4"` is set in your workflow inputs (it is the task's default, but may have been overridden).

---

### `ln: failed to create symbolic link '/cromwell_root': Permission denied`

**Cause:** Standard Cromwell Batch localization behavior. Cromwell tries to create a symlink in `/cromwell_root` and fails gracefully when the directory is not writable (normal in the micromamba base image).

**Action:** None. This message is safe to ignore.

---

### `BuildNTCBackground` receives zero files

**Cause:** No sample in the batch has `sample_type = NTC` or `sample_type = NC`. The `select_all()` call returns an empty array.

**Fix:** Include at least one NTC or NC sample in every batch submission. Verify the `sample_type` column uses exactly `NTC` or `NC` (case-sensitive, no trailing spaces).

---

### `run_pc8_valid = no_pc8_in_run` in run summary

**Cause:** No `PC_MIX8` sample was included in the batch, or the `sample_type` was not exactly `PC_MIX8`.

**Action:** If a PC_MIX8 was intended, check the Terra data table for the correct spelling. If the run legitimately had no PC_MIX8, this value is expected.

---

### All calls show `Negative` or `Not_Detected` despite expected signal

**Cause 1:** The Centrifuge task failed or produced an empty kreport — check the `centrifuger_kreports` output for content.

**Cause 2:** The FASTQ files are empty or the scrubber removed most reads — check `clean_r1`/`clean_r2` file sizes.

**Cause 3:** The database does not contain the target organism's sequences.

**Diagnostic approach:** Download `align_metrics.tsv` and `genus_counts.tsv` from Terra. If raw read counts are present but calls are negative, the threshold or NTC background is the issue. If raw counts are zero, the problem is upstream (alignment/classification failure).

---

### Validation result is `Discordant` for an expected organism

**Likely causes:**

1. The organism was detected but fell below the confirmation threshold. Check `calls.tsv` for the genus — it may have a `Probable` or `Not_Confirmed` call rather than `Confirmed`. `Probable` is a positive call and should produce `Concordant`. `Not_Confirmed` is not positive.
2. The genus name in `expected_taxa` does not exactly match the genus name in `calls.tsv`. Names are case-sensitive.
3. NTC background for that genus is high. Check `ntc_background.tsv` for the genus's `align_ntc_reads` value.

---

### WDL parse error in Terra

Terra uses **WDL 1.0** (Cromwell). If you have edited the WDL, avoid these WDL 2.0 patterns:

| Pattern | Issue | WDL 1.0 fix |
|---|---|---|
| `samples[*].field` | Splat operator not supported | Use `range(length(arr))` and index access |
| Referencing scatter variable outside scatter block | Not in scope | Use conditional declarations (`if (cond) { File f = task.out }`) inside the scatter and `select_all()` outside |
| `Optional[T]` used directly without unwrapping | Type mismatch | Wrap with `select_first([val, default])` or `select_all()` |

---

### Centrifuge database prefix not found after extraction

**Cause:** The `centrifuger_db` prefix string does not match the actual filenames inside the archive.

**Fix:** Extract the archive locally and check what files are present:

```bash
tar -tzf centrifuger_index.tar.gz | grep ".1.cfr\|.1.cf"
```

The value of `centrifuger_db` must exactly match the filename prefix (without the `.1.cfr` or `.1.cf` suffix and without any path). Example: if the archive contains `centrifuger_bact_arch_plus_rickettsiales.1.cfr`, set:

```json
"AFI_Rickettsiales_Batch.centrifuger_db": "centrifuger_bact_arch_plus_rickettsiales"
```

---

## 12. Reference: Sample Types and Modes

### Sample types

| `sample_type` | Description | Contributes to NTC background | PC8 validity evaluated |
|---|---|---|---|
| `NTC` | No-template control | Yes | No |
| `NC` | Negative control (alias for NTC) | Yes | No |
| `PC_MIX8` | 8-organism positive control (Zymo) | No | Yes — must be `mode = validation` |
| `PC_SINGLE` | Single-organism positive control | No | No |
| `MIXED4` | 4-organism mixed positive control | No | No |
| `PC` | Generic positive control (routine monitoring) | No | No |
| `clinical` | Patient specimen | No | No |

### Modes

| `mode` | Output file | What it checks |
|---|---|---|
| `routine` | `routine_summary.tsv` | Lists detected taxa; flags routine positive controls |
| `validation` | `validation_summary.tsv` | Compares detected taxa against `expected_taxa`; computes `pc8_pass` for PC_MIX8 |

A single batch run can contain a mix of `routine` and `validation` samples.

### sample_type / mode combinations in a typical run

| Sample | sample_type | mode | expected_taxa |
|---|---|---|---|
| Negative control | `NTC` | `routine` | (blank) |
| 8-plex positive control | `PC_MIX8` | `validation` | `Orientia;Rickettsia;Escherichia;...` |
| Single-target positive control | `PC_SINGLE` | `validation` | `Orientia` |
| Clinical specimen | `clinical` | `routine` | (blank) |

---

## 13. Reference: Docker Images

| Image | Version | Used by task(s) | Who maintains it |
|---|---|---|---|
| `phemarajata614/afi-terra` | `0.4.1` | ParseCentrifugerKreport, ExtractMetrics, MinimapRick16S, InterpretCalls, CompareExpectedConcordance, SummarizeRoutineTaxa, BuildNTCBackground, BuildRunSummary | This repo — rebuild when scripts or Dockerfile change |
| `staphb/fastp` | `0.23.4` | FastpClean | StaPH-B (public) |
| `phemarajata614/centrifuger` | `1.1.0` | RunCentrifuger | This repo — rebuild if Centrifuge version changes |
| `us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber` | `2.2.1` | ncbi_scrub_pe | Theiagen / NCBI (public) |

### What's in each image

**`afi-terra:0.4.1`** (built from `Dockerfile` in this repo):
- Python 3.11 + pandas
- samtools (for BAM metrics)
- minimap2 (for 16S alignment)
- All scripts in `/opt/afi/scripts/`

**`staphb/fastp:0.23.4`** (public, StaPH-B):
- fastp only — lean single-tool container

**`phemarajata614/centrifuger:1.1.0`** (built separately):
- Centrifuge/Centrifuger classifier and kreport tools

**`sra-human-scrubber:2.2.1`** (public, Theiagen/NCBI):
- NCBI SRA human read scrubbing tool

### Overriding image versions

All image inputs except the NCBI scrubber are exposed as workflow-level inputs. To use a different version, set the appropriate input in your Terra configuration:

```json
{
  "AFI_Rickettsiales_Batch.afi_core_docker":    "phemarajata614/afi-terra:0.4.1",
  "AFI_Rickettsiales_Batch.fastp_docker":       "staphb/fastp:0.23.4",
  "AFI_Rickettsiales_Batch.minimap_docker":     "phemarajata614/afi-terra:0.4.1",
  "AFI_Rickettsiales_Batch.centrifuger_docker": "phemarajata614/centrifuger:1.1.0"
}
```

> **Note:** `minimap_docker` defaults to `afi-terra:0.4.1` because minimap2 is bundled in the core image. If you wish to use a separate minimap2 container, set `minimap_docker` to `staphb/minimap2:2.28` and ensure that image includes samtools (needed for the sort and index steps in the same task).

### Rebuilding the core image

```bash
# Increment the tag any time you change scripts/ or Dockerfile
bash scripts/build_push_afi_core_image.sh phemarajata614 0.4.2 linux/amd64

# Then update afi_core_docker in Terra:
# "AFI_Rickettsiales_Batch.afi_core_docker": "phemarajata614/afi-terra:0.4.2"
```
