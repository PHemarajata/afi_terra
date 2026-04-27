# AFI Terra Pipeline: User Guide

**Version:** 0.4.1  
**Workflows:** `AFI_16S_Main` / `AFI_16S_Batch`  
**Platform:** Terra (Broad Institute)  
**Last updated:** 2026-04-06

---

## Table of Contents

1. [Overview](#1-overview)
2. [Architecture and Data Flow](#2-architecture-and-data-flow)
3. [Reference Files and Databases](#3-reference-files-and-databases)
4. [Importing Workflows into Terra](#4-importing-workflows-into-terra)
5. [Running the Single-Sample Workflow](#5-running-the-single-sample-workflow)
6. [Running the Batch Workflow](#6-running-the-batch-workflow)
7. [The Two-Pass Pattern](#7-the-two-pass-pattern)
8. [Terra Sheet Builder (GUI)](#8-terra-sheet-builder-gui)
9. [Helper Scripts](#9-helper-scripts)
10. [Output Files Reference](#10-output-files-reference)
11. [Interpreting Calls and Summaries](#11-interpreting-calls-and-summaries)
12. [Detection Thresholds](#12-detection-thresholds)
13. [Docker Images](#13-docker-images)
14. [Troubleshooting](#14-troubleshooting)

---

## 1. Overview

The `afi_terra` pipeline detects and identifies Rickettsiales organisms — primarily *Orientia*, *Rickettsia*, and related genera — from 16S rRNA metagenomic sequencing data. It runs on Terra, the Broad Institute's cloud execution platform, using WDL 1.0 workflows backed by Docker images.

### Target organisms

The pipeline is designed to detect genera in the order Rickettsiales. For *Orientia* and *Rickettsia* specifically, a confirmatory 16S rRNA alignment step provides higher specificity. All other genera detected by the classifier (e.g., *Leptospira*, *Burkholderia*, *Anaplasma*) are reported from Centrifuger read counts alone.

### Two workflows

| Workflow | WDL file | Use case |
|---|---|---|
| `AFI_16S_Main` | `wdl/AFI_16S_Main.wdl` | Single sample |
| `AFI_16S_Batch` | `wdl/AFI_16S_Batch.wdl` | Full sequencing run (scatter/gather) |

For clinical use, run the batch workflow. It processes all samples in a run together, automatically computing a per-run NTC (no-template control) background from the NTC/NC samples in that run, and producing a consolidated `run_summary.tsv`. The **Terra Sheet Builder** desktop app (Section 8) can generate the required sample sheet TSV without manual editing.

### Two sample modes

| Mode | Use case | Primary output |
|---|---|---|
| `routine` | Clinical samples with unknown organisms | `routine_summary.tsv` |
| `validation` | Samples with known expected organisms | `validation_summary.tsv` |

A single batch run can contain both modes simultaneously.

### Classifier

The pipeline uses **Centrifuger** as its sole classifier, backed by a combined bacteria/archaea + Rickettsiales database. Centrifuge read counts provide evidence for all genera; *Orientia* and *Rickettsia* additionally receive a confirmatory 16S rRNA alignment step for higher specificity.

---

## 2. Architecture and Data Flow

### Processing steps

Each sample passes through up to nine steps. Steps 1–5 run in Phase 1 of the batch scatter; steps 6–8 run in Phase 2 (after per-run NTC backgrounds are computed and matched); step 9 is the final batch-level gather.

```
Step 1  Human read dehosting        NCBI SRA Scrubber (optional, default on)
Step 2  QC / adapter trimming       fastp v0.23.4
Step 3  Taxonomic classification    Centrifuger -> kreport -> genus counts TSV
Step 4  16S alignment               minimap2 -> BAM
Step 5  Alignment metrics           extract_rick16s_metrics.py -> align_metrics.tsv
        ---- batch gather: BuildNTCBackground (NTC/NC samples only; per run_id) ----
        ---- batch gather: MatchNTCBackground (assigns each sample its run's NTC) ----
Step 6  NTC-aware interpretation    call_taxa.py -> calls.tsv + taxa_evidence.tsv
Step 7a Validation summary          compare_expected (validation mode only)
Step 7b Routine summary             summarize_routine (routine mode only)
        ---- batch gather: BuildRunSummary ----
Step 8  Run summary                 run_summary.tsv (batch only)
```

### ASCII data flow diagram

```
 FASTQ R1 + R2
      |
      v
 [Step 1] ncbi_scrub_pe          (if use_human_scrub=true)
      |
      v
 [Step 2] fastp                  clean_r1 / clean_r2
      |
      +------------------------------------+
      |                                    |
      v                                    v
 [Step 3] Centrifuger            [Step 4] minimap2 (16S panel)
      |                                    |
      v                                    v
 parse_centrifuge_kreport.py      [Step 5] extract_rick16s_metrics.py
      |                                    |
  genus_counts.tsv              align_metrics.tsv
      |                                    |
      +----[ BuildNTCBackground ]-----------+
      |      (NTC/NC samples only;
      |       one TSV per run_id)
      v
 ntc_background_<run_id>.tsv (per run)
      |
      v
 [ MatchNTCBackground ]          (maps each sample to its run's NTC file)
      |
      v
 per_sample ntc_background.tsv
      |
      v
 [Step 6] call_taxa.py           calls.tsv
                                  taxa_evidence.tsv  (NEW)
      |
      +---------------------+
      |                     |
      v                     v
 [Step 7a]               [Step 7b]
 validation mode         routine mode
 CompareExpected         SummarizeRoutine
      |                     |
 validation_summary.tsv  routine_summary.tsv
      |                     |
      +-----[ BuildRunSummary (batch only) ]-----+
                                                 |
                                          run_summary.tsv
```

### Batch workflow phases

The batch workflow (`AFI_16S_Batch`) uses a two-scatter, four-gather design:

**Phase 1 scatter** — runs steps 1–5 for all samples in parallel. As samples complete, each NTC/NC sample conditionally exposes its `align_metrics.tsv` and `genus_counts.tsv`.

**BuildNTCBackground gather** — collects all NTC/NC metrics files and computes one `ntc_background_<run_id>.tsv` **per distinct `run_id`**, ensuring samples from different sequencing runs are never cross-contaminated by each other's NTC signal.

**MatchNTCBackground gather** — maps each sample in the batch to the NTC background file for its `run_id`, producing one matched NTC file per sample. This task fails if any `run_id` has no NTC/NC sample — every run submitted in a batch must include at least one negative control.

**Phase 2 scatter** — runs steps 6–7 for all samples in parallel, each using its matched `ntc_background.tsv`.

**BuildRunSummary gather** — merges all per-sample summaries and `calls.tsv` files into a single `run_summary.tsv` annotated with per-run IDs and run-level PC8 validity.

---

## 3. Reference Files and Databases

### 16S Rickettsiales panel

The pipeline includes a curated 16S rRNA reference panel in the repository under `align_rickettsiales_16S/refs/`. Three files are available:

| File | Description |
|---|---|
| `rickettsiales_panel_16S.fa` | 16S reference sequences, original (31 KB) |
| `rickettsiales_panel_16S.clean.fa` | Cleaned version (recommended for use) |
| `rickettsiales_panel_16S.mmi` | Pre-built minimap2 index (226 KB) |

For Terra runs, upload one of these files to your GCS bucket and provide its `gs://` path as `rickettsiales_panel`. You may supply either the FASTA (`.fa`) or the pre-built index (`.mmi`); minimap2 will build the index at runtime from FASTA if needed.

**Recommendation:** Use `rickettsiales_panel_16S.clean.fa` or `rickettsiales_panel_16S.mmi` for all runs. The pre-built `.mmi` index avoids index-build overhead at the cost of a small additional file transfer.

### Centrifuger database

The Centrifuger database is a combined bacteria/archaea + Rickettsiales index. Because the database is large (~67 GiB compressed), it is stored in GCS as a TAR.GZ archive and extracted at runtime by the `RunCentrifuger` task.

You must provide two inputs:

| Input | Example value | Description |
|---|---|---|
| `centrifuger_db` | `centrifuger_bact_arch_plus_rickettsiales` | Index prefix string (no path) |
| `centrifuger_db_archives` | `["gs://bucket/centrifuger_index.tar.gz"]` | Array of GCS paths to TAR.GZ archives |

The task extracts all archives into a local directory, searches for a file matching `<prefix>.1.cfr` or `<prefix>.1.cf`, and derives the full local prefix from that path before invoking `centrifuger -x`.

The database must include both the standard bacteria/archaea content and the Rickettsiales-specific sequences so that a single kreport covers all organisms needed for PC8 validity checking.

**VM sizing:** The Centrifuger task requires substantial compute resources. The recommended defaults are:

```
centrifuger_memory = "96G"
centrifuger_disks  = "local-disk 375 HDD"
classify_threads   = 8
```

The default Terra VM (typically 1 CPU / 2 GB RAM) is far too small and will cause the task to fail or segfault. The 96 GB RAM allocation comfortably covers the ~90–110 GB in-memory footprint of the combined database, and 375 GB disk provides ample space for the ~67 GiB compressed archive plus extracted index. These parameters are exposed at the workflow level so they can be set in the input JSON.

### NTC background file

The NTC background (`ntc_background.tsv`) encodes the maximum read counts observed in no-template control samples for each genus, separately for the alignment and Centrifuger evidence sources. These thresholds prevent false positive calls caused by contamination or bleed-through.

In the **batch workflow**, one NTC background file is computed automatically per distinct `run_id` from all NTC/NC samples in that run. You do not need to supply or manage this file manually — it is an intermediate output of `BuildNTCBackground` and consumed internally by `MatchNTCBackground`.

In the **single-sample workflow**, you must supply an `ntc_background` file manually.

#### Format

**Dual-column format** (produced by `BuildNTCBackground` in the batch workflow):

```
genus	align_ntc_reads	cfr_ntc_reads
Orientia	0	0
Rickettsia	12	45
Leptospira	0	230
```

**Legacy single-column format** (produced by `build_ntc_background_from_metrics.py` — single-sample use):

```
genus	mapped_reads
Orientia	0
Rickettsia	12
```

`call_taxa.py` accepts both formats. When the legacy format is provided, `mapped_reads` is used as `align_ntc_reads` and `cfr_ntc_reads` defaults to 0.

#### Placeholder file

For single-sample runs where no NTC data is available yet, use the placeholder file at `wdl/inputs/ntc_background.placeholder.tsv`. It contains only the header with no data rows, meaning all NTC thresholds default to zero and no NTC correction is applied:

```
genus	mapped_reads
```

Upload this file to GCS and use its path as `ntc_background` in the single-sample workflow input JSON.

---

## 4. Importing Workflows into Terra

### Option 1: Dockstore (recommended)

Importing via Dockstore avoids import resolution errors because Terra pulls the full descriptor set (main WDL + all task imports) automatically.

1. In Dockstore, link the GitHub repository `PHemarajata/afi_terra`. The `.dockstore.yml` at the repository root registers both workflows:
   - `/wdl/AFI_16S_Main.wdl`
   - `/wdl/AFI_16S_Batch.wdl`
2. Create or refresh a workflow version from the desired Git tag or branch.
3. In your Terra workspace, navigate to **Workflows** and click **Find a Workflow**.
4. Select **Dockstore** as the source, search for the workflow name, and import via TRS (Tool Registry Service).
5. Terra resolves all imports (task WDLs, NCBI scrub WDL) automatically from the repository.

### Option 2: Direct ZIP upload (fallback)

If Dockstore is not available, upload a ZIP archive to Terra that preserves relative directory paths. The ZIP must contain at minimum:

```
wdl/AFI_16S_Batch.wdl
wdl/AFI_16S_Main.wdl
wdl/tasks/align.wdl
wdl/tasks/classify.wdl
wdl/tasks/interpret.wdl
wdl/tasks/metrics.wdl
wdl/tasks/preprocess.wdl
wdl/tasks/validate.wdl
NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl
```

> **Important:** Uploading only `AFI_16S_Batch.wdl` alone will fail because it imports `AFI_16S_Main.wdl` and the task WDLs using relative paths.

To create the ZIP from the repository root:

```bash
zip -r afi_terra_wdl.zip \
  wdl/AFI_16S_Batch.wdl \
  wdl/AFI_16S_Main.wdl \
  wdl/tasks/ \
  NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl
```

Then in Terra, navigate to **Workflows**, click **+**, choose **Upload WDL**, and upload the ZIP.

---

## 5. Running the Single-Sample Workflow

The single-sample workflow (`AFI_16S_Main`) processes one pair of FASTQ files through the full pipeline and produces per-sample outputs. Use this workflow for ad-hoc testing or when running individual samples outside of a sequencing run context.

> **Note:** For production clinical runs, use the batch workflow instead. The single-sample workflow requires a pre-computed `ntc_background.tsv` to be supplied manually, whereas the batch workflow computes this automatically.

### All inputs

#### Required inputs

| Input | Type | Description |
|---|---|---|
| `sample_id` | String | Unique sample identifier; used in output filenames |
| `r1_fastq` | File | GCS path to R1 FASTQ (.fastq.gz) |
| `r2_fastq` | File | GCS path to R2 FASTQ (.fastq.gz) |
| `rickettsiales_panel` | File | GCS path to 16S reference FASTA or .mmi index |
| `ntc_background` | File | GCS path to NTC background TSV; use placeholder for first pass |
| `centrifuger_db` | String | Centrifuger index prefix |
| `centrifuger_db_archives` | Array[File] | GCS paths to TAR.GZ archives containing the Centrifuger index |

#### Mode and type inputs

| Input | Type | Default | Allowed values | Description |
|---|---|---|---|---|
| `mode` | String | `"routine"` | `routine`, `validation` | Determines summary output type |
| `sample_type` | String | `"clinical"` | `NTC`, `NC`, `PC_MIX8`, `PC_SINGLE`, `MIXED4`, `clinical`, `PC` | Sample classification |
| `expected_taxon` | String | `""` | Delimited list of genera | Required for validation mode; delimiters: `;`, `,`, or `|` |
| `use_human_scrub` | Boolean | `true` | `true`, `false` | Whether to run NCBI SRA Scrubber for human read removal |

#### Detection threshold inputs

| Input | Type | Default | Description |
|---|---|---|---|
| `align_confirm_reads` | Int | `100` | Minimum mapped reads for Confirmed alignment call |
| `align_confirm_breadth` | Float | `0.25` | Minimum breadth of coverage for Confirmed alignment call |
| `align_fold` | Float | `5.0` | Minimum fold over NTC reads for Confirmed alignment call |
| `cfr_floor` | Int | `500` | Minimum Centrifuger reads for Detected call |
| `cfr_fold` | Float | `5.0` | Minimum fold over NTC reads for Detected call |

#### Compute resource inputs

| Input | Type | Default | Description |
|---|---|---|---|
| `centrifuger_memory` | String | `"96G"` | Memory allocation for Centrifuger task |
| `centrifuger_disks` | String | `"local-disk 375 HDD"` | Disk allocation for Centrifuger task |
| `classify_threads` | Int | `8` | CPU threads for Centrifuger |

#### Docker image inputs

| Input | Default image | Description |
|---|---|---|
| `afi_core_docker` | `phemarajata614/afi-terra:0.4.1` | Core analysis image (Python, samtools, minimap2, scripts) |
| `fastp_docker` | `staphb/fastp:0.23.4` | QC trimming |
| `minimap_docker` | `phemarajata614/afi-terra:0.4.1` | Alignment (minimap2 bundled in afi-terra image) |
| `centrifuger_docker` | `phemarajata614/centrifuger:1.1.0` | Centrifuger classifier |

### Outputs

| Output | Type | Description |
|---|---|---|
| `clean_r1`, `clean_r2` | File | QC-trimmed reads |
| `centrifuger_classification` | File | Raw Centrifuger classification TSV |
| `centrifuger_kreport` | File | Centrifuger kreport |
| `centrifuger_genus_counts` | File | Parsed genus-level counts |
| `minimap_bam`, `minimap_bai` | File | Sorted, indexed alignment BAM |
| `align_metrics` | File | Per-genus alignment metrics |
| `calls` | File | Final taxa calls TSV |
| `taxa_evidence` | File | Per-taxon evidence justification TSV |
| `validation_summary` | File? | Concordance summary (validation mode only) |
| `routine_summary` | File? | Detected taxa list (routine mode only) |
| `scrubbed_r1`, `scrubbed_r2` | File? | Dehosted reads (if `use_human_scrub=true`) |

### Example input JSON (routine)

```json
{
  "AFI_16S_Main.sample_id": "SAMPLE001",
  "AFI_16S_Main.sample_type": "clinical",
  "AFI_16S_Main.mode": "routine",
  "AFI_16S_Main.use_human_scrub": true,
  "AFI_16S_Main.r1_fastq": "gs://YOUR_BUCKET/fastq/SAMPLE001_R1.fastq.gz",
  "AFI_16S_Main.r2_fastq": "gs://YOUR_BUCKET/fastq/SAMPLE001_R2.fastq.gz",
  "AFI_16S_Main.rickettsiales_panel": "gs://YOUR_BUCKET/ref/rickettsiales_panel_16S.clean.fa",
  "AFI_16S_Main.ntc_background": "gs://YOUR_BUCKET/ref/ntc_background.placeholder.tsv",
  "AFI_16S_Main.centrifuger_db": "centrifuger_bact_arch_plus_rickettsiales",
  "AFI_16S_Main.centrifuger_db_archives": [
    "gs://YOUR_BUCKET/db/centrifuger_index.tar.gz"
  ],
  "AFI_16S_Main.centrifuger_memory": "96G",
  "AFI_16S_Main.centrifuger_disks": "local-disk 375 HDD",
  "AFI_16S_Main.classify_threads": 8,
  "AFI_16S_Main.afi_core_docker": "phemarajata614/afi-terra:0.4.1",
  "AFI_16S_Main.centrifuger_docker": "phemarajata614/centrifuger:1.1.0"
}
```

### Example input JSON (validation mode)

```json
{
  "AFI_16S_Main.sample_id": "PC_MIX8_001",
  "AFI_16S_Main.sample_type": "PC_MIX8",
  "AFI_16S_Main.mode": "validation",
  "AFI_16S_Main.expected_taxon": "Orientia;Rickettsia;Leptospira;Burkholderia;Anaplasma;Ehrlichia;Coxiella;Bartonella",
  "AFI_16S_Main.use_human_scrub": true,
  "AFI_16S_Main.r1_fastq": "gs://YOUR_BUCKET/fastq/PC001_R1.fastq.gz",
  "AFI_16S_Main.r2_fastq": "gs://YOUR_BUCKET/fastq/PC001_R2.fastq.gz",
  "AFI_16S_Main.rickettsiales_panel": "gs://YOUR_BUCKET/ref/rickettsiales_panel_16S.clean.fa",
  "AFI_16S_Main.ntc_background": "gs://YOUR_BUCKET/ref/ntc_background.tsv",
  "AFI_16S_Main.centrifuger_db": "centrifuger_bact_arch_plus_rickettsiales",
  "AFI_16S_Main.centrifuger_db_archives": [
    "gs://YOUR_BUCKET/db/centrifuger_index.tar.gz"
  ],
  "AFI_16S_Main.centrifuger_memory": "96G",
  "AFI_16S_Main.centrifuger_disks": "local-disk 375 HDD"
}
```

### Launching in Terra

1. In your Terra workspace, navigate to **Workflows** and select `AFI_16S_Main`.
2. Choose **Run workflow with inputs defined by file paths**.
3. Upload or paste your input JSON.
4. Click **Run Analysis**.

---

## 6. Running the Batch Workflow

The batch workflow (`AFI_16S_Batch`) processes one or more sequencing runs in a single Terra submission. It scatters each sample through the same processing steps as the single-sample workflow, automatically builds a per-run NTC background from the NTC/NC samples in each run, and produces a consolidated run summary.

> **Multi-run support:** A single batch submission can include samples from multiple sequencing runs. Assign each sample a `run_id` matching its run; NTC backgrounds are computed independently per `run_id` so high contamination in one run cannot inflate thresholds for another.

### Inputs

The batch workflow takes parallel arrays of per-sample values (WDL 1.0 compatible). All arrays must be the same length and in the same sample order.

#### Run-level inputs

| Input | Type | Default | Description |
|---|---|---|---|
| `rickettsiales_panel` | File | — | GCS path to 16S reference FASTA or .mmi |
| `centrifuger_db` | String | `""` | Centrifuger index prefix |
| `centrifuger_db_archives` | Array[File] | `[]` | GCS paths to Centrifuger TAR.GZ archives |
| `use_human_scrub` | Boolean | `true` | Run-wide human read removal switch |
| `classify_threads` | Int | `8` | CPU threads for Centrifuger tasks |

#### Per-sample arrays (one value per sample, in matching order)

| Input | Type | Description |
|---|---|---|
| `run_ids` | Array[String] | Sequencing run identifier per sample (used for per-run NTC grouping) |
| `sample_ids` | Array[String] | Sample identifiers |
| `r1_fastqs` | Array[File] | R1 FASTQ GCS paths |
| `r2_fastqs` | Array[File] | R2 FASTQ GCS paths |
| `sample_types` | Array[String] | Sample type for each sample |
| `modes` | Array[String] | `routine` or `validation` for each sample |
| `expected_taxa` | Array[String] | Expected taxa string for each sample; `""` for routine samples |

#### Threshold and resource inputs

These are identical to the single-sample workflow. See Section 5 for the full table. Key inputs:

- `align_confirm_reads` (default: 100)
- `align_confirm_breadth` (default: 0.25)
- `align_fold` (default: 5.0)
- `cfr_floor` (default: 500)
- `cfr_fold` (default: 5.0)
- `centrifuger_memory` (default: `"96G"`)
- `centrifuger_disks` (default: `"local-disk 375 HDD"`)

#### Docker image inputs

Same as single-sample workflow (see Section 5). Override at the workflow level to change the image used by all tasks.

### Outputs

| Output | Type | Description |
|---|---|---|
| `calls` | Array[File] | Final taxa calls per sample |
| `taxa_evidence_files` | Array[File] | Per-taxon evidence justification per sample |
| `align_metrics` | Array[File] | Alignment metrics per sample |
| `centrifuger_genus_counts` | Array[File] | Centrifuger genus counts per sample |
| `ntc_backgrounds_per_run` | Array[File] | Auto-computed NTC background, one per distinct `run_id` |
| `ntc_background_run_ids` | Array[String] | Run IDs corresponding to `ntc_backgrounds_per_run` |
| `validation_summaries` | Array[File?] | Per-sample concordance summaries (validation mode) |
| `routine_summaries` | Array[File?] | Per-sample taxa lists (routine mode) |
| `run_summary` | File | Batch-level summary with per-run PC8 validity |

### Example batch input JSON

```json
{
  "AFI_16S_Batch.run_ids": [
    "run1",
    "run1",
    "run1",
    "run1"
  ],
  "AFI_16S_Batch.sample_ids": [
    "SAMPLE001_S1",
    "SAMPLE002_S2",
    "NTC_S11",
    "PC_S12"
  ],
  "AFI_16S_Batch.r1_fastqs": [
    "gs://YOUR_BUCKET/fastq/SAMPLE001_S1_R1.fastq.gz",
    "gs://YOUR_BUCKET/fastq/SAMPLE002_S2_R1.fastq.gz",
    "gs://YOUR_BUCKET/fastq/NTC_S11_R1.fastq.gz",
    "gs://YOUR_BUCKET/fastq/PC_S12_R1.fastq.gz"
  ],
  "AFI_16S_Batch.r2_fastqs": [
    "gs://YOUR_BUCKET/fastq/SAMPLE001_S1_R2.fastq.gz",
    "gs://YOUR_BUCKET/fastq/SAMPLE002_S2_R2.fastq.gz",
    "gs://YOUR_BUCKET/fastq/NTC_S11_R2.fastq.gz",
    "gs://YOUR_BUCKET/fastq/PC_S12_R2.fastq.gz"
  ],
  "AFI_16S_Batch.sample_types": ["clinical", "clinical", "NTC", "PC_MIX8"],
  "AFI_16S_Batch.modes": ["routine", "routine", "routine", "routine"],
  "AFI_16S_Batch.expected_taxa": ["", "", "", ""],
  "AFI_16S_Batch.rickettsiales_panel": "gs://YOUR_BUCKET/ref/rickettsiales_panel_16S.clean.fa",
  "AFI_16S_Batch.centrifuger_db": "centrifuger_bact_arch_plus_rickettsiales",
  "AFI_16S_Batch.centrifuger_db_archives": [
    "gs://YOUR_BUCKET/db/centrifuger_index.tar.gz"
  ],
  "AFI_16S_Batch.centrifuger_memory": "96G",
  "AFI_16S_Batch.centrifuger_disks": "local-disk 375 HDD",
  "AFI_16S_Batch.use_human_scrub": true,
  "AFI_16S_Batch.classify_threads": 8,
  "AFI_16S_Batch.afi_core_docker": "phemarajata614/afi-terra:0.4.1",
  "AFI_16S_Batch.centrifuger_docker": "phemarajata614/centrifuger:1.1.0"
}
```

> **NTC background:** The batch workflow computes the NTC background automatically. Every distinct `run_id` in the batch **must** have at least one sample with `sample_type` of `NTC` or `NC`; otherwise `MatchNTCBackground` will fail.

### Sample types in a batch run

| sample_type | Description | Contributes to NTC background |
|---|---|---|
| `NTC` | No-template control | Yes |
| `NC` | Negative control (alias for NTC) | Yes |
| `PC_MIX8` | 8-organism positive control | No |
| `PC_SINGLE` | Single-organism positive control | No |
| `MIXED4` | 4-organism mixed control | No |
| `clinical` | Clinical patient sample | No |
| `PC` | Generic positive control | No |

### Launching the batch workflow in Terra

1. In your Terra workspace, navigate to **Workflows** and select `AFI_16S_Batch`.
2. Choose **Run workflow with inputs defined by file paths**.
3. Upload or paste your input JSON.
4. Click **Run Analysis**.

Alternatively, use a Terra data model with a `sample_set` entity and map workflow inputs to `this.samples.*` attributes per the column mapping in the batch WDL header comments.

---

## 7. The Two-Pass Pattern

The NTC background thresholds reflect contamination levels unique to each sequencing run. The **batch workflow** (`AFI_16S_Batch`) handles this entirely automatically — no manual intervention is needed. The two-pass manual procedure below applies only to the **single-sample workflow** or exceptional reprocessing scenarios.

### Batch workflow: fully automatic NTC handling

`AFI_16S_Batch` implements multi-phase processing in a single Terra submission:

1. **Phase 1** processes all samples through alignment and classification.
2. **BuildNTCBackground** computes one `ntc_background_<run_id>.tsv` per distinct `run_id`, using only the NTC/NC samples from that run.
3. **MatchNTCBackground** assigns each sample its run's NTC background file.
4. **Phase 2** applies NTC-aware interpretation to all samples in parallel.

No manual steps are required. Every `run_id` in the submission must have at least one NTC or NC sample.

### Single-sample workflow: manual two-pass

When using `AFI_16S_Main` or reprocessing a sample with a custom NTC background:

**Step 1: First pass — run with placeholder NTC background**

Supply `wdl/inputs/ntc_background.placeholder.tsv` (all-zero thresholds) as `ntc_background`. Include all clinical, control, and NTC samples in the same submission. Collect the `align_metrics` output for each NTC/NC sample from Terra.

**Step 2: Build run-specific NTC background**

```bash
python3 scripts/build_ntc_background_from_metrics.py \
  --ntc-metrics path/to/NTC_S11.align_metrics.tsv \
  --ntc-metrics path/to/NTC_S12.align_metrics.tsv \
  --out wdl/inputs/run3.ntc_background.tsv
```

This takes the per-genus maximum across all NTC metrics files. The output uses the legacy `genus, mapped_reads` format accepted by `call_taxa.py`.

**Step 3: Upload and re-run**

```bash
gsutil cp wdl/inputs/run3.ntc_background.tsv gs://YOUR_BUCKET/ref/run3_ntc_background.tsv
```

Re-submit the same sample with `ntc_background` pointing to the real file. The second pass produces NTC-corrected calls.

### Notes

- Always include NTC/NC and PC_MIX8 in each run so controls are processed within the same run context.
- In `auto` mode-policy (via `build_batch_inputs_json.py`), rows with `expected_results` become `validation` mode; rows without become `routine`.
- MIXED4 and PC_SINGLE are always treated as validation mode in auto mode.
- PC_MIX8 can be used in both validation mode (with `expected_taxon`) and routine mode (without).

---

## 8. Terra Sheet Builder (GUI)

The **Terra Sheet Builder** is a desktop GUI application that generates a Terra-compatible sample-set TSV for `AFI_16S_Batch` without manually editing JSON or TSV files. It is the recommended way to prepare batch inputs for routine clinical use.

### Getting the application

**Pre-built binaries** are built automatically by GitHub Actions on every push to `tools/terra_sheet_builder/`. Download the artifact for your platform from **Actions → Build Terra Sheet Builder** on the repository page:

| Platform | Artifact |
|---|---|
| Linux (amd64) | Single-file ELF executable |
| Windows (x64) | `.exe`, no console window |
| macOS Intel | `.app` bundle, zipped |
| macOS Apple Silicon (arm64) | `.app` bundle, zipped |

**Run from source** (requires Python 3.10+):

```bash
pip install PySide6
python3 tools/terra_sheet_builder/terra_sheet_builder.py
```

### Screen 1 — FASTQ selection

A three-step wizard guides you through entering run information before any sample-level editing.

**Step 1 — Run count.** Enter the number of sequencing runs to include in this batch (spinner, 1–50).

**Step 2 — Run names.** Enter a unique name for each run. Names are used as the `run_id` value for all samples in that run. No spaces are allowed; underscores and hyphens are fine (e.g., `run_2024_11`).

**Step 3 — FASTQ files per run.** For each run, choose one of three methods to add samples:

| Method | How to use |
|---|---|
| **Browse folder** | Select a directory. The app auto-discovers all R1/R2 FASTQ pairs using the `_R1_`/`_R2_` naming convention. |
| **Add files…** | Multi-file dialog. Select any number of FASTQ files. The app auto-pairs R1+R2 by name pattern. If you select exactly two files that do not match the auto-pair pattern, the app falls back to prompting you to assign R1 and R2 manually. |
| **Import from TSV…** | Global import (applies to all runs at once). Upload a mapping TSV with at minimum `run_id` and `sample_id` columns. Optional `r1_fastq`/`r2_fastq` columns can provide file paths directly. If path columns are absent, the app prompts you for a FASTQ folder and performs fuzzy matching to assign files: exact prefix first, then substring matching, then difflib sequence similarity. |

### Screen 2 — Sample metadata

An editable table with columns: `sample_id`, `run_id`, `sample_type`, `mode`, `expected_taxa`.

**Auto-fill:** Changing `sample_type` via the dropdown automatically populates `mode` and `expected_taxa`:

| sample_type | mode | expected_taxa (auto-filled) |
|---|---|---|
| `PC_MIX8` | `validation` | `Bacillus;Listeria;Staphylococcus;Enterococcus;Limosilactobacillus;Salmonella;Escherichia;Pseudomonas` |
| `MIXED4` | `validation` | 4-organism string (configurable) |
| `PC_SINGLE` | `validation` | First organism option |
| `clinical` | `routine` | (empty) |
| `NTC` / `NC` / `PC` | `routine` | (empty) |

**Header fields** (above the table):
- **Analysis Date** — date string added to the output TSV comment field
- **Table Name** — lowercase alphanumeric name (≤32 chars, must start with a letter); used as the Terra entity type name
- **Operator Initials** — recorded in `analysis_comments`

**Validate button.** Checks the following before enabling export:
- No duplicate `sample_id` values across all runs
- Each run contains at least one `NTC` or `NC` sample
- Each run contains at least one positive control (`PC_MIX8`, `PC_SINGLE`, `MIXED4`, or `PC`)
- All validation-mode rows have a non-empty `expected_taxa` field
- Table name matches the required format (`^[a-z][a-z0-9_]{0,31}$`)

Any failed check is highlighted with a descriptive message. Fix the issues and click **Validate** again.

**Export TSV** (enabled only after validation passes). Saves a Terra-compatible TSV to disk.

### Output TSV format

The exported file uses `entity:<table_name>_id` as the first column (required by Terra) followed by all sample fields. An `analysis_comments` column is appended as the last column, containing `<table_name> | <date> | <initials>` for traceability.

**Importing into Terra:**

1. In your Terra workspace, go to **Data** → **Import Data** → **Upload TSV**.
2. Select the exported TSV file.
3. Terra will create or update a sample set entity table with the name you provided.

---

## 9. Helper Scripts

All helper scripts are in the `scripts/` directory. They run locally (not in Terra) and are used to prepare inputs, post-process outputs, and manage Docker images.

### `build_batch_inputs_json.py`

Converts a sample sheet TSV or an AFI mapping TSV into a Terra-ready batch input JSON. This is the primary way to prepare inputs for `AFI_16S_Batch`.

**Usage with a standard sample sheet:**

```bash
python3 scripts/build_batch_inputs_json.py \
  --sample-sheet wdl/inputs/batch_samples.template.tsv \
  --out-json wdl/inputs/batch_run1.json \
  --rickettsiales-panel gs://YOUR_BUCKET/ref/rickettsiales_panel_16S.clean.fa \
  --centrifuger-db centrifuger_bact_arch_plus_rickettsiales
```

**Usage with an AFI mapping TSV (e.g., `AFI_optimizeProtocol.tsv`):**

```bash
python3 scripts/build_batch_inputs_json.py \
  --mapping-tsv AFI_optimizeProtocol.tsv \
  --run-id 3 \
  --out-json wdl/inputs/run3.batch.json \
  --fastq-uri-prefix gs://YOUR_BUCKET/fastq \
  --rickettsiales-panel gs://YOUR_BUCKET/ref/rickettsiales_panel_16S.clean.fa \
  --centrifuger-db centrifuger_bact_arch_plus_rickettsiales \
  --mode-policy auto \
  --per-sample-use-human-scrub true
```

**All options:**

| Option | Required | Description |
|---|---|---|
| `--sample-sheet PATH` | One of these two | Standard TSV sample sheet |
| `--mapping-tsv PATH` | One of these two | AFI_optimizeProtocol.tsv format |
| `--out-json PATH` | Yes | Output JSON path |
| `--rickettsiales-panel URI` | No | GCS URI for 16S panel |
| `--centrifuger-db STRING` | No | Centrifuger index prefix |
| `--mode-policy auto|all_validation|all_routine` | No (default: `auto`) | How to assign mode to samples |
| `--fastq-uri-prefix URI` | Required with `--mapping-tsv` | GCS prefix for FASTQ files |
| `--run-id ID` | No (repeatable) | Filter rows by run_id column |
| `--per-sample-use-human-scrub true|false` | No | Set per-sample use_human_scrub |
| `--no-default-use-human-scrub` | No | Disable run-wide human scrub default |
| `--classify-threads N` | No (default: 8) | Centrifuger thread count |

**Sample sheet format** (`--sample-sheet`):

The template is at `wdl/inputs/batch_samples.template.tsv`. Columns:

```
sample_id    sample_type    mode    classifier_mode    r1_fastq    r2_fastq    expected_taxon    expected_taxa    use_human_scrub
```

- `expected_taxon` and `expected_taxa` are mutually exclusive; provide only one.
- `expected_taxa` takes a semicolon-delimited list (e.g., `Orientia;Rickettsia`).
- `use_human_scrub` accepts `true`, `false`, `1`, `0`, `yes`, `no`, `t`, `f`.
- Omit `use_human_scrub` to inherit the run-wide default.

Example rows:

```
PC001	PC_MIX8	validation	single	gs://bucket/PC001_R1.fastq.gz	gs://bucket/PC001_R2.fastq.gz		Orientia;Rickettsia;Leptospira;Burkholderia	true
SAMPLE001	clinical	routine	single	gs://bucket/SAMPLE001_R1.fastq.gz	gs://bucket/SAMPLE001_R2.fastq.gz			true
NTC001	NTC	routine	single	gs://bucket/NTC001_R1.fastq.gz	gs://bucket/NTC001_R2.fastq.gz			true
```

### `build_ntc_background_from_metrics.py`

Builds a run-specific NTC background TSV from one or more NTC `align_metrics.tsv` output files. Use this for manual (single-sample) two-pass workflows when you collect NTC metrics from Terra and want to create a custom background.

```bash
python3 scripts/build_ntc_background_from_metrics.py \
  --ntc-metrics path/to/NTC1.align_metrics.tsv \
  --ntc-metrics path/to/NTC2.align_metrics.tsv \
  --out wdl/inputs/run_ntc_background.tsv
```

- `--ntc-metrics` is repeatable; pass one flag per NTC metrics file.
- Input columns required: `genus`, `mapped_reads`.
- Output columns: `genus`, `mapped_reads` (legacy format; compatible with `call_taxa.py`).
- The script takes the maximum `mapped_reads` value per genus across all input files.

### `build_push_afi_core_image.sh`

Builds and pushes the AFI core Docker image to Docker Hub.

```bash
bash scripts/build_push_afi_core_image.sh phemarajata614 0.4.1 linux/amd64
```

Arguments:
1. Docker Hub username / organization
2. Image tag
3. Platform (optional; defaults to `linux/amd64`)

Equivalent manual commands:

```bash
docker build --platform linux/amd64 -t phemarajata614/afi-terra:0.4.1 .
docker push phemarajata614/afi-terra:0.4.1
# Optional latest tag:
docker tag phemarajata614/afi-terra:0.4.1 phemarajata614/afi-terra:latest
docker push phemarajata614/afi-terra:latest
```

To verify the image after building:

```bash
docker run --rm phemarajata614/afi-terra:0.4.1 bash -lc \
  "minimap2 --version && samtools --version | head -n 1 && python3 -c 'import pandas; print(pandas.__version__)'"
```

### `make_terra_import_sheet.py`

Generates or validates a Terra-compatible sample import TSV, and optionally writes an annotated Excel workbook. This is the CLI companion to the Terra Sheet Builder GUI (Section 8).

```bash
# Generate a blank template with two example runs:
python3 scripts/make_terra_import_sheet.py --template --output my_run.tsv

# Validate and convert an existing CSV/TSV:
python3 scripts/make_terra_import_sheet.py --input samples.csv --output terra_import.tsv

# Also write an annotated Excel workbook:
python3 scripts/make_terra_import_sheet.py --input samples.csv --output terra_import.tsv --excel
```

**Options:**

| Option | Description |
|---|---|
| `--template` | Write a template file with two pre-filled example runs |
| `--input PATH` | Input CSV or TSV (delimiter auto-detected) |
| `--output PATH` | Output Terra TSV path |
| `--excel` | Also write an annotated Excel workbook (requires openpyxl) |

**Validation rules enforced:**
1. Each `run_id` must contain at least one `NTC` or `NC` sample.
2. Each `run_id` must contain at least one positive control (`PC_MIX8`, `PC_SINGLE`, `MIXED4`, or `PC`).
3. Validation-mode samples must have a non-empty `expected_taxa` field.
4. `sample_type` and `mode` values must be from the allowed sets.

**Required input columns:** `sample_id`, `run_id`, `r1_fastq`, `r2_fastq`, `sample_type`, `mode`, `expected_taxa`.

**Output:** Terra TSV with `entity:sample_id` as the first column. Import via **Data → Import Data → Upload TSV** in your Terra workspace.

### `compare_single_double_outputs.py`

> **Archive only.** This script was used to compare single-classifier (Centrifuger) vs. double-classifier (Kraken2 + Centrifuger) runs during pipeline validation. The double-classifier mode has been removed; this script is retained for reference only.

Accepts either direct paths to summary TSV files or parent directories (searched recursively). Outputs:

| File | Description |
|---|---|
| `<prefix>.summary.tsv` | Per-sample comparison summary |
| `<prefix>.validation_differences.tsv` | Samples where validation results differ |
| `<prefix>.routine_differences.tsv` | Samples where routine taxa lists differ |

---

## 10. Output Files Reference

### Per-sample outputs (both workflows)

These files are produced for every sample regardless of mode.

#### `<sample_id>.centrifuger.classification.tsv`

Raw Centrifuger classification output. One row per read pair. This is the direct output of `centrifuger -x` and uses the standard Centrifuger column format. For most use cases, the parsed genus counts file is more useful.

#### `<sample_id>.centrifuger.kreport.tsv`

Kreport-format output from `centrifuger-kreport`. Follows the Kraken2 kreport column format (percentage, clade reads, taxon reads, rank, taxid, name). Used as input to `parse_centrifuge_kreport.py`.

#### `<sample_id>.genus_counts.tsv`

Genus-level aggregation of Centrifuger kreport reads. Columns:

```
genus    reads
```

Only rank-G (genus) rows are retained. Clade read counts are used so that reads assigned to child species are included under the genus.

#### `<sample_id>.bam` and `<sample_id>.bai`

minimap2 alignment of clean reads against the 16S Rickettsiales panel, sorted and indexed with samtools. Used by `extract_rick16s_metrics.py` to compute per-genus alignment metrics.

#### `<sample_id>.align_metrics.tsv`

Per-genus alignment metrics extracted from the BAM file. Columns:

```
genus    mapped_reads    max_breadth
```

- `mapped_reads`: number of reads mapping to reference sequences for this genus
- `max_breadth`: maximum fraction of any reference sequence covered by aligned reads (0.0–1.0)

This file is the input to `call_taxa.py` for the alignment evidence source. For NTC/NC samples, this file is also collected to build the NTC background.

#### `calls.tsv`

NTC-aware taxa calls combining alignment and Centrifuger evidence. This is the core interpretive output. See Section 11 for column descriptions and call value meanings.

#### `taxa_evidence.tsv`

Detailed per-genus evidence table produced by `InterpretCalls`. Provides the full data behind each call with human-readable justification strings. One row per genus per evidence source.

| Column | Description |
|---|---|
| `sample` | Sample identifier |
| `genus` | Genus name |
| `source` | Evidence source: `alignment` or `centrifuge` |
| `call` | Detection call (same values as `calls.tsv`) |
| `reads` | Reads supporting this call |
| `breadth` | Alignment breadth fraction (alignment rows only; blank for centrifuge) |
| `ntc_reads` | NTC background reads for this genus and source |
| `ntc_fold` | `reads / ntc_reads` (fold above NTC background) |
| `cfr_reads` | Centrifuger reads for this genus (populated on alignment rows for cross-reference) |
| `align_confirmed` | `true` if reads ≥ confirmation threshold AND breadth ≥ confirmation threshold, ignoring NTC fold; used by rescue logic |
| `rescued` | `true` if this alignment call is Confirmed or Probable via rescue when it would otherwise have been Non_Confirmed |
| `evidence_summary` | Human-readable justification string, e.g. `ALIGNMENT CONFIRMED: 150 reads (threshold >=100), breadth 0.31 (threshold >=0.25), 7.5x NTC` |

### Mode-specific per-sample outputs

#### `routine_summary.tsv` (routine mode only)

Produced by `SummarizeRoutineTaxa`. One row per sample. Columns:

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `sample_type` | Sample type string |
| `taxa_present` | Comma-separated list of detected genera (positive calls only) |
| `n_taxa_present` | Count of detected genera |
| `routine_positive_control` | `true` if sample_type is PC_MIX8 or PC; `false` otherwise |

#### `validation_summary.tsv` (validation mode only)

Produced by `CompareExpectedConcordance`. One row per sample. Columns:

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `sample_type` | Sample type string |
| `expected_taxon` | Raw `expected_taxon` input string |
| `expected_taxa` | Comma-separated parsed expected taxa list |
| `expected_taxa_count` | Number of expected taxa |
| `detected_taxa` | Comma-separated list of all detected genera |
| `detected_expected_taxa` | Intersection of detected and expected taxa |
| `detected_expected_count` | Count of detected expected taxa |
| `missing_expected_taxa` | Expected taxa not detected |
| `unexpected_detected_taxa` | Detected taxa not in expected list |
| `validation_result` | `Concordant`, `Discordant`, or `Not_applicable` |
| `pc8_pass` | `true`, `false`, or `not_applicable` |
| `validation_control_class` | Control type classification or `none` |
| `order_rescue_taxa` | Expected taxa that were rescued via order-level alignment evidence (see Section 11) |
| `rescue_mechanisms` | Semicolon-separated list of rescue mechanisms applied, e.g. `order_alignment_rescue;centrifuge_rescue` |

### Batch-only outputs

#### `ntc_background.tsv` (batch workflow only)

Auto-computed NTC background, produced by `BuildNTCBackground`. Columns:

```
genus    align_ntc_reads    cfr_ntc_reads
```

Both Orientia and Rickettsia are always present, even if not seen in NTC samples (with zero values). This file can be downloaded from Terra and reused in future single-sample runs for the same NTC batch context.

#### `ntc_backgrounds_per_run` (batch workflow only)

An array of NTC background TSV files — one file per distinct `run_id` in the batch. Produced by `BuildNTCBackground`. Each file covers only the NTC/NC samples from that run. The companion output `ntc_background_run_ids` (Array[String]) lists the run identifiers in the same order.

These files can be downloaded and reused as the `ntc_background` input for future single-sample runs from the same sequencing run.

#### `taxa_evidence_files` (batch workflow only)

Array of `taxa_evidence.tsv` files, one per sample in the batch. See the per-sample `taxa_evidence.tsv` description above for column details.

#### `run_summary.tsv` (batch workflow only)

Run-level summary produced by `BuildRunSummary`. One row per sample. Columns:

| Column | Description |
|---|---|
| `run_id` | Run identifier for this sample (from `run_ids` array input) |
| `sample_id` | Sample identifier |
| `sample_type` | Sample type |
| `expected_taxa` | Expected taxa string for validation samples; empty for routine samples |
| `detected_taxa` | Comma-separated list of all detected genera (positive calls) |
| `n_detected` | Count of detected genera |
| `validation_result` | From per-sample validation summary, or empty for routine samples |
| `pc8_pass` | From per-sample validation summary, or empty for routine samples |
| `run_pc8_valid` | Run-level PC8 validity: value of `pc8_pass` from the PC_MIX8 sample, or `no_pc8_in_run` if none |
| `rescue_mechanisms` | Rescue mechanisms applied to this sample (see Section 11); empty if none |

### Optional outputs (conditionally present)

| Output | Condition |
|---|---|
| `scrubbed_r1`, `scrubbed_r2` | Only when `use_human_scrub=true` |
| `validation_summary` | Only when `mode=validation` |
| `routine_summary` | Only when `mode=routine` |

---

## 11. Interpreting Calls and Summaries

### The `calls.tsv` file

The `calls.tsv` output is the central result for each sample. It has one row per genus evaluated, combining alignment and Centrifuger evidence with NTC-corrected thresholds.

**Columns:**

| Column | Description |
|---|---|
| `sample` | Sample identifier |
| `genus` | Genus name |
| `source` | Evidence source: `alignment` or `centrifuge` |
| `reads` | Mapped reads (alignment) or classified reads (centrifuge) |
| `breadth` | Fraction of reference covered (alignment only; empty for centrifuge) |
| `ntc_reads` | NTC threshold for this genus and source |
| `cfr_reads` | Centrifuger read count for this genus (on alignment rows; blank on centrifuge rows) |
| `align_confirmed` | `true` if alignment reads and breadth meet confirmation thresholds regardless of NTC fold; blank on centrifuge rows |
| `rescued` | `true` if this call was elevated to Confirmed/Probable via a rescue mechanism; `false` otherwise |
| `call` | Detection call (see below) |

**Call values by source:**

*Alignment source (Orientia and Rickettsia only):*

| Call | Meaning |
|---|---|
| `Confirmed` | Strong positive: reads and breadth meet confirmation thresholds and are well above NTC |
| `Probable` | Moderate positive: reads and breadth meet lower thresholds and reads exceed NTC |
| `Not_Confirmed` | Equivocal: sufficient reads but at or below NTC level; cannot distinguish signal from NTC noise |
| `Negative` | Negative: insufficient reads |

*Centrifuge source (all genera except Orientia and Rickettsia):*

| Call | Meaning |
|---|---|
| `Detected` | Positive: reads exceed floor threshold and are well above NTC |
| `Not_Detected` | Negative or below threshold |

**Note on Orientia/Rickettsia:** Both genera always appear in `calls.tsv` from the alignment source, even if zero reads were observed (they receive a `Negative` call). This ensures consistent output structure across all samples.

**Note on genus prioritization:** If Centrifuger detects Orientia or Rickettsia, those rows are skipped in the centrifuge section because alignment is the authoritative evidence source for these genera. Alignment supersedes classifier calls for these two genera.

### Positive calls summary

"Positive" calls across both sources are `Confirmed`, `Probable`, and `Detected`. These are the calls that appear in `taxa_present` / `detected_taxa` in the summary files.

`Not_Confirmed` is not a positive call. It indicates that reads were observed but they are not distinguishable from NTC noise.

### Validation results

In validation mode, `validation_result` compares the set of detected genera against the expected taxa list:

| Value | Condition |
|---|---|
| `Concordant` | All expected taxa were detected (no missing expected taxa) |
| `Discordant` | At least one expected taxon was not detected |
| `Not_applicable` | Sample type not in validation types, or no expected taxa provided |

Unexpected detections (genera detected but not expected) do not affect the `validation_result`. They are recorded in `unexpected_detected_taxa` for review.

### PC8 validity

For `PC_MIX8` samples, a PC8 pass/fail flag is computed:

| Value | Condition |
|---|---|
| `true` | At least 6 of the expected taxa were detected |
| `false` | Fewer than 6 of the expected taxa were detected |
| `not_applicable` | Sample is not PC_MIX8 or has no expected taxa |

A typical PC_MIX8 control has 8 expected organisms. Detecting at least 6 is considered a passing run.

In the batch run summary, `run_pc8_valid` carries the `pc8_pass` value from the PC_MIX8 sample in the run (or `no_pc8_in_run` if no PC_MIX8 sample was included). This flag annotates every row in `run_summary.tsv` so users can immediately assess run validity regardless of which sample row they are examining.

### Routine positive controls

In routine mode summaries, `routine_positive_control` is `true` for samples with `sample_type` of `PC_MIX8` or `PC`. This allows downstream processing to distinguish positive controls from clinical samples when filtering the routine summary.

### Control type handling

The pipeline distinguishes several control types:

| sample_type | Validation control class | Routine PC flag |
|---|---|---|
| `PC_MIX8` | `PC_MIX8` (8-organism mix) | `true` |
| `MIXED4` | `MIXED4` (4-organism mix) | `false` |
| `PC_SINGLE` | `PC_SINGLE` (single organism) | `false` |
| `PC` | `PC` (generic positive control) | `true` |
| `NTC` | `none` | `false` |
| `NC` | `none` | `false` |
| `clinical` | `none` | `false` |

### Rickettsiales rescue mechanisms

High NTC contamination or low-abundance signals can cause a Concordant result to appear Discordant when the NTC fold gate suppresses an otherwise valid alignment call. Two rescue sub-cases address this for expected *Orientia* or *Rickettsia* taxa.

#### Order-level alignment rescue

If an expected *Orientia* or *Rickettsia* genus is missing from the detected set (would cause a Discordant result), but **any** Rickettsiales alignment row in `calls.tsv` has `align_confirmed=true` (reads ≥ confirmation threshold AND breadth ≥ confirmation threshold, ignoring the NTC fold gate), that expected genus is moved from "missing" to `order_rescue_taxa` in `validation_summary.tsv`. The sample is then called `Concordant` if all other expected taxa are also detected or rescued.

`align_confirmed` represents the raw signal strength independent of NTC noise. A high-NTC run may suppress the final call to `Not_Confirmed`, but if the read and breadth evidence is strong, rescue recovers the concordance.

#### Centrifuge rescue

If Centrifuger independently detects Rickettsiales-adjacent genera (e.g., *Anaplasma*, *Ehrlichia*, *Neorickettsia*, *Wolbachia*) with cfr_reads ≥ 500, or detects *Orientia*/*Rickettsia* reads at that level, this corroborating evidence can also rescue an expected genus from the missing list.

#### Rescue columns in output files

| Column | File | Description |
|---|---|---|
| `align_confirmed` | `calls.tsv`, `taxa_evidence.tsv` | `true` if raw alignment evidence meets confirmation thresholds (ignores NTC fold) |
| `rescued` | `calls.tsv`, `taxa_evidence.tsv` | `true` if call was elevated via rescue |
| `order_rescue_taxa` | `validation_summary.tsv` | Expected taxa rescued by order-level alignment evidence |
| `rescue_mechanisms` | `validation_summary.tsv`, `run_summary.tsv` | Semicolon-separated list of mechanisms applied |

Rescue is only applied in **validation mode**. Routine samples are not subject to rescue logic since there is no expected taxa list to compare against.

---

## 12. Detection Thresholds

### Alignment thresholds (Orientia and Rickettsia)

The alignment-based call for a genus is determined by comparing the mapped reads and breadth of coverage against both absolute thresholds and the NTC baseline. The thresholds below are defaults; all can be overridden via workflow inputs.

**Confirmed** (strong positive):
```
reads >= align_confirm_reads (default: 100)
AND breadth >= align_confirm_breadth (default: 0.25)
AND reads >= align_fold * align_ntc_reads (default: 5.0 x NTC)
```

**Probable** (moderate positive):
```
reads >= 50
AND breadth >= 0.20
AND reads > align_ntc_reads
```

**Not_Confirmed** (equivocal):
```
reads >= 50
AND reads <= align_ntc_reads
```

**Negative**:
```
reads < 50
```

The tier logic is evaluated top-to-bottom; the first matching tier is assigned.

### Centrifuger thresholds (all other genera)

**Detected** (positive):
```
reads >= cfr_floor (default: 500)
AND reads >= cfr_fold * cfr_ntc_reads (default: 5.0 x NTC)
```

**Not_Detected** (negative): all other cases.

### NTC read values

`align_ntc_reads` and `cfr_ntc_reads` come from the `ntc_background.tsv` file. For genera not present in the NTC background file, these values default to 0, which means the NTC fold requirement has no practical effect (any read count satisfies 5.0 x 0 = 0).

In the batch workflow, NTC background values represent the maximum reads observed for each genus across all NTC/NC samples in the run, separately for alignment and Centrifuger evidence. Taking the maximum is conservative: it uses the worst observed NTC level as the threshold.

### Threshold tuning guidance

The defaults are appropriate for most clinical 16S metagenomic applications. Consider adjusting thresholds in the following situations:

| Situation | Suggested adjustment |
|---|---|
| High NTC contamination in runs | The NTC background mechanism handles this automatically; no threshold change needed |
| Very low-abundance clinical samples | Lower `align_confirm_reads` to 50 and `cfr_floor` to 250 |
| High specificity required | Raise `align_confirm_reads` to 200 or `align_confirm_breadth` to 0.30 |
| Evaluation / validation runs | Use default thresholds to compare against expected results |

All threshold parameters are exposed at the workflow level in both `AFI_16S_Main` and `AFI_16S_Batch`.

---

## 13. Docker Images

The pipeline uses four Docker images. All images are pulled from public registries at runtime by Terra.

| Image | Registry | Used by |
|---|---|---|
| `phemarajata614/afi-terra:0.4.1` | Docker Hub | Core analysis: metrics extraction, interpretation, validation/routine summaries, NTC background build, kreport parsing, minimap2 alignment |
| `staphb/fastp:0.23.4` | Docker Hub | QC trimming (Step 2) |
| `phemarajata614/centrifuger:1.1.0` | Docker Hub | Centrifuger classification (Step 3) |
| `us-docker.pkg.dev/general-theiagen/ncbi/sra-human-scrubber:2.2.1` | Google Artifact Registry | Human read removal (Step 1) |

> **Note:** minimap2 is now bundled inside `afi-terra:0.4.1`. The separate `staphb/minimap2:2.28` image is no longer used. Both `afi_core_docker` and `minimap_docker` workflow inputs now default to `phemarajata614/afi-terra:0.4.1`.

### AFI core image (`phemarajata614/afi-terra:0.4.1`)

The core image is built from the `Dockerfile` at the repository root. It contains:

- Python 3.11
- pandas
- samtools
- minimap2
- All Python scripts from `scripts/` installed to `/opt/afi/scripts/`

The following scripts are embedded in the image and invoked directly by WDL task `command` blocks:

- `/opt/afi/scripts/parse_centrifuge_kreport.py`
- `/opt/afi/scripts/extract_rick16s_metrics.py`
- `/opt/afi/scripts/call_taxa.py`
- `/opt/afi/scripts/build_ntc_background.py`

### Overriding image versions

All image inputs are exposed at the workflow level with default values. Override them in your input JSON to use a different version without modifying the WDL:

```json
{
  "AFI_16S_Batch.afi_core_docker": "phemarajata614/afi-terra:0.4.1",
  "AFI_16S_Batch.minimap_docker": "phemarajata614/afi-terra:0.4.1",
  "AFI_16S_Batch.fastp_docker": "staphb/fastp:0.23.4",
  "AFI_16S_Batch.centrifuger_docker": "phemarajata614/centrifuger:1.1.0"
}
```

> **Note:** The NCBI scrubber image is not exposed as a workflow-level input; it is hardcoded in the imported `task_ncbi_scrub.wdl`. To use a different scrubber image, modify that task file before importing.

### Building and pushing a custom core image

If you have made local changes to Python scripts or the Dockerfile:

```bash
# Using the helper script
bash scripts/build_push_afi_core_image.sh YOUR_DOCKERHUB_USER 0.5.0 linux/amd64

# Or manually
docker build --platform linux/amd64 -t YOUR_DOCKERHUB_USER/afi-terra:0.5.0 .
docker push YOUR_DOCKERHUB_USER/afi-terra:0.5.0
```

Then update your input JSON to reference the new image tag.

---

## 14. Troubleshooting

### Centrifuger task fails with "segfault" or exits abnormally

**Cause:** Terra launched the Centrifuger task on a VM that is too small for the database. The Centrifuger database requires approximately 96 GB RAM to load.

**Fix:** Ensure your input JSON sets the Centrifuger resource parameters explicitly:

```json
{
  "AFI_16S_Batch.centrifuger_memory": "96G",
  "AFI_16S_Batch.centrifuger_disks": "local-disk 375 HDD"
}
```

These are workflow-level inputs that override the task defaults. Without them, Terra may select a default VM class (as small as 1 CPU / 2 GB RAM) which is insufficient.

### Centrifuger task fails with "Could not find localized Centrifuger index"

**Cause:** The TAR.GZ archive was extracted but the index files do not match the `centrifuger_db` prefix you provided.

**Fix:** Verify that the extracted archive contains files named `<your_prefix>.1.cfr` or `<your_prefix>.1.cf`. The value of `centrifuger_db` must exactly match the index prefix inside the archive (without path or extension). For example, if the archive extracts to `centrifuger_bact_arch_plus_rickettsiales.1.cfr`, then:

```json
{
  "AFI_16S_Batch.centrifuger_db": "centrifuger_bact_arch_plus_rickettsiales"
}
```

### `fastp: command not found` in FastpClean task

**Cause:** The `fastp_docker` image tag in your JSON does not resolve to a valid image.

**Fix:** Override the fastp docker image explicitly to use the correct staphb image:

```json
{
  "AFI_16S_Batch.fastp_docker": "staphb/fastp:0.23.4"
}
```

Note that `fastp` is run by the `FastpClean` task which uses `fastp_docker`, not `afi_core_docker`. The AFI core image provides Python scripts, samtools, and minimap2; fastp has its own dedicated image.

### Import errors when uploading WDL to Terra

**Cause:** The batch WDL imports `AFI_16S_Main.wdl`, which in turn imports task WDLs and the NCBI scrub WDL using relative paths. Uploading a single WDL file breaks these relative imports.

**Fix:** Use the Dockstore import method (Section 4, Option 1). If uploading a ZIP, ensure all imported files are included with the correct relative directory structure. The ZIP must preserve:

```
wdl/AFI_16S_Batch.wdl
wdl/AFI_16S_Main.wdl
wdl/tasks/*.wdl
NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl
```

### All calls show `Negative` or `Not_Detected` despite expected signal

**Cause 1:** The NTC background placeholder has non-zero values that are suppressing signals. Check the `ntc_background.tsv` used in the run.

**Cause 2:** The FASTQ files are empty or the human scrubber removed too many reads. Check the `clean_r1`/`clean_r2` outputs to see how many reads remain after fastp.

**Cause 3:** The Centrifuger database does not include the target organism. Verify that the database was built with Rickettsiales-specific sequences.

**Fix:** Inspect the `align_metrics.tsv` and `genus_counts.tsv` intermediate outputs to see raw read counts before the calling thresholds are applied. If raw counts are present but calls are negative, adjust thresholds or review the NTC background values.

### Validation result is `Discordant` for expected organisms

**Possible causes:**

1. The expected organism is present in `align_metrics.tsv` or `genus_counts.tsv` but fell below the detection threshold. Check the raw read counts.
2. The `expected_taxon` string uses a genus name that does not exactly match the genus name in the calls output. Genus names are case-sensitive.
3. The NTC background threshold for that genus is too high relative to the signal in this sample.

**Fix:** Cross-reference `calls.tsv` with the expected taxa list. Check that the genus spelling and capitalization match exactly. For NTC threshold issues, review `ntc_background.tsv` and consider whether the NTC sample had unusually high contamination for that genus.

### `Not_Confirmed` calls instead of `Confirmed`

**Cause:** The genus had sufficient reads but the NTC background for that genus was high, causing the reads-to-NTC fold requirement to fail.

**This is intentional behavior.** `Not_Confirmed` means the signal cannot be distinguished from NTC background noise. Investigate whether the NTC sample had genuine contamination or whether the sample was processed in a batch with an unusually contaminated NTC.

If you believe the NTC background is inflated (e.g., from a single outlier NTC run), consider reprocessing with a different NTC background file or manually editing the `ntc_background.tsv` to use more representative values.

### Terra run is very slow or times out

**Cause:** The Centrifuger classification step is the bottleneck. It requires substantial CPU and RAM.

**Fix:** Ensure `classify_threads` is set to at least 8 (default) and that the VM has enough RAM. The default resource allocation in the task is:

```
memory: 96G
disks:  local-disk 375 HDD
cpu:    8
```

Increasing `classify_threads` to 16 can speed up classification if you have a larger VM available.

### `run_pc8_valid` is `no_pc8_in_run`

**Cause:** No sample with `sample_type = PC_MIX8` was included in the batch run.

**Fix:** Always include a PC_MIX8 positive control sample in every batch run submission. The run summary's `run_pc8_valid` field depends on the presence of this control to assess run quality.

### Building NTC background from metrics produces unexpected genera

**Cause:** Centrifuger or alignment detected genera in the NTC sample that are not real organisms. This can happen with low-complexity NTC samples due to index noise.

**Note:** The NTC background is intentionally conservative — it takes the maximum observed reads per genus, so noisy NTC detections translate directly into higher thresholds. This reduces false positives in clinical samples at the cost of potentially masking very low-abundance true positives.

If a specific genus has an unreasonably high NTC background, you can manually edit the downloaded `ntc_background.tsv` file to use more representative values, then supply the edited file as the `ntc_background` input in a single-sample rerun.

### `MatchNTCBackground` fails: "no NTC/NC found for run_id X"

**Cause:** Every `run_id` that appears in the batch's `run_ids` array must have at least one `NTC` or `NC` sample with the same `run_id`. If a run_id appears only in clinical samples, `BuildNTCBackground` cannot compute a background for it and `MatchNTCBackground` will abort.

**Fix:** Check your input sample list and ensure that every distinct `run_id` value has at least one sample with `sample_type = NTC` or `sample_type = NC`. The Terra Sheet Builder (Section 8) and `make_terra_import_sheet.py` (Section 9) both validate this condition before export.

### What does `rescue_mechanisms` mean in `run_summary.tsv`?

The `rescue_mechanisms` column records which (if any) rescue logic was applied during concordance evaluation for validation-mode samples. Values are semicolon-separated:

| Value | Meaning |
|---|---|
| (empty) | No rescue was needed or applied |
| `order_alignment_rescue` | An expected *Orientia* or *Rickettsia* genus was rescued because `align_confirmed=true` on a Rickettsiales alignment row, even though the NTC-gated call was `Not_Confirmed` |
| `centrifuge_rescue` | Centrifuge read evidence (cfr_reads ≥ 500 for Rickettsiales-adjacent genera) corroborated the expected detection |

Rescue mechanisms only apply in validation mode. Routine samples always have an empty `rescue_mechanisms` field. For more details on rescue logic, see Section 11.

