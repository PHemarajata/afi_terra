version 1.0

# Batch workflow for the AFI Rickettsiales pipeline.
#
# Terra data-model integration
# ────────────────────────────
# Launch this workflow from a Terra **sample_set** entity.  Each member sample
# in the set must have the following columns in the sample table:
#
#   sample_id    (String)  — entity name column
#   run_id       (String)  — sequencing run identifier (can differ per sample
#                            to process multiple runs in one batch submission)
#   r1_fastq     (File)
#   r2_fastq     (File)
#   sample_type  (String)  — NTC | NC | PC_MIX8 | PC_SINGLE | MIXED4 | clinical | PC
#   mode         (String)  — routine | validation
#   expected_taxa (String) — semicolon-delimited genera for validation samples; "" otherwise
#
# Typical Terra input mapping
#   AFI_Rickettsiales_Batch.run_ids       → this.samples.run_id
#   AFI_Rickettsiales_Batch.sample_ids    → this.samples.sample_id
#   AFI_Rickettsiales_Batch.r1_fastqs     → this.samples.r1_fastq
#   AFI_Rickettsiales_Batch.r2_fastqs     → this.samples.r2_fastq
#   AFI_Rickettsiales_Batch.sample_types  → this.samples.sample_type
#   AFI_Rickettsiales_Batch.modes         → this.samples.mode
#   AFI_Rickettsiales_Batch.expected_taxa → this.samples.expected_taxa
#   AFI_Rickettsiales_Batch.use_human_scrub → workspace.use_human_scrub  (or hardcode)
#
# Multi-run support: include samples from multiple run_ids in a single set.
# run_summary.tsv will be grouped by run_id; pc8_valid is computed per run.
#
# Design: two-scatter with automatic NTC background computation.
#
#   Phase 1 scatter  — dehosting, QC, Centrifuge, 16S alignment, metrics extraction
#   BuildNTCBackground — single gather task; NTC files pre-filtered via conditional
#                        declarations inside Phase 1 scatter + select_all (WDL 1.0 safe)
#   Phase 2 scatter  — interpretation + per-sample validation/routine summary
#   BuildRunSummary  — final gather; annotates all samples with run_pc8_valid flag

import "tasks/preprocess.wdl" as prep
import "tasks/classify.wdl"   as cls
import "tasks/align.wdl"      as aln
import "tasks/metrics.wdl"    as met
import "tasks/interpret.wdl"  as ipt
import "tasks/validate.wdl"   as vld
import "../NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as scrub

workflow AFI_Rickettsiales_Batch {

  input {
    # ── Per-sample run identifier (one per sample, same length as sample_ids) ──
    # Map from the sample table run_id column.  Samples may belong to different
    # runs; run_summary.tsv will group rows and compute pc8_valid per run_id.
    Array[String] run_ids

    # ── Per-sample inputs (parallel arrays, same length) ──────────────────────
    # Map from Terra sample table columns.  All arrays must have equal length.
    Array[String] sample_ids    # entity key column
    Array[File]   r1_fastqs
    Array[File]   r2_fastqs
    Array[String] sample_types  # NTC, NC, PC_MIX8, PC_SINGLE, MIXED4, clinical, PC
    Array[String] modes         # routine | validation
    Array[String] expected_taxa # "" for non-validation samples; "Genus1;Genus2" for validation

    # ── Reference files ────────────────────────────────────────────────────────
    File   rickettsiales_panel      # 16S reference FASTA / pre-built .mmi

    String       centrifuger_db          = ""
    Array[File]  centrifuger_db_archives = []

    # ── Alignment (Module 3) thresholds — Orientia / Rickettsia ───────────────
    Int   align_confirm_reads   = 100
    Float align_confirm_breadth = 0.25
    Float align_fold            = 5.0

    # ── Centrifuge (Module 1) thresholds — all other genera ───────────────────
    Int   cfr_floor = 500
    Float cfr_fold  = 5.0

    # ── Run-wide options ───────────────────────────────────────────────────────
    # Single switch applies to every sample; set via workspace attribute or JSON.
    Boolean use_human_scrub   = true
    Int     classify_threads  = 16

    # ── Docker images ──────────────────────────────────────────────────────────
    String afi_core_docker    = "phemarajata614/afi-terra:0.4.1"  # python + samtools + scripts
    String fastp_docker       = "staphb/fastp:0.23.4"             # QC trimming
    String minimap_docker     = "phemarajata614/afi-terra:0.4.1"  # alignment + samtools sort/index
    String centrifuger_docker = "phemarajata614/centrifuger:1.1.0"
    String centrifuger_memory = "128G"
    String centrifuger_disks  = "local-disk 500 HDD"
  }

  # ===========================================================================
  # Phase 1 scatter: preprocessing + classification + alignment + metrics
  # ===========================================================================
  scatter (i in range(length(sample_ids))) {

    String p1_id   = sample_ids[i]
    File   p1_in_r1 = r1_fastqs[i]
    File   p1_in_r2 = r2_fastqs[i]

    if (use_human_scrub) {
      call scrub.ncbi_scrub_pe as HumanScrub {
        input:
          read1      = p1_in_r1,
          read2      = p1_in_r2,
          samplename = p1_id
      }
    }

    File p1_r1 = select_first([HumanScrub.read1_dehosted, p1_in_r1])
    File p1_r2 = select_first([HumanScrub.read2_dehosted, p1_in_r2])

    call prep.FastpClean as P1_Fastp {
      input:
        r1           = p1_r1,
        r2           = p1_r2,
        docker_image = fastp_docker
    }

    call cls.RunCentrifuger as P1_Centrifuger {
      input:
        sample_id               = p1_id,
        r1_fastq                = P1_Fastp.clean_r1,
        r2_fastq                = P1_Fastp.clean_r2,
        centrifuger_db          = centrifuger_db,
        centrifuger_db_archives = centrifuger_db_archives,
        threads                 = classify_threads,
        docker_image            = centrifuger_docker,
        memory                  = centrifuger_memory,
        disks                   = centrifuger_disks
    }

    call cls.ParseCentrifugerKreport as P1_ParseKreport {
      input:
        sample_id    = p1_id,
        kreport      = P1_Centrifuger.classifier_report_tsv,
        docker_image = afi_core_docker
    }

    call aln.MinimapRick16S as P1_Minimap {
      input:
        r1           = P1_Fastp.clean_r1,
        r2           = P1_Fastp.clean_r2,
        panel        = rickettsiales_panel,
        docker_image = minimap_docker
    }

    call met.ExtractMetrics as P1_Metrics {
      input:
        bam          = P1_Minimap.bam,
        panel        = rickettsiales_panel,
        docker_image = afi_core_docker
    }

    # Expose metrics only for NTC/NC samples so select_all() can filter them
    # outside the scatter without referencing the scatter variable (WDL 1.0 safe).
    # Also expose the run_id so BuildNTCBackground can compute per-run backgrounds.
    Boolean p1_is_ntc = (sample_types[i] == "NTC") || (sample_types[i] == "NC")
    if (p1_is_ntc) {
      String ntc_run_id_conditional  = run_ids[i]
      File   ntc_align_conditional   = P1_Metrics.metrics
      File   ntc_cfr_conditional     = P1_ParseKreport.genus_counts
    }

  } # end Phase 1 scatter

  # ===========================================================================
  # BuildNTCBackground: compute one background TSV per distinct run_id.
  # NTC samples from different runs are never pooled together, so a high-
  # signal NTC in run2 cannot raise the detection threshold for run1 samples.
  # ===========================================================================
  call vld.BuildNTCBackground {
    input:
      ntc_run_ids          = select_all(ntc_run_id_conditional),
      ntc_align_metrics    = select_all(ntc_align_conditional),
      ntc_cfr_genus_counts = select_all(ntc_cfr_conditional),
      docker_image         = afi_core_docker
  }

  # Map each sample to the NTC background for its run_id.
  call vld.MatchNTCBackground {
    input:
      all_sample_run_ids  = run_ids,
      per_run_ids         = BuildNTCBackground.per_run_ids,
      per_run_backgrounds = BuildNTCBackground.per_run_backgrounds,
      docker_image        = afi_core_docker
  }

  # ===========================================================================
  # Phase 2 scatter: interpretation + per-sample summaries
  # ===========================================================================
  scatter (i in range(length(sample_ids))) {

    call ipt.InterpretCalls as P2_Interpret {
      input:
        sample_id             = sample_ids[i],
        align_metrics         = P1_Metrics.metrics[i],
        cfr_genus_counts      = P1_ParseKreport.genus_counts[i],
        ntc_background        = MatchNTCBackground.per_sample_backgrounds[i],
        align_confirm_reads   = align_confirm_reads,
        align_confirm_breadth = align_confirm_breadth,
        align_fold            = align_fold,
        cfr_floor             = cfr_floor,
        cfr_fold              = cfr_fold,
        docker_image          = afi_core_docker
    }

    if (modes[i] == "validation") {
      call vld.CompareExpectedConcordance as P2_Validate {
        input:
          sample_id      = sample_ids[i],
          sample_type    = sample_types[i],
          expected_taxon = expected_taxa[i],
          final_calls    = P2_Interpret.calls,
          docker_image   = afi_core_docker
      }
    }

    if (modes[i] == "routine") {
      call vld.SummarizeRoutineTaxa as P2_Routine {
        input:
          sample_id    = sample_ids[i],
          sample_type  = sample_types[i],
          final_calls  = P2_Interpret.calls,
          docker_image = afi_core_docker
      }
    }

  } # end Phase 2 scatter

  # ===========================================================================
  # BuildRunSummary: run-level summary with PC8 validity flag
  # ===========================================================================
  call vld.BuildRunSummary {
    input:
      run_ids              = run_ids,
      calls_files          = P2_Interpret.calls,
      validation_summaries = P2_Validate.validation_summary,
      routine_summaries    = P2_Routine.routine_summary,
      docker_image         = afi_core_docker
  }

  # ===========================================================================
  # Outputs
  # ===========================================================================
  output {
    # Phase 1 — preprocessing
    Array[File?] scrubbed_r1 = HumanScrub.read1_dehosted
    Array[File?] scrubbed_r2 = HumanScrub.read2_dehosted
    Array[File]  clean_r1    = P1_Fastp.clean_r1
    Array[File]  clean_r2    = P1_Fastp.clean_r2

    # Phase 1 — classification
    Array[File] centrifuger_kreports     = P1_Centrifuger.classifier_report_tsv
    Array[File] centrifuger_genus_counts = P1_ParseKreport.genus_counts

    # Phase 1 — alignment
    Array[File] minimap_bam   = P1_Minimap.bam
    Array[File] minimap_bai   = P1_Minimap.bai
    Array[File] align_metrics = P1_Metrics.metrics

    # NTC backgrounds — one per run_id (computed separately to prevent cross-run pooling)
    Array[File]   ntc_backgrounds_per_run = BuildNTCBackground.per_run_backgrounds
    Array[String] ntc_background_run_ids  = BuildNTCBackground.per_run_ids

    # Phase 2 — interpretation
    Array[File] calls = P2_Interpret.calls

    # Phase 2 — per-sample summaries
    Array[File?] validation_summaries = P2_Validate.validation_summary
    Array[File?] routine_summaries    = P2_Routine.routine_summary

    # Run-level summary with run_id + PC8 validity
    File run_summary = BuildRunSummary.run_summary
  }
}
