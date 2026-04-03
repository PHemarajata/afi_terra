version 1.0

# Batch workflow for the AFI Rickettsiales pipeline.
#
# Design: two-scatter with automatic NTC background computation.
#
#   Phase 1 scatter  — dehosting, QC, Centrifuge, 16S alignment, metrics extraction
#   BuildNTCBackground — single gather task; filters NTC/NC samples; builds ntc_background.tsv
#   Phase 2 scatter  — interpretation + per-sample validation/routine summary
#   BuildRunSummary  — final gather; annotates all samples with run_pc8_valid flag
#
# The ntc_background.tsv is fully automatic — no manual 2-pass Terra submit needed.
# All NTC/NC samples in the batch are used to derive run-specific thresholds.

import "tasks/preprocess.wdl" as prep
import "tasks/classify.wdl"   as cls
import "tasks/align.wdl"      as aln
import "tasks/metrics.wdl"    as met
import "tasks/interpret.wdl"  as ipt
import "tasks/validate.wdl"   as vld
import "../NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as scrub

struct SampleSpec {
  String  sample_id
  String  sample_type       # NTC, NC, PC_MIX8, PC_SINGLE, MIXED4, clinical, PC
  String  mode              # validation | routine
  File    r1_fastq
  File    r2_fastq
  String? expected_taxon    # single genus or delimited list for validation mode
  String? expected_taxa     # alias; takes precedence over expected_taxon if both set
  Boolean? use_human_scrub  # per-sample override; defaults to default_use_human_scrub
}

workflow AFI_Rickettsiales_Batch {

  input {
    Array[SampleSpec] samples

    File   rickettsiales_panel      # 16S reference FASTA / pre-built .mmi

    String centrifuger_db = ""
    Array[File] centrifuger_db_archives = []

    # Alignment (Module 3) thresholds — Orientia / Rickettsia
    Int   align_confirm_reads   = 100
    Float align_confirm_breadth = 0.25
    Float align_fold            = 5.0

    # Centrifuge (Module 1) thresholds — all other genera
    Int   cfr_floor = 500
    Float cfr_fold  = 5.0

    Boolean default_use_human_scrub = true
    Int     classify_threads        = 16

    String afi_core_docker    = "phemarajata614/afi-terra:0.4.0"
    String centrifuger_docker  = "phemarajata614/centrifuger:1.1.0"
    String centrifuger_memory  = "128G"
    String centrifuger_disks   = "local-disk 500 HDD"
  }

  # ===========================================================================
  # Phase 1 scatter: preprocessing + classification + alignment + metrics
  # ===========================================================================
  scatter (sample in samples) {

    Boolean do_scrub = select_first([sample.use_human_scrub, default_use_human_scrub])

    if (do_scrub) {
      call scrub.ncbi_scrub_pe as HumanScrub {
        input:
          read1      = sample.r1_fastq,
          read2      = sample.r2_fastq,
          samplename = sample.sample_id
      }
    }

    File p1_r1 = select_first([HumanScrub.read1_dehosted, sample.r1_fastq])
    File p1_r2 = select_first([HumanScrub.read2_dehosted, sample.r2_fastq])

    call prep.FastpClean as P1_Fastp {
      input:
        r1           = p1_r1,
        r2           = p1_r2,
        docker_image = afi_core_docker
    }

    call cls.RunCentrifuger as P1_Centrifuger {
      input:
        sample_id               = sample.sample_id,
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
        sample_id    = sample.sample_id,
        kreport      = P1_Centrifuger.classifier_report_tsv,
        docker_image = afi_core_docker
    }

    call aln.MinimapRick16S as P1_Minimap {
      input:
        r1           = P1_Fastp.clean_r1,
        r2           = P1_Fastp.clean_r2,
        panel        = rickettsiales_panel,
        docker_image = afi_core_docker
    }

    call met.ExtractMetrics as P1_Metrics {
      input:
        bam          = P1_Minimap.bam,
        panel        = rickettsiales_panel,
        docker_image = afi_core_docker
    }

    # Expose metrics files only for NTC/NC samples so we can select_all them
    # outside the scatter without needing the scatter variable (WDL 1.0 safe).
    Boolean p1_is_ntc = (sample.sample_type == "NTC") || (sample.sample_type == "NC")
    if (p1_is_ntc) {
      File ntc_align_conditional = P1_Metrics.metrics
      File ntc_cfr_conditional   = P1_ParseKreport.genus_counts
    }

  } # end Phase 1 scatter

  # ===========================================================================
  # BuildNTCBackground: gather NTC outputs → ntc_background.tsv
  # ===========================================================================
  # ntc_align_conditional and ntc_cfr_conditional are Array[File?] after the
  # scatter; select_all filters to only the NTC/NC samples' files.
  call vld.BuildNTCBackground {
    input:
      ntc_align_metrics    = select_all(ntc_align_conditional),
      ntc_cfr_genus_counts = select_all(ntc_cfr_conditional),
      docker_image         = afi_core_docker
  }

  # ===========================================================================
  # Phase 2 scatter: interpretation + per-sample summaries
  # ===========================================================================
  scatter (i in range(length(samples))) {

    String p2_expected = select_first([
      samples[i].expected_taxa,
      samples[i].expected_taxon,
      ""
    ])

    call ipt.InterpretCalls as P2_Interpret {
      input:
        sample_id             = samples[i].sample_id,
        align_metrics         = P1_Metrics.metrics[i],
        cfr_genus_counts      = P1_ParseKreport.genus_counts[i],
        ntc_background        = BuildNTCBackground.ntc_background,
        align_confirm_reads   = align_confirm_reads,
        align_confirm_breadth = align_confirm_breadth,
        align_fold            = align_fold,
        cfr_floor             = cfr_floor,
        cfr_fold              = cfr_fold,
        docker_image          = afi_core_docker
    }

    if (samples[i].mode == "validation") {
      call vld.CompareExpectedConcordance as P2_Validate {
        input:
          sample_id      = samples[i].sample_id,
          sample_type    = samples[i].sample_type,
          expected_taxon = p2_expected,
          final_calls    = P2_Interpret.calls,
          docker_image   = afi_core_docker
      }
    }

    if (samples[i].mode == "routine") {
      call vld.SummarizeRoutineTaxa as P2_Routine {
        input:
          sample_id    = samples[i].sample_id,
          sample_type  = samples[i].sample_type,
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
    Array[File] centrifuger_kreports      = P1_Centrifuger.classifier_report_tsv
    Array[File] centrifuger_genus_counts  = P1_ParseKreport.genus_counts

    # Phase 1 — alignment
    Array[File] minimap_bam   = P1_Minimap.bam
    Array[File] minimap_bai   = P1_Minimap.bai
    Array[File] align_metrics = P1_Metrics.metrics

    # NTC background (auto-computed)
    File ntc_background = BuildNTCBackground.ntc_background

    # Phase 2 — interpretation
    Array[File] calls = P2_Interpret.calls

    # Phase 2 — per-sample summaries
    Array[File?] validation_summaries = P2_Validate.validation_summary
    Array[File?] routine_summaries    = P2_Routine.routine_summary

    # Run-level summary with PC8 validity
    File run_summary = BuildRunSummary.run_summary
  }
}
