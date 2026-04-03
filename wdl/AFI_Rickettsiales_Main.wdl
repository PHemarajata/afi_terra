version 1.0

# Single-sample workflow.
# For full batch runs with automatic NTC background computation,
# use AFI_Rickettsiales_Batch.wdl instead.
#
# When running standalone, supply a pre-computed ntc_background.tsv.
# A placeholder file with all-zero rows is acceptable for first-pass runs
# where you just want to see raw calls before NTC thresholds are applied.

import "tasks/preprocess.wdl"  as prep
import "tasks/classify.wdl"    as cls
import "tasks/align.wdl"       as aln
import "tasks/metrics.wdl"     as met
import "tasks/interpret.wdl"   as ipt
import "tasks/validate.wdl"    as vld
import "../NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as scrub

workflow AFI_Rickettsiales_Main {

  input {
    String  sample_id
    String  sample_type   = "clinical"  # NTC, NC, PC_MIX8, PC_SINGLE, MIXED4, clinical, PC
    String  mode          = "routine"   # routine | validation
    Boolean use_human_scrub = true

    File r1_fastq
    File r2_fastq

    String expected_taxon = ""   # single genus or delimited list (;,|) for validation mode

    File rickettsiales_panel   # 16S reference FASTA (or pre-built .mmi index)
    File ntc_background        # pre-computed ntc_background.tsv; use zeros for first-pass

    String centrifuger_db = ""
    Array[File] centrifuger_db_archives = []

    # Alignment (Module 3) thresholds — Orientia / Rickettsia
    Int   align_confirm_reads   = 100
    Float align_confirm_breadth = 0.25
    Float align_fold            = 5.0

    # Centrifuge (Module 1) thresholds — all other genera
    Int   cfr_floor = 500
    Float cfr_fold  = 5.0

    String afi_core_docker    = "phemarajata614/afi-terra:0.4.1"  # python + samtools + scripts
    String fastp_docker       = "staphb/fastp:0.23.4"             # QC trimming
    String minimap_docker     = "phemarajata614/afi-terra:0.4.1"  # alignment + samtools sort/index
    String centrifuger_docker = "phemarajata614/centrifuger:1.1.0"
    String centrifuger_memory = "128G"
    String centrifuger_disks  = "local-disk 500 HDD"
    Int    classify_threads   = 16
  }

  # -------------------------------------------------------------------------
  # Step 1: Optional human read dehosting
  # -------------------------------------------------------------------------
  if (use_human_scrub) {
    call scrub.ncbi_scrub_pe as HumanScrub {
      input:
        read1      = r1_fastq,
        read2      = r2_fastq,
        samplename = sample_id
    }
  }

  File effective_r1 = select_first([HumanScrub.read1_dehosted, r1_fastq])
  File effective_r2 = select_first([HumanScrub.read2_dehosted, r2_fastq])

  # -------------------------------------------------------------------------
  # Step 2: Quality control / adapter trimming
  # -------------------------------------------------------------------------
  call prep.FastpClean {
    input:
      r1           = effective_r1,
      r2           = effective_r2,
      docker_image = fastp_docker
  }

  # -------------------------------------------------------------------------
  # Step 3: Taxonomic classification (Centrifuge — single classifier)
  # -------------------------------------------------------------------------
  call cls.RunCentrifuger {
    input:
      sample_id             = sample_id,
      r1_fastq              = FastpClean.clean_r1,
      r2_fastq              = FastpClean.clean_r2,
      centrifuger_db        = centrifuger_db,
      centrifuger_db_archives = centrifuger_db_archives,
      threads               = classify_threads,
      docker_image          = centrifuger_docker,
      memory                = centrifuger_memory,
      disks                 = centrifuger_disks
  }

  call cls.ParseCentrifugerKreport {
    input:
      sample_id    = sample_id,
      kreport      = RunCentrifuger.classifier_report_tsv,
      docker_image = afi_core_docker
  }

  # -------------------------------------------------------------------------
  # Step 4: 16S Rickettsiales alignment (confirmatory for Orientia/Rickettsia)
  # -------------------------------------------------------------------------
  call aln.MinimapRick16S {
    input:
      r1           = FastpClean.clean_r1,
      r2           = FastpClean.clean_r2,
      panel        = rickettsiales_panel,
      docker_image = minimap_docker
  }

  call met.ExtractMetrics {
    input:
      bam          = MinimapRick16S.bam,
      panel        = rickettsiales_panel,
      docker_image = afi_core_docker
  }

  # -------------------------------------------------------------------------
  # Step 5: NTC-aware taxa interpretation
  # -------------------------------------------------------------------------
  call ipt.InterpretCalls {
    input:
      sample_id           = sample_id,
      align_metrics       = ExtractMetrics.metrics,
      cfr_genus_counts    = ParseCentrifugerKreport.genus_counts,
      ntc_background      = ntc_background,
      align_confirm_reads = align_confirm_reads,
      align_confirm_breadth = align_confirm_breadth,
      align_fold          = align_fold,
      cfr_floor           = cfr_floor,
      cfr_fold            = cfr_fold,
      docker_image        = afi_core_docker
  }

  # -------------------------------------------------------------------------
  # Step 6: Mode-specific output
  # -------------------------------------------------------------------------
  if (mode == "validation") {
    call vld.CompareExpectedConcordance {
      input:
        sample_id     = sample_id,
        sample_type   = sample_type,
        expected_taxon = expected_taxon,
        final_calls   = InterpretCalls.calls,
        docker_image  = afi_core_docker
    }
  }

  if (mode == "routine") {
    call vld.SummarizeRoutineTaxa {
      input:
        sample_id    = sample_id,
        sample_type  = sample_type,
        final_calls  = InterpretCalls.calls,
        docker_image = afi_core_docker
    }
  }

  output {
    # Dehosting
    File? scrubbed_r1 = HumanScrub.read1_dehosted
    File? scrubbed_r2 = HumanScrub.read2_dehosted

    # QC
    File clean_r1 = FastpClean.clean_r1
    File clean_r2 = FastpClean.clean_r2

    # Classification
    File centrifuger_classification = RunCentrifuger.classification_tsv
    File centrifuger_kreport        = RunCentrifuger.classifier_report_tsv
    File centrifuger_genus_counts   = ParseCentrifugerKreport.genus_counts

    # Alignment
    File minimap_bam = MinimapRick16S.bam
    File minimap_bai = MinimapRick16S.bai
    File align_metrics = ExtractMetrics.metrics

    # Interpretation
    File calls = InterpretCalls.calls

    # Summaries
    File? validation_summary = CompareExpectedConcordance.validation_summary
    File? routine_summary    = SummarizeRoutineTaxa.routine_summary
  }
}
