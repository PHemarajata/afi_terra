version 1.0

import "tasks/preprocess.wdl" as prep
import "tasks/classify.wdl" as cls
import "tasks/align.wdl" as aln
import "tasks/metrics.wdl" as met
import "tasks/interpret.wdl" as ipt
import "tasks/validate.wdl" as vld
import "../NCBI_scrub_PE/tasks/quality_control/read_filtering/task_ncbi_scrub.wdl" as scrub

workflow AFI_Rickettsiales_Main {

  input {
    String sample_id
    String sample_type = "clinical"
    String mode = "routine"
    String classifier_mode = "single"
    Boolean use_human_scrub = true
    File r1_fastq
    File r2_fastq

    String expected_taxon = ""

    File rickettsiales_panel
    File ntc_background

    String? kraken_db_16g
    String? kraken_db_rick
    String? centrifuger_db

    Int classify_threads = 16
  }

  if (use_human_scrub) {
    call scrub.ncbi_scrub_pe as HumanScrub {
      input:
        read1 = r1_fastq,
        read2 = r2_fastq,
        samplename = sample_id
    }
  }

  File effective_r1 = select_first([HumanScrub.read1_dehosted, r1_fastq])
  File effective_r2 = select_first([HumanScrub.read2_dehosted, r2_fastq])

  call prep.FastpClean as FastpClean {
    input:
      r1 = effective_r1,
      r2 = effective_r2
  }

  if (classifier_mode == "double") {
    call cls.RunKraken2_16G {
      input:
        sample_id = sample_id,
        r1_fastq = FastpClean.clean_r1,
        r2_fastq = FastpClean.clean_r2,
        kraken_db = select_first([kraken_db_16g]),
        threads = classify_threads
    }

    call cls.RunKraken2_Rick {
      input:
        sample_id = sample_id,
        r1_fastq = FastpClean.clean_r1,
        r2_fastq = FastpClean.clean_r2,
        kraken_db = select_first([kraken_db_rick]),
        threads = classify_threads
    }

    call cls.MergeDoubleDatabaseReports {
      input:
        sample_id = sample_id,
        report_16g = RunKraken2_16G.kraken_report,
        report_rick = RunKraken2_Rick.kraken_report
    }
  }

  if (classifier_mode == "single") {
    call cls.RunCentrifuger {
      input:
        sample_id = sample_id,
        r1_fastq = FastpClean.clean_r1,
        r2_fastq = FastpClean.clean_r2,
        centrifuger_db = select_first([centrifuger_db]),
        threads = classify_threads
    }
  }

  call aln.MinimapRick16S {
    input:
      r1 = FastpClean.clean_r1,
      r2 = FastpClean.clean_r2,
      panel = rickettsiales_panel
  }

  call met.ExtractMetrics {
    input:
      bam = MinimapRick16S.bam,
      panel = rickettsiales_panel
  }

  call ipt.InterpretCalls {
    input:
      sample_id = sample_id,
      metrics = ExtractMetrics.metrics,
      ntc_table = ntc_background
  }

  if (mode == "validation") {
    call vld.CompareExpectedConcordance {
      input:
        sample_id = sample_id,
        sample_type = sample_type,
        expected_taxon = expected_taxon,
        final_calls = InterpretCalls.calls
    }
  }

  if (mode == "routine") {
    call vld.SummarizeRoutineTaxa {
      input:
        sample_id = sample_id,
        sample_type = sample_type,
        final_calls = InterpretCalls.calls
    }
  }

  output {
    File? scrubbed_r1 = HumanScrub.read1_dehosted
    File? scrubbed_r2 = HumanScrub.read2_dehosted
    File clean_r1 = FastpClean.clean_r1
    File clean_r2 = FastpClean.clean_r2
    File minimap_bam = MinimapRick16S.bam
    File minimap_bai = MinimapRick16S.bai
    File metrics = ExtractMetrics.metrics
    File calls = InterpretCalls.calls
    File? classifier_16g_report = RunKraken2_16G.kraken_report
    File? classifier_rick_report = RunKraken2_Rick.kraken_report
    File? merged_classifier_report = MergeDoubleDatabaseReports.merged_report
    File? centrifuger_report = RunCentrifuger.classifier_report_tsv
    File? validation_summary = CompareExpectedConcordance.validation_summary
    File? routine_summary = SummarizeRoutineTaxa.routine_summary
  }
}
