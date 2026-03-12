version 1.0

import "AFI_Rickettsiales_Main.wdl" as single

struct SampleSpec {
  String sample_id
  String sample_type
  String mode
  String classifier_mode
  File r1_fastq
  File r2_fastq
  String? expected_taxon
  Boolean? use_human_scrub
}

workflow AFI_Rickettsiales_Batch {

  input {
    Array[SampleSpec] samples

    File rickettsiales_panel
    File ntc_background

    String? kraken_db_16g
    String? kraken_db_rick
    String? centrifuger_db

    Boolean default_use_human_scrub = true
    Int classify_threads = 16
  }

  scatter (sample in samples) {
    call single.AFI_Rickettsiales_Main as RunSample {
      input:
        sample_id = sample.sample_id,
        sample_type = sample.sample_type,
        mode = sample.mode,
        classifier_mode = sample.classifier_mode,
        use_human_scrub = select_first([sample.use_human_scrub, default_use_human_scrub]),
        r1_fastq = sample.r1_fastq,
        r2_fastq = sample.r2_fastq,
        expected_taxon = select_first([sample.expected_taxon, ""]),
        rickettsiales_panel = rickettsiales_panel,
        ntc_background = ntc_background,
        kraken_db_16g = kraken_db_16g,
        kraken_db_rick = kraken_db_rick,
        centrifuger_db = centrifuger_db,
        classify_threads = classify_threads
    }
  }

  output {
    Array[File?] scrubbed_r1 = RunSample.scrubbed_r1
    Array[File?] scrubbed_r2 = RunSample.scrubbed_r2
    Array[File] clean_r1 = RunSample.clean_r1
    Array[File] clean_r2 = RunSample.clean_r2
    Array[File] minimap_bam = RunSample.minimap_bam
    Array[File] minimap_bai = RunSample.minimap_bai
    Array[File] metrics = RunSample.metrics
    Array[File] calls = RunSample.calls
    Array[File?] classifier_16g_report = RunSample.classifier_16g_report
    Array[File?] classifier_rick_report = RunSample.classifier_rick_report
    Array[File?] merged_classifier_report = RunSample.merged_classifier_report
    Array[File?] centrifuger_report = RunSample.centrifuger_report
    Array[File?] validation_summary = RunSample.validation_summary
    Array[File?] routine_summary = RunSample.routine_summary
  }
}
