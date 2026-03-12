version 1.0

task RunKraken2_16G {
  input {
    String sample_id
    File r1_fastq
    File r2_fastq
    String kraken_db
    Int threads = 16
    String docker_image = "afi_pipeline:latest"
  }

  command <<<
  kraken2 \
    --db ~{kraken_db} \
    --paired ~{r1_fastq} ~{r2_fastq} \
    --report ~{sample_id}.kk2_16g.report.tsv \
    --output ~{sample_id}.kk2_16g.output.tsv \
    --threads ~{threads}
  >>>

  output {
    File kraken_report = "~{sample_id}.kk2_16g.report.tsv"
    File kraken_output = "~{sample_id}.kk2_16g.output.tsv"
  }

  runtime {
    docker: docker_image
  }
}

task RunKraken2_Rick {
  input {
    String sample_id
    File r1_fastq
    File r2_fastq
    String kraken_db
    Int threads = 16
    String docker_image = "afi_pipeline:latest"
  }

  command <<<
  kraken2 \
    --db ~{kraken_db} \
    --paired ~{r1_fastq} ~{r2_fastq} \
    --report ~{sample_id}.kk2_rick.report.tsv \
    --output ~{sample_id}.kk2_rick.output.tsv \
    --threads ~{threads}
  >>>

  output {
    File kraken_report = "~{sample_id}.kk2_rick.report.tsv"
    File kraken_output = "~{sample_id}.kk2_rick.output.tsv"
  }

  runtime {
    docker: docker_image
  }
}

task MergeDoubleDatabaseReports {
  input {
    String sample_id
    File report_16g
    File report_rick
    String docker_image = "afi_pipeline:latest"
  }

  command <<<
  {
    echo -e "classifier\treport_file"
    echo -e "kraken2_16g\t~{report_16g}"
    echo -e "kraken2_rick\t~{report_rick}"
  } > ~{sample_id}.double_classifier_reports.tsv
  >>>

  output {
    File merged_report = "~{sample_id}.double_classifier_reports.tsv"
  }

  runtime {
    docker: docker_image
  }
}

task RunCentrifuger {
  input {
    String sample_id
    File r1_fastq
    File r2_fastq
    String centrifuger_db
    Int threads = 16
    String docker_image = "phemarajata614/centrifuger:1.1.0"
  }

  command <<<
  centrifuger \
    -x ~{centrifuger_db} \
    -1 ~{r1_fastq} \
    -2 ~{r2_fastq} \
    -t ~{threads} \
    > ~{sample_id}.centrifuger.classification.tsv

  centrifuger-kreport \
    -x ~{centrifuger_db} \
    ~{sample_id}.centrifuger.classification.tsv \
    > ~{sample_id}.centrifuger.kreport.tsv
  >>>

  output {
    File classification_tsv = "~{sample_id}.centrifuger.classification.tsv"
    File classifier_report_tsv = "~{sample_id}.centrifuger.kreport.tsv"
  }

  runtime {
    docker: docker_image
  }
}
