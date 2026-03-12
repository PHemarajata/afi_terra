version 1.0

task run_centrifuger_pe {
  input {
    File read1
    File read2
    File index_tar_gz
    String index_prefix
    String sample_id

    String docker_image = "phemarajata614/centrifuger:1.1.0"

    Int threads = 4
    String memory = "16G"
    Int disk_gb = 200
    Int boot_disk_gb = 30
  }

  command <<<
    set -euo pipefail

    echo "===== ENVIRONMENT ====="
    echo "Hostname: $(hostname)"
    echo "Working directory: $(pwd)"
    echo "PATH=$PATH"
    echo "Threads requested: ~{threads}"
    echo

    echo "===== LOCALIZED INPUTS ====="
    ls -lah
    echo
    df -h
    echo

    echo "===== TOOL CHECK ====="
    which centrifuger
    which centrifuger-quant
    which centrifuger-kreport
    echo

    echo "===== UNPACK INDEX ====="
    mkdir -p db
    tar -xzf "~{index_tar_gz}" -C db
    ls -lah db
    echo

    echo "===== RUN CENTRIFUGER ====="
    centrifuger \
      -x "db/~{index_prefix}" \
      -1 "~{read1}" \
      -2 "~{read2}" \
      -t ~{threads} \
      > "~{sample_id}.centrifuger.classification.tsv"

    echo "===== RUN CENTRIFUGER-QUANT ====="
    centrifuger-quant \
      -x "db/~{index_prefix}" \
      -c "~{sample_id}.centrifuger.classification.tsv" \
      > "~{sample_id}.centrifuger.report.tsv"

    echo "===== RUN CENTRIFUGER-KREPORT ====="
    centrifuger-kreport \
      -x "db/~{index_prefix}" \
      "~{sample_id}.centrifuger.classification.tsv" \
      > "~{sample_id}.kreport.tsv"

    echo "===== FINAL OUTPUTS ====="
    ls -lah "~{sample_id}".*
  >>>

  output {
    File classification_tsv = "~{sample_id}.centrifuger.classification.tsv"
    File report_tsv = "~{sample_id}.centrifuger.report.tsv"
    File kreport_tsv = "~{sample_id}.kreport.tsv"
  }

  runtime {
    docker: docker_image
    cpu: threads
    memory: memory
    disks: "local-disk " + disk_gb + " SSD"
    bootDiskSizeGb: boot_disk_gb
  }
}

workflow centrifuger_pe_wf {
  input {
    File read1
    File read2
    File index_tar_gz
    String index_prefix
    String sample_id

    String docker_image = "phemarajata614/centrifuger:1.1.0"

    Int threads = 4
    String memory = "16G"
    Int disk_gb = 200
    Int boot_disk_gb = 30
  }

  call run_centrifuger_pe {
    input:
      read1 = read1,
      read2 = read2,
      index_tar_gz = index_tar_gz,
      index_prefix = index_prefix,
      sample_id = sample_id,
      docker_image = docker_image,
      threads = threads,
      memory = memory,
      disk_gb = disk_gb,
      boot_disk_gb = boot_disk_gb
  }

  output {
    File classification_tsv = run_centrifuger_pe.classification_tsv
    File report_tsv = run_centrifuger_pe.report_tsv
    File kreport_tsv = run_centrifuger_pe.kreport_tsv
  }
}