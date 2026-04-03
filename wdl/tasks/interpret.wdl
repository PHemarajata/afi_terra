version 1.0

task InterpretCalls {

  input {
    String sample_id
    File   align_metrics    # from ExtractMetrics  — genus, mapped_reads, max_breadth
    File   cfr_genus_counts # from ParseCentrifugerKreport — genus, reads
    File   ntc_background   # from BuildNTCBackground — genus, align_ntc_reads, cfr_ntc_reads

    # Alignment (Module 3) thresholds — Orientia / Rickettsia
    Int   align_confirm_reads   = 100
    Float align_confirm_breadth = 0.25
    Float align_fold            = 5.0

    # Centrifuge (Module 1) thresholds — all other genera
    Int   cfr_floor = 500
    Float cfr_fold  = 5.0

    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 /opt/afi/scripts/call_taxa.py \
    --sample ~{sample_id} \
    --align-metrics ~{align_metrics} \
    --cfr-genus ~{cfr_genus_counts} \
    --ntc ~{ntc_background} \
    --align-confirm-reads ~{align_confirm_reads} \
    --align-confirm-breadth ~{align_confirm_breadth} \
    --align-fold ~{align_fold} \
    --cfr-floor ~{cfr_floor} \
    --cfr-fold ~{cfr_fold} \
    --out calls.tsv
  >>>

  output {
    File calls = "calls.tsv"
  }

  runtime {
    docker: docker_image
    memory: "4G"
    disks:  "local-disk 20 HDD"
  }
}
