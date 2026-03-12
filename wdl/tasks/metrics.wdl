version 1.0

task ExtractMetrics {

  input {
    File bam
    File panel
  }

  command <<<
  python3 /opt/afi/scripts/extract_rick16s_metrics.py \
    --bam ~{bam} \
    --panel ~{panel} \
    --out metrics.tsv
  >>>

  output {
    File metrics = "metrics.tsv"
  }

  runtime {
    docker: "phemarajata614/afi-terra:0.1.0"
  }
}
