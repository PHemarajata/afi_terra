version 1.0

task InterpretCalls {

  input {
    String sample_id
    File metrics
    File ntc_table
    String docker_image = "phemarajata614/afi-terra:0.1.0"
  }

  command <<<
  python3 /opt/afi/scripts/call_taxa.py \
    --sample ~{sample_id} \
    --metrics ~{metrics} \
    --ntc ~{ntc_table} \
    --out calls.tsv
  >>>

  output {
    File calls = "calls.tsv"
  }

  runtime {
    docker: docker_image
  }
}
