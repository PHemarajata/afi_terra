version 1.0

task FastpClean {

  input {
    File r1
    File r2
  }

  command <<<
  fastp \
    -i ~{r1} \
    -I ~{r2} \
    -o clean_R1.fastq.gz \
    -O clean_R2.fastq.gz
  >>>

  output {
    File clean_r1 = "clean_R1.fastq.gz"
    File clean_r2 = "clean_R2.fastq.gz"
  }

  runtime {
    docker: "phemarajata614/afi-terra:0.1.0"
  }
}
