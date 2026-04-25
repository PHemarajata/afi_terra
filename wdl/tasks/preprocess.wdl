version 1.0

task FastpClean {

  input {
    File r1
    File r2
    # staphb/fastp is a lean public container maintained by StaPH-B.
    # It is separate from afi-terra so the core image stays python/samtools only.
    String docker_image = "staphb/fastp:0.23.4"
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
    docker: docker_image
    memory: "8G"
    disks:  "local-disk 100 HDD"
  }
}
