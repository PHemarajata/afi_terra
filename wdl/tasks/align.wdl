version 1.0

task MinimapRick16S {

  input {
    File r1
    File r2
    File panel
    Int threads = 8
    String sort_memory_per_thread = "1G"
    # afi-terra image bundles minimap2 + samtools (required for sort/index).
    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  set -euo pipefail
  minimap2 -t ~{threads} -ax sr ~{panel} ~{r1} ~{r2} \
  | samtools sort -@ ~{threads} -m ~{sort_memory_per_thread} -o align.bam
  samtools index align.bam
  >>>

  output {
    File bam = "align.bam"
    File bai = "align.bam.bai"
  }

  runtime {
    docker: docker_image
    cpu:    threads
    memory: "16G"
    disks:  "local-disk 100 HDD"
  }
}
