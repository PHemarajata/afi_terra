version 1.0

task MinimapRick16S {

  input {
    File r1
    File r2
    File panel
    # Use the staphb/minimap2 image with the -samtools tag, which bundles
    # both minimap2 and samtools (required for sort/index).
    String docker_image = "staphb/minimap2:2.28-samtools"
  }

  command <<<
  minimap2 -ax sr ~{panel} ~{r1} ~{r2} \
  | samtools sort -o align.bam
  samtools index align.bam
  >>>

  output {
    File bam = "align.bam"
    File bai = "align.bam.bai"
  }

  runtime {
    docker: docker_image
    memory: "16G"
    disks:  "local-disk 100 HDD"
  }
}
