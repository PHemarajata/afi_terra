version 1.0

task MinimapRick16S {

  input {
    File r1
    File r2
    File panel
    # staphb/minimap2 bundles minimap2 + samtools — no need to carry those
    # tools in the afi-terra image.  Maintained by StaPH-B.
    String docker_image = "staphb/minimap2:2.28"
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
