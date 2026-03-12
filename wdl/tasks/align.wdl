version 1.0

task MinimapRick16S {

  input {
    File r1
    File r2
    File panel
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
    docker: "phemarajata614/afi-terra:0.1.0"
  }
}
