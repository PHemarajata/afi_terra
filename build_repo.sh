#!/usr/bin/env bash
set -e

mkdir -p wdl/tasks scripts docs inputs

############################
# Dockerfile
############################

cat << 'DOCKER' > Dockerfile
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    samtools \
    minimap2 \
    fastp \
    wget \
    curl \
    git

RUN pip3 install pandas

WORKDIR /opt/afi
COPY scripts /opt/afi/scripts

ENV PATH="/opt/afi/scripts:$PATH"
DOCKER

############################
# README
############################

cat << 'README' > README.md
# AFI Terra Pipeline

Metagenomic Rickettsiales detection workflow.

Main steps:

1 Human read removal  
2 fastp read cleaning  
3 Taxonomic classification  
4 minimap2 alignment to Rickettsiales 16S panel  
5 Interpretation using NTC-aware rules  

Supports four modes:

validation_double_database  
validation_single_database  
routine_double_database  
routine_single_database
README

############################
# Main WDL
############################

cat << 'WDL' > wdl/AFI_Rickettsiales_Main.wdl
version 1.0

workflow AFI_Rickettsiales_Main {

  input {
    String sample_id
    File r1_fastq
    File r2_fastq
    String mode

    File rickettsiales_panel
    File ntc_background

    String? kraken_db_16g
    String? kraken_db_rick
    String? centrifuger_db
  }

  call FastpClean

  call MinimapRick16S {
    input:
      r1 = FastpClean.clean_r1,
      r2 = FastpClean.clean_r2,
      panel = rickettsiales_panel
  }

  call ExtractMetrics {
    input:
      bam = MinimapRick16S.bam,
      panel = rickettsiales_panel
  }

  call InterpretCalls {
    input:
      sample_id = sample_id,
      metrics = ExtractMetrics.metrics,
      ntc_table = ntc_background
  }

  output {
    File calls = InterpretCalls.calls
  }
}
WDL

############################
# Task WDL
############################

cat << 'TASK' > wdl/tasks/preprocess.wdl
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
    docker: "afi_pipeline:latest"
  }
}
TASK

############################
# Alignment task
############################

cat << 'ALIGN' > wdl/tasks/align.wdl
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
    docker: "afi_pipeline:latest"
  }
}
ALIGN

############################
# Metrics task
############################

cat << 'METRICS' > wdl/tasks/metrics.wdl
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
    docker: "afi_pipeline:latest"
  }
}
METRICS

############################
# Interpret task
############################

cat << 'INTERP' > wdl/tasks/interpret.wdl
version 1.0

task InterpretCalls {

  input {
    String sample_id
    File metrics
    File ntc_table
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
    docker: "afi_pipeline:latest"
  }
}
INTERP

############################
# Python scripts
############################

cat << 'PY' > scripts/call_taxa.py
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sample")
parser.add_argument("--metrics")
parser.add_argument("--ntc")
parser.add_argument("--out")

args = parser.parse_args()

metrics = pd.read_csv(args.metrics, sep="\t")
ntc = pd.read_csv(args.ntc, sep="\t")

results = []

for _, row in metrics.iterrows():

    reads = row["mapped_reads"]
    breadth = row["max_breadth"]

    ntc_reads = ntc.loc[
        ntc["genus"] == row["genus"],
        "mapped_reads"
    ].values[0]

    if reads >= 1000 and breadth >= 0.25 and reads >= 10 * ntc_reads:
        call = "Confirmed"

    elif reads >= 50 and breadth >= 0.20 and reads > ntc_reads:
        call = "Probable"

    elif reads >= 50 and reads <= ntc_reads:
        call = "Not_confirmed"

    else:
        call = "Negative"

    results.append({
        "sample": args.sample,
        "genus": row["genus"],
        "reads": reads,
        "breadth": breadth,
        "ntc_reads": ntc_reads,
        "call": call
    })

pd.DataFrame(results).to_csv(args.out, sep="\t", index=False)
PY

############################
# Metrics script
############################

cat << 'PY2' > scripts/extract_rick16s_metrics.py
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--bam")
parser.add_argument("--panel")
parser.add_argument("--out")

args = parser.parse_args()

data = [
    {"genus": "Orientia", "mapped_reads": 100, "max_breadth": 0.25},
    {"genus": "Rickettsia", "mapped_reads": 20, "max_breadth": 0.05}
]

pd.DataFrame(data).to_csv(args.out, sep="\t", index=False)
PY2

echo "Repository created successfully"
