version 1.0

# Centrifuge classification (single classifier — replaces Kraken2 dual-DB mode).
# The database must include standard bacteria/archaea AND Rickettsiales so that
# a single kreport covers all organisms needed for PC8 validity checking and
# general taxonomic screening.

task RunCentrifuger {
  input {
    String sample_id
    File r1_fastq
    File r2_fastq
    String centrifuger_db = ""          # path prefix when archives are NOT used
    Array[File] centrifuger_db_archives = []  # preferred: tar.gz archive(s); Terra localizes these
    Int threads = 8
    String docker_image = "phemarajata614/centrifuger:1.1.0"
    String memory = "96G"
    String disks = "local-disk 375 HDD"
  }

  command <<<
  set -euo pipefail

  db_prefix="~{centrifuger_db}"
  if [[ -n "~{sep=' ' centrifuger_db_archives}" ]]; then
    mkdir -p centrifuger_db
    for archive in ~{sep=' ' centrifuger_db_archives}; do
      tar -xzf "$archive" -C centrifuger_db
    done

    # Auto-detect the index prefix from any .1.cfr / .1.cf file in the
    # extracted directory — no need to know the internal naming convention.
    prefix_file="$(find centrifuger_db -type f \( -name "*.1.cfr" -o -name "*.1.cf" \) | sort | head -n 1 || true)"
    if [[ -z "$prefix_file" ]]; then
      echo "Could not find a Centrifuger index (.1.cfr or .1.cf) after extracting archive(s)." >&2
      exit 1
    fi

    db_prefix="${prefix_file%.1.cfr}"
    db_prefix="${db_prefix%.1.cf}"
  fi

  centrifuger \
    -x "$db_prefix" \
    -1 ~{r1_fastq} \
    -2 ~{r2_fastq} \
    -t ~{threads} \
    > ~{sample_id}.centrifuger.classification.tsv

  centrifuger-kreport \
    -x "$db_prefix" \
    ~{sample_id}.centrifuger.classification.tsv \
    > ~{sample_id}.centrifuger.kreport.tsv
  >>>

  output {
    File classification_tsv    = "~{sample_id}.centrifuger.classification.tsv"
    File classifier_report_tsv = "~{sample_id}.centrifuger.kreport.tsv"
  }

  runtime {
    docker: docker_image
    cpu:    threads
    memory: memory
    disks:  disks
  }
}

# Parse a Centrifuge/Kraken2-style kreport to a simple genus-level TSV.
# Outputs: genus_counts.tsv  (columns: genus, reads)
# Only rank-G rows are kept; clade read count (column 1) is used so that
# reads assigned to child species are included under the genus.
task ParseCentrifugerKreport {
  input {
    String sample_id
    File kreport
    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 /opt/afi/scripts/parse_centrifuge_kreport.py \
    --kreport ~{kreport} \
    --out ~{sample_id}.genus_counts.tsv
  >>>

  output {
    File genus_counts = "~{sample_id}.genus_counts.tsv"
  }

  runtime {
    docker: docker_image
    memory: "4G"
    disks:  "local-disk 20 HDD"
  }
}
