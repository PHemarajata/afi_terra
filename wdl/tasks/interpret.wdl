version 1.0

task InterpretCalls {

  input {
    String sample_id
    File   align_metrics    # from ExtractMetrics  — genus, mapped_reads, max_breadth
    File   cfr_genus_counts # from ParseCentrifugerKreport — genus, reads
    File   ntc_background   # from BuildNTCBackground — genus, align_ntc_reads, cfr_ntc_reads

    # Alignment (Module 3) thresholds — Orientia / Rickettsia
    Int   align_confirm_reads   = 100
    Float align_confirm_breadth = 0.25
    Float align_fold            = 5.0

    # Centrifuge (Module 1) thresholds — all other genera
    Int   cfr_floor = 500
    Float cfr_fold  = 5.0

    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 /opt/afi/scripts/call_taxa.py \
    --sample ~{sample_id} \
    --align-metrics ~{align_metrics} \
    --cfr-genus ~{cfr_genus_counts} \
    --ntc ~{ntc_background} \
    --align-confirm-reads ~{align_confirm_reads} \
    --align-confirm-breadth ~{align_confirm_breadth} \
    --align-fold ~{align_fold} \
    --cfr-floor ~{cfr_floor} \
    --cfr-fold ~{cfr_fold} \
    --out ~{sample_id}.calls.tsv

  python3 - <<'PY'
import csv

CONFIRM_READS   = ~{align_confirm_reads}
CONFIRM_BREADTH = ~{align_confirm_breadth}
ALIGN_FOLD      = ~{align_fold}
CFR_FLOOR       = ~{cfr_floor}
CFR_FOLD        = ~{cfr_fold}
PRESENT_CALLS   = {"Confirmed", "Probable", "Detected"}


def ntc_fold_str(reads, ntc):
    return f"{reads / ntc:.1f}x" if ntc > 0 else "no NTC background"


def build_summary(row):
    source = row["source"]
    call   = row["call"]
    reads  = int(row["reads"]) if row["reads"] else 0
    ntc    = int(row["ntc_reads"]) if row["ntc_reads"] else 0
    fold   = ntc_fold_str(reads, ntc)

    if source == "alignment":
        brd   = float(row["breadth"]) if row["breadth"] else 0.0
        cfr_r = int(row["cfr_reads"]) if row["cfr_reads"] else 0
        a_cfm = str(row.get("align_confirmed", "")).lower()

        if call == "Confirmed":
            s = (
                f"ALIGNMENT CONFIRMED: {reads} reads (threshold >={CONFIRM_READS}), "
                f"breadth {brd:.4f} (threshold >={CONFIRM_BREADTH}). "
                f"NTC background: {ntc} reads; fold above NTC: {fold} "
                f"(threshold >={ALIGN_FOLD}x). "
            )
        else:  # Probable
            s = (
                f"ALIGNMENT PROBABLE: {reads} reads (threshold >=50; "
                f"below confirmed threshold of {CONFIRM_READS}), "
                f"breadth {brd:.4f} (threshold >=0.20). "
                f"NTC background: {ntc} reads; fold above NTC: {fold}. "
            )

        if cfr_r >= CFR_FLOOR:
            s += (f"Centrifuge: {cfr_r} reads (>={CFR_FLOOR} threshold) -- "
                  "centrifuge independently confirms.")
        elif cfr_r > 0:
            s += (f"Centrifuge: {cfr_r} reads (below {CFR_FLOOR} threshold) -- "
                  "below centrifuge threshold; call relies on alignment.")
        else:
            s += ("Centrifuge: 0 reads -- "
                  "not detected by centrifuge; call relies on alignment.")

        if a_cfm == "true":
            s += (" Rickettsiales order rescue eligible "
                  "(align_confirmed=true: breadth/reads thresholds met without NTC gate).")

    else:  # centrifuge
        s = (
            f"CENTRIFUGE DETECTED: {reads} reads (threshold >={CFR_FLOOR}). "
            f"NTC background: {ntc} reads; fold above NTC: {fold} "
            f"(threshold >={CFR_FOLD}x)."
        )

    return s


rows = []
with open("~{sample_id}.calls.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["call"] not in PRESENT_CALLS:
            continue
        reads = int(row["reads"]) if row["reads"] else 0
        ntc   = int(row["ntc_reads"]) if row["ntc_reads"] else 0
        rows.append({
            "sample":           row["sample"],
            "genus":            row["genus"],
            "source":           row["source"],
            "call":             row["call"],
            "reads":            reads,
            "breadth":          row["breadth"],
            "ntc_reads":        ntc,
            "ntc_fold":         ntc_fold_str(reads, ntc),
            "cfr_reads":        row.get("cfr_reads", ""),
            "align_confirmed":  row.get("align_confirmed", ""),
            "rescued":          row.get("rescued", ""),
            "evidence_summary": build_summary(row),
        })

fieldnames = [
    "sample", "genus", "source", "call",
    "reads", "breadth", "ntc_reads", "ntc_fold",
    "cfr_reads", "align_confirmed", "rescued",
    "evidence_summary",
]
with open("~{sample_id}.taxa_evidence.tsv", "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
PY
  >>>

  output {
    File calls        = "~{sample_id}.calls.tsv"
    File taxa_evidence = "~{sample_id}.taxa_evidence.tsv"
  }

  runtime {
    docker: docker_image
    memory: "4G"
    disks:  "local-disk 20 HDD"
  }
}
