version 1.0

# ---------------------------------------------------------------------------
# CompareExpectedConcordance
#   Validation mode: compare detected taxa against expected_taxon list.
#   Also computes pc8_pass for PC_MIX8 samples (TP_rows >= 6).
# ---------------------------------------------------------------------------
task CompareExpectedConcordance {
  input {
    String sample_id
    String sample_type
    String expected_taxon = ""
    File   final_calls
    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 - <<'PY'
import re
import pandas as pd

calls       = pd.read_csv("~{final_calls}", sep="\t")
sample_id   = "~{sample_id}"
sample_type = "~{sample_type}"
expected    = "~{expected_taxon}".strip()

VALIDATION_TYPES = {"PC_MIX8", "MIXED4", "PC_SINGLE", "CLINICAL", "PC"}
POSITIVE_CALLS   = {"Confirmed", "Probable", "Detected"}

# Detected genera (any positive call from either source)
detected = sorted(set(
    calls.loc[calls["call"].isin(POSITIVE_CALLS), "genus"].astype(str)
))

# Parse expected taxa (delimited by ; , |)
expected_list = [x.strip() for x in re.split(r"[;,|]", expected) if x.strip()]

detected_set   = set(detected)
expected_set   = set(expected_list)
detected_exp   = sorted(expected_set & detected_set)
missing_exp    = sorted(expected_set - detected_set)
unexpected_det = sorted(detected_set - expected_set)
tp_rows        = len(detected_exp)

if sample_type.upper() in VALIDATION_TYPES and expected_list:
    validation_result = "Concordant" if len(missing_exp) == 0 else "Discordant"
else:
    validation_result = "Not_applicable"

# PC8 validity: PC_MIX8 must detect >= 6 of its expected organisms
if sample_type.upper() == "PC_MIX8" and expected_list:
    pc8_pass = "true" if tp_rows >= 6 else "false"
else:
    pc8_pass = "not_applicable"

out = pd.DataFrame([{
    "sample_id":               sample_id,
    "sample_type":             sample_type,
    "expected_taxon":          expected,
    "expected_taxa":           ",".join(expected_list),
    "expected_taxa_count":     len(expected_list),
    "detected_taxa":           ",".join(detected),
    "detected_expected_taxa":  ",".join(detected_exp),
    "detected_expected_count": tp_rows,
    "missing_expected_taxa":   ",".join(missing_exp),
    "unexpected_detected_taxa":",".join(unexpected_det),
    "validation_result":       validation_result,
    "pc8_pass":                pc8_pass,
    "validation_control_class": (
        sample_type if sample_type.upper() in {"PC_MIX8", "MIXED4", "PC_SINGLE", "PC"}
        else "none"
    ),
}])

out.to_csv("~{sample_id}.validation_summary.tsv", sep="\t", index=False)
PY
  >>>

  output {
    File   validation_summary = "~{sample_id}.validation_summary.tsv"
  }

  runtime {
    docker: docker_image
    memory: "4G"
    disks:  "local-disk 20 HDD"
  }
}

# ---------------------------------------------------------------------------
# SummarizeRoutineTaxa
#   Routine mode: list detected taxa; flag known positive control types.
# ---------------------------------------------------------------------------
task SummarizeRoutineTaxa {
  input {
    String sample_id
    String sample_type
    File   final_calls
    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 - <<'PY'
import pandas as pd

calls       = pd.read_csv("~{final_calls}", sep="\t")
sample_id   = "~{sample_id}"
sample_type = "~{sample_type}"

POSITIVE_CALLS = {"Confirmed", "Probable", "Detected"}
ROUTINE_PC_TYPES = {"PC_MIX8", "PC"}

detected = sorted(set(
    calls.loc[calls["call"].isin(POSITIVE_CALLS), "genus"].astype(str)
))

out = pd.DataFrame([{
    "sample_id":              sample_id,
    "sample_type":            sample_type,
    "taxa_present":           ",".join(detected),
    "n_taxa_present":         len(detected),
    "routine_positive_control": (
        "true" if sample_type.upper() in ROUTINE_PC_TYPES else "false"
    ),
}])

out.to_csv("~{sample_id}.routine_summary.tsv", sep="\t", index=False)
PY
  >>>

  output {
    File routine_summary = "~{sample_id}.routine_summary.tsv"
  }

  runtime {
    docker: docker_image
    memory: "4G"
    disks:  "local-disk 20 HDD"
  }
}

# ---------------------------------------------------------------------------
# BuildNTCBackground
#   Gather task — called once after the Phase 1 scatter.
#   Receives ONLY the NTC/NC sample files (pre-filtered in the scatter via
#   a conditional declaration + select_all), so no sample_type array needed.
#   Computes per-genus max NTC reads from alignment and centrifuge sources.
# ---------------------------------------------------------------------------
task BuildNTCBackground {
  input {
    Array[File] ntc_align_metrics     # align_metrics.tsv for each NTC/NC sample
    Array[File] ntc_cfr_genus_counts  # genus_counts.tsv  for each NTC/NC sample
    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 /opt/afi/scripts/build_ntc_background.py \
    --align-metrics-file ~{write_lines(ntc_align_metrics)} \
    --cfr-genus-file ~{write_lines(ntc_cfr_genus_counts)} \
    --out ntc_background.tsv
  >>>

  output {
    File ntc_background = "ntc_background.tsv"
  }

  runtime {
    docker: docker_image
    memory: "4G"
    disks:  "local-disk 20 HDD"
  }
}

# ---------------------------------------------------------------------------
# BuildRunSummary
#   Final gather task — produces a single run_summary.tsv that lists every
#   sample and annotates all rows with the run's PC8 validity status.
#
#   sample_id and sample_type are derived from the file contents:
#     calls_files        — sample column gives sample_id
#     validation/routine — sample_id + sample_type columns
#   This avoids WDL 1.0 incompatible splat syntax (samples[*].field).
# ---------------------------------------------------------------------------
task BuildRunSummary {
  input {
    Array[String]  run_ids              # per-sample run IDs (parallel to calls_files)
    Array[File]    calls_files
    Array[File?]   validation_summaries
    Array[File?]   routine_summaries
    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 - <<'PY'
import csv
import sys

def read_tsv(path: str) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows

def load_lines(path: str) -> list[str]:
    with open(path, encoding="utf-8") as fh:
        return [ln.rstrip("\n") for ln in fh if ln.strip()]

calls_paths = load_lines("~{write_lines(calls_files)}")
val_paths   = load_lines("~{write_lines(select_all(validation_summaries))}")
rout_paths  = load_lines("~{write_lines(select_all(routine_summaries))}")

POSITIVE_CALLS = {"Confirmed", "Probable", "Detected"}

# Build lookup dicts keyed by sample_id from summaries
# (gives us sample_type, validation_result, pc8_pass)
summary_by_id: dict[str, dict] = {}
for p in val_paths + rout_paths:
    for row in read_tsv(p):
        sid = row.get("sample_id", "").strip()
        if sid:
            summary_by_id[sid] = row

# Associate each calls file with its run_id (parallel arrays)
run_ids_list = load_lines("~{write_lines(run_ids)}")

# Map sample_id → run_id from the parallel arrays
sample_run_map: dict[str, str] = {}
for run_id_val, calls_path in zip(run_ids_list, calls_paths):
    calls_peek = read_tsv(calls_path)
    sid = calls_peek[0]["sample"] if calls_peek else ""
    if sid:
        sample_run_map[sid] = run_id_val

# Determine PC8 validity per run_id
run_pc8_valid: dict[str, str] = {}
for sid, row in summary_by_id.items():
    if row.get("sample_type", "").upper() == "PC_MIX8":
        rid = sample_run_map.get(sid, "unknown")
        run_pc8_valid[rid] = row.get("pc8_pass", "not_applicable")
# Default for run_ids that had no PC_MIX8
for rid in set(run_ids_list):
    run_pc8_valid.setdefault(rid, "no_pc8_in_run")

out_rows = []
for calls_path in calls_paths:
    calls = read_tsv(calls_path)
    sid = calls[0]["sample"] if calls else ""
    rid = sample_run_map.get(sid, run_ids_list[0] if run_ids_list else "")

    # Collect reads for positive-call genera, sort by reads descending
    detected_reads: list[tuple[str, int]] = []
    for r in calls:
        if r.get("call") in POSITIVE_CALLS:
            try:
                reads = int(r.get("reads", 0) or 0)
            except (ValueError, TypeError):
                reads = 0
            detected_reads.append((r["genus"], reads))
    detected_reads.sort(key=lambda x: x[1], reverse=True)

    total_reads = sum(rd for _, rd in detected_reads)
    if total_reads > 0:
        detected_taxa_str = ",".join(
            f"{g} ({rd/total_reads*100:.1f}%)" for g, rd in detected_reads
        )
    else:
        detected_taxa_str = ",".join(g for g, _ in detected_reads)

    summary = summary_by_id.get(sid, {})
    out_rows.append({
        "run_id":            rid,
        "sample_id":         sid,
        "sample_type":       summary.get("sample_type", ""),
        "detected_taxa":     detected_taxa_str,
        "n_detected":        len(detected_reads),
        "validation_result": summary.get("validation_result", ""),
        "pc8_pass":          summary.get("pc8_pass", ""),
        "run_pc8_valid":     run_pc8_valid.get(rid, "no_pc8_in_run"),
    })

fieldnames = [
    "run_id", "sample_id", "sample_type", "detected_taxa", "n_detected",
    "validation_result", "pc8_pass", "run_pc8_valid",
]
with open("run_summary.tsv", "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(out_rows)

n_runs = len(set(run_ids_list))
print(f"run_summary.tsv: {len(out_rows)} samples across {n_runs} run(s)", file=sys.stderr)
PY
  >>>

  output {
    File run_summary = "run_summary.tsv"
  }

  runtime {
    docker: docker_image
    memory: "8G"
    disks:  "local-disk 20 HDD"
  }
}
