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
    String docker_image = "phemarajata614/afi-terra:0.4.0"
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

out.to_csv("validation_summary.tsv", sep="\t", index=False)
PY
  >>>

  output {
    File   validation_summary = "validation_summary.tsv"
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
    String docker_image = "phemarajata614/afi-terra:0.4.0"
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

out.to_csv("routine_summary.tsv", sep="\t", index=False)
PY
  >>>

  output {
    File routine_summary = "routine_summary.tsv"
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
#   Accepts parallel arrays (sample_types, align_metrics, cfr_genus_counts),
#   filters NTC/NC samples, computes per-genus max NTC reads from each source,
#   and writes ntc_background.tsv for use in the Phase 2 scatter.
# ---------------------------------------------------------------------------
task BuildNTCBackground {
  input {
    Array[String] sample_types
    Array[File]   align_metrics
    Array[File]   cfr_genus_counts
    String docker_image = "phemarajata614/afi-terra:0.4.0"
  }

  command <<<
  python3 /opt/afi/scripts/build_ntc_background.py \
    --sample-types-file ~{write_lines(sample_types)} \
    --align-metrics-file ~{write_lines(align_metrics)} \
    --cfr-genus-file ~{write_lines(cfr_genus_counts)} \
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
#   Inputs:
#     sample_ids       — parallel arrays from Phase 2 scatter
#     sample_types
#     calls_files      — each sample's calls.tsv
#     validation_files — optional; present only for validation-mode samples
#     routine_files    — optional; present only for routine-mode samples
# ---------------------------------------------------------------------------
task BuildRunSummary {
  input {
    Array[String]  sample_ids
    Array[String]  sample_types
    Array[File]    calls_files
    Array[File?]   validation_summaries
    Array[File?]   routine_summaries
    String docker_image = "phemarajata614/afi-terra:0.4.0"
  }

  command <<<
  python3 - <<'PY'
import csv
import re
import sys
from pathlib import Path

# ---- helper: read TSV safely -----------------------------------------------
def read_tsv(path: str) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows

# ---- load parallel arrays written by WDL write_lines() --------------------
def load_lines(path: str) -> list[str]:
    with open(path, encoding="utf-8") as fh:
        return [l.rstrip("\n") for l in fh if l.strip()]

sample_ids   = load_lines("~{write_lines(sample_ids)}")
sample_types = load_lines("~{write_lines(sample_types)}")
calls_paths  = load_lines("~{write_lines(calls_files)}")
val_paths    = load_lines("~{write_lines(select_all(validation_summaries))}")
rout_paths   = load_lines("~{write_lines(select_all(routine_summaries))}")

POSITIVE_CALLS = {"Confirmed", "Probable", "Detected"}

# Index validation summaries by sample_id
val_by_id: dict[str, dict] = {}
for p in val_paths:
    for row in read_tsv(p):
        sid = row.get("sample_id", "").strip()
        if sid:
            val_by_id[sid] = row

rout_by_id: dict[str, dict] = {}
for p in rout_paths:
    for row in read_tsv(p):
        sid = row.get("sample_id", "").strip()
        if sid:
            rout_by_id[sid] = row

# Determine run-level PC8 validity
run_pc8_valid = "no_pc8_in_run"
for sid, row in val_by_id.items():
    if row.get("sample_type", "").upper() == "PC_MIX8":
        pc8_pass = row.get("pc8_pass", "not_applicable")
        run_pc8_valid = pc8_pass  # "true" or "false"
        break

out_rows = []
for sid, stype, calls_path in zip(sample_ids, sample_types, calls_paths):
    calls = read_tsv(calls_path)
    detected = sorted({r["genus"] for r in calls if r.get("call") in POSITIVE_CALLS})

    # Pull per-sample details
    val_row  = val_by_id.get(sid, {})
    rout_row = rout_by_id.get(sid, {})
    validation_result = val_row.get("validation_result", "")
    pc8_pass_this     = val_row.get("pc8_pass", "")

    out_rows.append({
        "sample_id":          sid,
        "sample_type":        stype,
        "detected_taxa":      ",".join(detected),
        "n_detected":         len(detected),
        "validation_result":  validation_result,
        "pc8_pass":           pc8_pass_this,
        "run_pc8_valid":      run_pc8_valid,
    })

fieldnames = [
    "sample_id", "sample_type", "detected_taxa", "n_detected",
    "validation_result", "pc8_pass", "run_pc8_valid",
]
with open("run_summary.tsv", "w", newline="", encoding="utf-8") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(out_rows)

print(f"run_summary.tsv: {len(out_rows)} samples, run_pc8_valid={run_pc8_valid}", file=sys.stderr)
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
