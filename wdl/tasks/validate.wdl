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

# Concordance is genus-level: the pipeline detects at genus resolution, so we
# extract the first word from each expected taxon string before comparing.
# Examples:
#   "Orientia tsutsugamushi" → "Orientia"
#   "Leptospira spp"        → "Leptospira"
#   "Rickettsia"            → "Rickettsia"  (already genus)
#   "Burkholderia pseudomallei" → "Burkholderia"
# dict.fromkeys preserves order while deduplicating (handles MIXED4 where
# Streptococcus pneumoniae + Streptococcus suis → one "Streptococcus" entry).
def to_genus(t: str) -> str:
    return t.strip().split()[0] if t.strip() else ""

expected_genera  = list(dict.fromkeys(to_genus(t) for t in expected_list if to_genus(t)))
detected_lower   = {g.lower() for g in detected}
expected_lower   = {g.lower() for g in expected_genera}

detected_exp   = sorted(g for g in expected_genera if g.lower() in detected_lower)
missing_exp    = sorted(g for g in expected_genera if g.lower() not in detected_lower)
unexpected_det = sorted(g for g in detected        if g.lower() not in expected_lower)
tp_rows        = len(detected_exp)

# ---------------------------------------------------------------------------
# Order-level Rickettsiales rescue
# ---------------------------------------------------------------------------
# The V1–V2 16S region cannot reliably resolve Orientia from Rickettsia at
# genus level.  If an expected genus is Rickettsiales (Orientia or Rickettsia)
# and ANY Rickettsiales genus has alignment evidence meeting the breadth/reads
# thresholds — regardless of NTC comparison — count as Concordant at order
# level.  This ports the "any_align" rescue from afi_validate_modular.py.
#
# align_confirmed in final_calls.tsv = reads >= confirm_reads AND
# breadth >= confirm_breadth WITHOUT the NTC fold gate (computed by call_taxa.py).
RICK_GENERA_SET = {"orientia", "rickettsia"}

rick_aln_rows = calls[
    (calls["source"] == "alignment") &
    (calls["genus"].str.lower().isin(RICK_GENERA_SET))
]
any_rick_align_confirmed = any(
    str(v).lower() == "true"
    for v in rick_aln_rows.get("align_confirmed", pd.Series(dtype=str))
)

# Broader Rickettsiales order set for centrifuge rescue
RICK_ORDER_SET = {
    "orientia", "rickettsia", "anaplasma", "ehrlichia",
    "neorickettsia", "neoehrlichia", "wolbachia"
}
CFR_FLOOR = 500

# Pass 1: alignment rescue (existing logic)
order_rescued = []
truly_missing = []
for g in missing_exp:
    if g.lower() in RICK_GENERA_SET and any_rick_align_confirmed:
        order_rescued.append(g)
    else:
        truly_missing.append(g)

# Pass 2: centrifuge rescue
# Sub-case A: cfr_reads on alignment rows covers Orientia/Rickettsia
# (centrifuge rows for these genera are suppressed in call_taxa.py)
any_rick_cfr = any(
    int(float(str(r.get("cfr_reads", 0) or 0))) >= CFR_FLOOR
    for _, r in rick_aln_rows.iterrows()
)
# Sub-case B: centrifuge Detected rows for broader Rickettsiales genera
if not any_rick_cfr:
    rick_cfr_rows = calls[
        (calls["source"] == "centrifuge") &
        (calls["genus"].str.lower().isin(RICK_ORDER_SET)) &
        (calls["call"] == "Detected")
    ]
    any_rick_cfr = not rick_cfr_rows.empty

cfr_rescued   = []
still_missing = []
for g in truly_missing:
    if g.lower() in RICK_ORDER_SET and any_rick_cfr:
        cfr_rescued.append(g)
    else:
        still_missing.append(g)

# Build rescue_mechanisms string (empty when no rescue was needed)
rescue_parts = []
if order_rescued:
    rescue_parts.append(f"Rickettsiales_order({','.join(order_rescued)})")
if cfr_rescued:
    rescue_parts.append(f"Rickettsiales_cfr_rescue({','.join(cfr_rescued)})")
rescue_mechanisms = ";".join(rescue_parts)

if sample_type.upper() in VALIDATION_TYPES and expected_list:
    validation_result = "Concordant" if not still_missing else "Discordant"
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
    "missing_expected_taxa":   ",".join(still_missing),
    "order_rescue_taxa":       ",".join(order_rescued + cfr_rescued),
    "unexpected_detected_taxa":",".join(unexpected_det),
    "validation_result":       validation_result,
    "rescue_mechanisms":       rescue_mechanisms,
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
#   Receives ONLY the NTC sample files (pre-filtered in the scatter via
#   a conditional declaration + select_all), so no sample_type array needed.
#   NC (buffer/extraction negative controls) are excluded from background
#   computation — they are processed through the pipeline but only NTC samples
#   define the per-run background thresholds.
#   Computes per-run, per-genus max NTC reads from alignment + centrifuge.
#
#   ntc_run_ids is parallel to ntc_align_metrics / ntc_cfr_genus_counts.
#   One ntc_background TSV is written per distinct run_id so that samples
#   from different runs are never cross-contaminated by each other's NTCs.
# ---------------------------------------------------------------------------
task BuildNTCBackground {
  input {
    Array[String] ntc_run_ids          # run_id for each NTC sample (parallel to below)
    Array[File]   ntc_align_metrics    # align_metrics.tsv for each NTC sample
    Array[File]   ntc_cfr_genus_counts # genus_counts.tsv for each NTC sample
    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 /opt/afi/scripts/build_ntc_background.py \
    --align-metrics-file ~{write_lines(ntc_align_metrics)} \
    --cfr-genus-file ~{write_lines(ntc_cfr_genus_counts)} \
    --run-ids-file ~{write_lines(ntc_run_ids)} \
    --out-dir per_run_backgrounds
  >>>

  output {
    Array[File]   per_run_backgrounds = glob("per_run_backgrounds/ntc_background_*.tsv")
    Array[String] per_run_ids         = read_lines("run_ids.txt")
  }

  runtime {
    docker: docker_image
    memory: "4G"
    disks:  "local-disk 20 HDD"
  }
}

# ---------------------------------------------------------------------------
# MatchNTCBackground
#   Maps each sample in the batch to the NTC background for its run_id.
#   Returns a per-sample Array[File] in the same order as sample_ids so
#   Phase 2 can index directly: per_sample_backgrounds[i].
# ---------------------------------------------------------------------------
task MatchNTCBackground {
  input {
    Array[String] all_sample_run_ids   # one per sample, same order as sample_ids
    Array[String] per_run_ids          # from BuildNTCBackground
    Array[File]   per_run_backgrounds  # parallel to per_run_ids
    String docker_image = "phemarajata614/afi-terra:0.4.1"
  }

  command <<<
  python3 - <<'PY'
import shutil, sys

def load_lines(path):
    with open(path, encoding="utf-8") as f:
        return [l.rstrip("\n") for l in f if l.strip()]

per_run_ids  = load_lines("~{write_lines(per_run_ids)}")
per_run_bgs  = load_lines("~{write_lines(per_run_backgrounds)}")
sample_rids  = load_lines("~{write_lines(all_sample_run_ids)}")

bg_map = dict(zip(per_run_ids, per_run_bgs))
missing = sorted({rid for rid in sample_rids if rid not in bg_map})
if missing:
    print(
        f"ERROR: no NTC background computed for run_id(s): {missing}. "
        "Ensure every run_id in the sample table has at least one NTC sample "
        "(sample_type='NTC'). NC samples are buffer/extraction controls and do "
        "not contribute to background computation.",
        file=sys.stderr,
    )
    sys.exit(1)

output_files = []
for i, rid in enumerate(sample_rids):
    dest = f"ntc_bg_{i:04d}.tsv"
    shutil.copy(bg_map[rid], dest)
    output_files.append(dest)

with open("per_sample_backgrounds.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(output_files) + "\n")

print(
    f"Matched {len(output_files)} samples to NTC backgrounds "
    f"across {len(bg_map)} run_id(s): {list(bg_map)}",
    file=sys.stderr,
)
PY
  >>>

  output {
    Array[File] per_sample_backgrounds = glob("ntc_bg_*.tsv")
  }

  runtime {
    docker: docker_image
    memory: "2G"
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
        "expected_taxa":     summary.get("expected_taxon", ""),
        "detected_taxa":     detected_taxa_str,
        "n_detected":        len(detected_reads),
        "validation_result": summary.get("validation_result", ""),
        "rescue_mechanisms": summary.get("rescue_mechanisms", ""),
        "pc8_pass":          summary.get("pc8_pass", ""),
        "run_pc8_valid":     run_pc8_valid.get(rid, "no_pc8_in_run"),
    })

fieldnames = [
    "run_id", "sample_id", "sample_type",
    "expected_taxa", "detected_taxa", "n_detected",
    "validation_result", "rescue_mechanisms", "pc8_pass", "run_pc8_valid",
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
