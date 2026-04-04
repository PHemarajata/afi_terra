#!/usr/bin/env python3
"""
Build NTC background file(s) from NTC/NC sample metrics.

Called from the BuildNTCBackground WDL task.  All files passed in are
already NTC/NC — filtering happens in the WDL scatter via a conditional
declaration + select_all, so no sample_type logic is needed here.

Input files (one path per line, written by WDL write_lines()):
  --align-metrics-file  : paths to align_metrics.tsv files (NTC samples only)
  --cfr-genus-file      : paths to centrifuge genus_counts.tsv files (NTC samples only)
  --run-ids-file        : run_id for each NTC sample (parallel to the above;
                          required for multi-run batches)

Output (per-run mode — always used when --run-ids-file is given):
  per_run_backgrounds/ntc_background_{run_id}.tsv  — one file per distinct run_id
  run_ids.txt                                       — run_ids, one per line
  backgrounds.txt                                   — matching file paths, one per line
  (Both manifests are read by WDL read_lines() to produce Array outputs.)

Output TSV columns (all modes):
  genus  align_ntc_reads  cfr_ntc_reads
"""
import argparse
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path


REQUIRED_GENERA = ("Orientia", "Rickettsia")


def read_lines(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as fh:
        return [line.rstrip("\n") for line in fh if line.strip()]


def load_align_metrics(path: str) -> dict[str, int]:
    out: dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genus = (row.get("genus") or "").strip()
            if not genus:
                continue
            try:
                reads = int(float(row.get("mapped_reads", 0) or 0))
            except (ValueError, TypeError):
                reads = 0
            out[genus] = max(out.get(genus, 0), reads)
    return out


def load_cfr_genus(path: str) -> dict[str, int]:
    out: dict[str, int] = {}
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            genus = (row.get("genus") or "").strip()
            if not genus:
                continue
            try:
                reads = int(float(row.get("reads", 0) or 0))
            except (ValueError, TypeError):
                reads = 0
            out[genus] = max(out.get(genus, 0), reads)
    return out


def compute_background(
    align_paths: list[str],
    cfr_paths: list[str],
) -> tuple[dict[str, int], dict[str, int]]:
    """Return (align_max, cfr_max) dicts for the given NTC file lists."""
    align_max: dict[str, int] = {}
    for p in align_paths:
        for genus, reads in load_align_metrics(p).items():
            align_max[genus] = max(align_max.get(genus, 0), reads)

    cfr_max: dict[str, int] = {}
    for p in cfr_paths:
        for genus, reads in load_cfr_genus(p).items():
            cfr_max[genus] = max(cfr_max.get(genus, 0), reads)

    return align_max, cfr_max


def write_background(
    align_max: dict[str, int],
    cfr_max: dict[str, int],
    out_path: str,
) -> int:
    """Write ntc_background.tsv; return number of genera written."""
    all_genera = sorted(set(list(align_max.keys()) + list(cfr_max.keys())))
    for g in REQUIRED_GENERA:
        if g not in all_genera:
            all_genera = sorted(all_genera + [g])

    p = Path(out_path)
    p.parent.mkdir(parents=True, exist_ok=True)

    with p.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["genus", "align_ntc_reads", "cfr_ntc_reads"],
            delimiter="\t",
        )
        writer.writeheader()
        for genus in all_genera:
            writer.writerow({
                "genus":           genus,
                "align_ntc_reads": align_max.get(genus, 0),
                "cfr_ntc_reads":   cfr_max.get(genus, 0),
            })
    return len(all_genera)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate NTC metrics into per-run ntc_background files."
    )
    parser.add_argument("--align-metrics-file", required=True,
                        help="File listing NTC align_metrics.tsv paths (one per line)")
    parser.add_argument("--cfr-genus-file", required=True,
                        help="File listing NTC centrifuge genus_counts.tsv paths (one per line)")
    parser.add_argument("--run-ids-file", required=True,
                        help="File with the run_id for each NTC sample (parallel to the "
                             "above two lists, one per line)")
    parser.add_argument("--out-dir", default="per_run_backgrounds",
                        help="Directory for per-run background TSV files "
                             "(default: per_run_backgrounds/)")
    args = parser.parse_args()

    align_paths = read_lines(args.align_metrics_file)
    cfr_paths   = read_lines(args.cfr_genus_file)
    run_ids     = read_lines(args.run_ids_file)

    if not (len(align_paths) == len(cfr_paths) == len(run_ids)):
        print(
            f"ERROR: mismatched list lengths — "
            f"align={len(align_paths)}, cfr={len(cfr_paths)}, run_ids={len(run_ids)}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Group NTC file paths by run_id (preserve insertion order of first seen)
    groups: dict[str, dict[str, list]] = {}
    for rid, ap, cp in zip(run_ids, align_paths, cfr_paths):
        if rid not in groups:
            groups[rid] = {"align": [], "cfr": []}
        groups[rid]["align"].append(ap)
        groups[rid]["cfr"].append(cp)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ordered_run_ids: list[str] = []
    bg_paths: list[str] = []

    for rid in sorted(groups.keys()):
        files = groups[rid]
        align_max, cfr_max = compute_background(files["align"], files["cfr"])
        out_path = os.path.abspath(str(out_dir / f"ntc_background_{rid}.tsv"))
        n = write_background(align_max, cfr_max, out_path)
        ordered_run_ids.append(rid)
        bg_paths.append(out_path)
        print(
            f"run_id={rid}: {len(files['align'])} NTC sample(s), "
            f"{n} genera → {out_path}",
            file=sys.stderr,
        )

    # Write WDL-readable manifests (used with read_lines() in the task output block)
    with open("run_ids.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(ordered_run_ids) + "\n")
    with open("backgrounds.txt", "w", encoding="utf-8") as fh:
        fh.write("\n".join(bg_paths) + "\n")

    print(
        f"Wrote {len(ordered_run_ids)} per-run background file(s) for "
        f"run_ids: {', '.join(ordered_run_ids)}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
