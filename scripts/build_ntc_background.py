#!/usr/bin/env python3
"""
Build a run-specific NTC background file from NTC/NC sample metrics.

Called from the BuildNTCBackground WDL task.  All files passed in are
already NTC/NC — filtering happens in the WDL scatter via a conditional
declaration + select_all, so no sample_type logic is needed here.

Input files (one path per line, written by WDL write_lines()):
  --align-metrics-file  : paths to align_metrics.tsv files (NTC samples only)
  --cfr-genus-file      : paths to centrifuge genus_counts.tsv files (NTC samples only)

Output TSV columns:
  genus  align_ntc_reads  cfr_ntc_reads
"""
import argparse
import csv
import sys
from pathlib import Path


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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate NTC metrics into ntc_background.tsv."
    )
    parser.add_argument("--align-metrics-file", required=True,
                        help="File listing NTC align_metrics.tsv paths (one per line)")
    parser.add_argument("--cfr-genus-file", required=True,
                        help="File listing NTC centrifuge genus_counts.tsv paths (one per line)")
    parser.add_argument("--out", required=True, help="Output ntc_background.tsv")
    args = parser.parse_args()

    align_paths = read_lines(args.align_metrics_file)
    cfr_paths   = read_lines(args.cfr_genus_file)

    align_max: dict[str, int] = {}
    for path in align_paths:
        for genus, reads in load_align_metrics(path).items():
            align_max[genus] = max(align_max.get(genus, 0), reads)

    cfr_max: dict[str, int] = {}
    for path in cfr_paths:
        for genus, reads in load_cfr_genus(path).items():
            cfr_max[genus] = max(cfr_max.get(genus, 0), reads)

    print(
        f"Processed {len(align_paths)} NTC align file(s), {len(cfr_paths)} NTC cfr file(s). "
        f"{len(align_max)} align genera, {len(cfr_max)} cfr genera.",
        file=sys.stderr,
    )

    all_genera = sorted(set(list(align_max.keys()) + list(cfr_max.keys())))
    for g in ("Orientia", "Rickettsia"):
        if g not in all_genera:
            all_genera = sorted(all_genera + [g])

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", encoding="utf-8", newline="") as fh:
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

    print(f"Wrote {len(all_genera)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
