#!/usr/bin/env python3
"""
Build a run-specific NTC background file from NTC/NC sample outputs.

Called from the BuildNTCBackground WDL task with three input files (one
path per line, written by WDL write_lines()):

  --sample-types-file   : one sample_type per line (parallel to metrics/kreport arrays)
  --align-metrics-file  : one path to align_metrics.tsv per line
  --cfr-genus-file      : one path to centrifuge genus_counts.tsv per line

NTC/NC rows are identified by sample_type (case-insensitive NTC or NC).

Output TSV columns:
  genus  align_ntc_reads  cfr_ntc_reads

  align_ntc_reads : max mapped_reads across NTC samples from 16S alignment
                    (meaningful for Orientia/Rickettsia only; 0 for others)
  cfr_ntc_reads   : max reads across NTC samples from centrifuge kreport
                    (covers all genera)

If no NTC samples are present, all values are 0 (no background) — samples
will still be called but thresholds reduce to the floor alone.
"""
import argparse
import csv
import sys
from pathlib import Path


def read_lines(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as fh:
        return [line.rstrip("\n") for line in fh if line.strip()]


def load_align_metrics(path: str) -> dict[str, int]:
    """Return {genus: mapped_reads} from an align_metrics.tsv file."""
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
    """Return {genus: reads} from a centrifuge genus_counts.tsv file."""
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


NTC_TYPES = {"NTC", "NC"}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build run-specific NTC background from WDL array inputs."
    )
    parser.add_argument("--sample-types-file", required=True,
                        help="File with one sample_type per line (WDL write_lines output)")
    parser.add_argument("--align-metrics-file", required=True,
                        help="File with one align_metrics.tsv path per line")
    parser.add_argument("--cfr-genus-file", required=True,
                        help="File with one centrifuge genus_counts.tsv path per line")
    parser.add_argument("--out", required=True, help="Output ntc_background.tsv")
    args = parser.parse_args()

    sample_types = read_lines(args.sample_types_file)
    align_paths  = read_lines(args.align_metrics_file)
    cfr_paths    = read_lines(args.cfr_genus_file)

    n = len(sample_types)
    if len(align_paths) != n or len(cfr_paths) != n:
        print(
            f"ERROR: array length mismatch — sample_types={n}, "
            f"align_metrics={len(align_paths)}, cfr_genus={len(cfr_paths)}",
            file=sys.stderr,
        )
        sys.exit(1)

    align_max: dict[str, int] = {}
    cfr_max:   dict[str, int] = {}
    ntc_count = 0

    for sample_type, align_path, cfr_path in zip(sample_types, align_paths, cfr_paths):
        if sample_type.strip().upper() not in NTC_TYPES:
            continue
        ntc_count += 1

        for genus, reads in load_align_metrics(align_path).items():
            align_max[genus] = max(align_max.get(genus, 0), reads)

        for genus, reads in load_cfr_genus(cfr_path).items():
            cfr_max[genus] = max(cfr_max.get(genus, 0), reads)

    print(
        f"Found {ntc_count} NTC/NC sample(s); "
        f"{len(align_max)} align genera, {len(cfr_max)} centrifuge genera.",
        file=sys.stderr,
    )

    all_genera = sorted(set(list(align_max.keys()) + list(cfr_max.keys())))

    # Always include Orientia and Rickettsia (even if 0) for alignment thresholding
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
                "genus": genus,
                "align_ntc_reads": align_max.get(genus, 0),
                "cfr_ntc_reads":   cfr_max.get(genus, 0),
            })

    print(f"Wrote {len(all_genera)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
