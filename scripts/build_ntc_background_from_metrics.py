#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build run-specific ntc_background.tsv from one or more NTC metrics TSV files."
    )
    parser.add_argument(
        "--ntc-metrics",
        action="append",
        required=True,
        help="Path to an NTC metrics TSV (repeatable). Expected columns: genus,mapped_reads",
    )
    parser.add_argument("--out", required=True, help="Output ntc_background.tsv path")
    args = parser.parse_args()

    max_reads_by_genus: dict[str, int] = {}

    for path_text in args.ntc_metrics:
        path = Path(path_text)
        if not path.exists():
            raise SystemExit(f"Missing NTC metrics file: {path}")

        with path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise SystemExit(f"Empty metrics file: {path}")

            required = {"genus", "mapped_reads"}
            missing = required - set(reader.fieldnames)
            if missing:
                raise SystemExit(f"{path} missing required columns: {', '.join(sorted(missing))}")

            for row in reader:
                genus = (row.get("genus") or "").strip()
                if not genus:
                    continue

                raw_reads = (row.get("mapped_reads") or "").strip()
                try:
                    mapped_reads = int(float(raw_reads)) if raw_reads else 0
                except ValueError:
                    mapped_reads = 0

                previous = max_reads_by_genus.get(genus, 0)
                if mapped_reads > previous:
                    max_reads_by_genus[genus] = mapped_reads

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genus", "mapped_reads"], delimiter="\t")
        writer.writeheader()
        for genus in sorted(max_reads_by_genus):
            writer.writerow({"genus": genus, "mapped_reads": max_reads_by_genus[genus]})

    print(f"Wrote {out_path} with {len(max_reads_by_genus)} taxa")


if __name__ == "__main__":
    main()
