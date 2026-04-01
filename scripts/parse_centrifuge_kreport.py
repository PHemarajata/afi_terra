#!/usr/bin/env python3
"""
Parse a Centrifuge/Kraken2-style kreport and extract genus-level read counts.

The kreport format (tab-separated, no header):
  col 0 : % of reads covered by clade
  col 1 : # reads covered (clade, including children)
  col 2 : # reads assigned directly
  col 3 : rank code  (G = genus, S = species, U = unclassified, ...)
  col 4 : NCBI taxon ID
  col 5 : scientific name (indented with spaces)

Only rows with rank == "G" are kept.  Genus name is stripped of leading
whitespace.  The clade count (col 1) is used as the read total for that
genus so that reads assigned to child species are included.

Output TSV: genus <TAB> reads
"""
import argparse
import csv
import sys


def parse_kreport(path: str) -> list[dict]:
    """Return list of {genus, reads} for rank-G rows, sorted by reads desc."""
    rows = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                continue
            rank = parts[3].strip()
            if rank != "G":
                continue
            try:
                reads = int(parts[1])
            except ValueError:
                reads = 0
            genus = parts[5].strip()
            if not genus:
                continue
            rows.append({"genus": genus, "reads": reads})
    # Deduplicate: if same genus name appears more than once (unusual), sum
    merged: dict[str, int] = {}
    for r in rows:
        merged[r["genus"]] = merged.get(r["genus"], 0) + r["reads"]
    return [{"genus": g, "reads": c} for g, c in sorted(merged.items())]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract genus-level read counts from a centrifuge kreport."
    )
    parser.add_argument("--kreport", required=True, help="Path to centrifuge kreport TSV")
    parser.add_argument("--out", required=True, help="Output TSV path (genus, reads)")
    args = parser.parse_args()

    rows = parse_kreport(args.kreport)

    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["genus", "reads"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} genus rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
