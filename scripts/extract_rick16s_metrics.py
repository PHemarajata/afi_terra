#!/usr/bin/env python3
"""
Extract per-genus alignment metrics from a minimap2 BAM file against the
rickettsiales 16S panel.

Reference-to-genus mapping (same logic as afi_validate_modular.py):
  AM494475.*, AP008981.*  -> Orientia
  NC_006142.*, CP004888.*, NC_009882.* -> Rickettsia
  everything else         -> Other (skipped)

Outputs TSV: genus  mapped_reads  max_breadth
"""
import argparse
import subprocess
import pandas as pd


def ref_to_genus(ref: str) -> str:
    if ref.startswith("AM494475") or ref.startswith("AP008981"):
        return "Orientia"
    if ref.startswith("NC_006142") or ref.startswith("CP004888") or ref.startswith("NC_009882"):
        return "Rickettsia"
    return "Other"


def get_idxstats(bam_path: str) -> list:
    """Return list of {ref, ref_len, mapped_reads} via samtools idxstats."""
    result = subprocess.run(
        ["samtools", "idxstats", bam_path],
        capture_output=True, text=True, check=True
    )
    rows = []
    for line in result.stdout.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) < 4 or parts[0] == "*":
            continue
        rows.append({
            "ref": parts[0],
            "ref_len": int(parts[1]),
            "mapped_reads": int(parts[2]),
        })
    return rows


def get_depth_coverage(bam_path: str) -> dict:
    """Return {ref: num_covered_positions} using samtools depth without -r.

    Avoids samtools misinterpreting reference names that contain ':' as
    region specifications (e.g. 'AM494475.1:1322610-1324109').
    """
    result = subprocess.run(
        ["samtools", "depth", bam_path],
        capture_output=True, text=True, check=True
    )
    coverage: dict = {}
    for line in result.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        ref = parts[0]
        coverage[ref] = coverage.get(ref, 0) + 1
    return coverage


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--panel", required=True)  # retained for interface compat
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    stats = get_idxstats(args.bam)
    # Run samtools depth once for the whole BAM (panel is small; avoids
    # colon-in-refname issue when passing -r to samtools depth).
    coverage = get_depth_coverage(args.bam)

    rows = []
    for s in stats:
        genus = ref_to_genus(s["ref"])
        if genus == "Other":
            continue
        breadth = 0.0
        if s["mapped_reads"] > 0 and s["ref_len"] > 0:
            covered = coverage.get(s["ref"], 0)
            breadth = covered / s["ref_len"]
        rows.append({
            "ref": s["ref"],
            "genus": genus,
            "mapped_reads": s["mapped_reads"],
            "ref_len": s["ref_len"],
            "breadth": breadth,
        })

    REQUIRED_GENERA = ["Orientia", "Rickettsia"]

    if rows:
        df = pd.DataFrame(rows)
        agg = df.groupby("genus").agg(
            mapped_reads=("mapped_reads", "sum"),
            max_breadth=("breadth", "max"),
        ).reset_index()
    else:
        agg = pd.DataFrame(columns=["genus", "mapped_reads", "max_breadth"])

    # Ensure both genera always appear (0 reads if absent)
    present = set(agg["genus"].tolist())
    extras = [
        {"genus": g, "mapped_reads": 0, "max_breadth": 0.0}
        for g in REQUIRED_GENERA if g not in present
    ]
    if extras:
        agg = pd.concat([agg, pd.DataFrame(extras)], ignore_index=True)

    agg[["genus", "mapped_reads", "max_breadth"]].to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
