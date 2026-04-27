#!/usr/bin/env python3
"""
NTC-aware taxa calling combining two evidence sources:

  Module 3 (alignment) — Orientia and Rickettsia only
    Source: 16S minimap2 alignment metrics (mapped_reads, max_breadth)
    Tiers : Confirmed / Probable / Not_Confirmed / Negative

  Module 1 (centrifuge) — all other genera
    Source: centrifuge kreport genus counts
    Tiers : Detected / Not_Detected

Thresholds match afi_validate_modular.py defaults:
  Alignment Confirmed  : reads >= ALIGN_CONFIRM_READS AND breadth >= ALIGN_CONFIRM_BREADTH
                         AND reads >= ALIGN_FOLD * align_ntc_reads
  Alignment Probable   : reads >= 50 AND breadth >= 0.20 AND reads > align_ntc_reads
  Alignment Not_Conf.  : reads >= 50 AND reads <= align_ntc_reads
  Alignment Negative   : reads < 50
  Centrifuge Detected  : reads >= CFR_FLOOR AND reads >= CFR_FOLD * cfr_ntc_reads
  Centrifuge Not_Det.  : otherwise

NTC background file format (from build_ntc_background.py):
  genus  align_ntc_reads  cfr_ntc_reads

Output TSV columns:
  sample  genus  source  reads  breadth  ntc_reads  call  cfr_reads  rescued

  cfr_reads — centrifuge read count for the genus (alignment rows only;
              blank for centrifuge rows where `reads` IS the cfr count).
  rescued   — true when source=alignment, call is Confirmed/Probable, AND
              centrifuge alone would NOT have reached the detection threshold
              (i.e. cfr_reads < cfr_floor OR cfr_reads < cfr_fold × cfr_ntc_reads).
              false for non-positive alignment calls; blank for centrifuge rows.
"""
import argparse
import csv
import sys

RICK_GENERA = {"Orientia", "Rickettsia"}


def load_ntc_background(path: str) -> dict[str, dict[str, int]]:
    """
    Returns {genus: {align_ntc_reads: int, cfr_ntc_reads: int}}.
    Handles both the new dual-column format and the old single-column
    format (mapped_reads only) for backwards compatibility with pre-computed files.
    """
    out: dict[str, dict[str, int]] = {}
    with open(path, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        for row in reader:
            genus = (row.get("genus") or "").strip()
            if not genus:
                continue
            if "align_ntc_reads" in fieldnames:
                # New dual-column format
                try:
                    align_r = int(float(row.get("align_ntc_reads", 0) or 0))
                except (ValueError, TypeError):
                    align_r = 0
                try:
                    cfr_r = int(float(row.get("cfr_ntc_reads", 0) or 0))
                except (ValueError, TypeError):
                    cfr_r = 0
            else:
                # Legacy single-column format (mapped_reads) — used for alignment only
                try:
                    legacy = int(float(row.get("mapped_reads", 0) or 0))
                except (ValueError, TypeError):
                    legacy = 0
                align_r = legacy
                cfr_r = 0
            out[genus] = {"align_ntc_reads": align_r, "cfr_ntc_reads": cfr_r}
    return out


def load_align_metrics(path: str) -> list[dict]:
    """Return rows from align_metrics.tsv: genus, mapped_reads, max_breadth."""
    rows = []
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
            try:
                breadth = float(row.get("max_breadth", 0) or 0)
            except (ValueError, TypeError):
                breadth = 0.0
            rows.append({"genus": genus, "reads": reads, "breadth": breadth})
    return rows


def load_cfr_genus(path: str) -> list[dict]:
    """Return rows from centrifuge genus_counts.tsv: genus, reads."""
    rows = []
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
            rows.append({"genus": genus, "reads": reads})
    return rows


def call_alignment(genus: str, reads: int, breadth: float, ntc: dict,
                   confirm_reads: int, confirm_breadth: float, fold: float) -> tuple[str, int, bool]:
    """Return (call, ntc_reads, align_confirmed) for an alignment-source genus.

    align_confirmed is True when reads and breadth meet the confirmation thresholds
    WITHOUT considering the NTC fold comparison.  This NTC-free flag is used by the
    order-level Rickettsiales rescue rule in CompareExpectedConcordance.
    """
    ntc_reads = ntc.get(genus, {}).get("align_ntc_reads", 0)
    align_confirmed = (reads >= confirm_reads and breadth >= confirm_breadth)
    if align_confirmed and reads >= fold * ntc_reads:
        call = "Confirmed"
    elif reads >= 50 and breadth >= 0.20 and reads > ntc_reads:
        call = "Probable"
    elif reads >= 50 and reads <= ntc_reads:
        call = "Not_Confirmed"
    else:
        call = "Negative"
    return call, ntc_reads, align_confirmed


def call_centrifuge(genus: str, reads: int, ntc: dict,
                    floor: int, fold: float) -> tuple[str, int]:
    """Return (call, ntc_reads) for a centrifuge-source genus."""
    ntc_reads = ntc.get(genus, {}).get("cfr_ntc_reads", 0)
    if reads >= floor and reads >= fold * ntc_reads:
        call = "Detected"
    else:
        call = "Not_Detected"
    return call, ntc_reads


def main() -> None:
    parser = argparse.ArgumentParser(
        description="NTC-aware taxa calling from alignment and centrifuge evidence."
    )
    parser.add_argument("--sample",         required=True)
    parser.add_argument("--align-metrics",  required=True,
                        help="align_metrics.tsv (genus, mapped_reads, max_breadth)")
    parser.add_argument("--cfr-genus",      required=True,
                        help="centrifuge genus_counts.tsv (genus, reads)")
    parser.add_argument("--ntc",            required=True,
                        help="ntc_background.tsv (genus, align_ntc_reads, cfr_ntc_reads)")
    parser.add_argument("--out",            required=True)

    # Alignment thresholds (Module 3)
    parser.add_argument("--align-confirm-reads",   type=int,   default=100)
    parser.add_argument("--align-confirm-breadth", type=float, default=0.25)
    parser.add_argument("--align-fold",            type=float, default=5.0)

    # Centrifuge thresholds (Module 1)
    parser.add_argument("--cfr-floor", type=int,   default=500)
    parser.add_argument("--cfr-fold",  type=float, default=5.0)

    args = parser.parse_args()

    ntc          = load_ntc_background(args.ntc)
    align_rows   = load_align_metrics(args.align_metrics)
    cfr_rows     = load_cfr_genus(args.cfr_genus)

    # Build a quick lookup of centrifuge reads for every genus (used for
    # cfr_reads + rescued columns on alignment rows).
    cfr_reads_by_genus: dict[str, int] = {r["genus"]: r["reads"] for r in cfr_rows}

    results = []

    # --- Module 3: alignment-based calls for Orientia / Rickettsia ---
    align_genera_seen = set()
    for row in align_rows:
        genus = row["genus"]
        if genus not in RICK_GENERA:
            continue
        align_genera_seen.add(genus)
        call, ntc_reads, align_confirmed = call_alignment(
            genus, row["reads"], row["breadth"], ntc,
            args.align_confirm_reads, args.align_confirm_breadth, args.align_fold,
        )
        cfr_r = cfr_reads_by_genus.get(genus, 0)
        cfr_ntc = ntc.get(genus, {}).get("cfr_ntc_reads", 0)
        cfr_would_detect = (cfr_r >= args.cfr_floor and cfr_r >= args.cfr_fold * cfr_ntc)
        if call in ("Confirmed", "Probable"):
            rescued = "false" if cfr_would_detect else "true"
        else:
            rescued = "false"
        results.append({
            "sample":          args.sample,
            "genus":           genus,
            "source":          "alignment",
            "reads":           row["reads"],
            "breadth":         round(row["breadth"], 4),
            "ntc_reads":       ntc_reads,
            "call":            call,
            "cfr_reads":       cfr_r,
            "rescued":         rescued,
            "align_confirmed": str(align_confirmed).lower(),
        })

    # Ensure both Rickettsiales genera always appear (Negative if absent from BAM)
    for genus in RICK_GENERA:
        if genus not in align_genera_seen:
            cfr_r = cfr_reads_by_genus.get(genus, 0)
            results.append({
                "sample":          args.sample,
                "genus":           genus,
                "source":          "alignment",
                "reads":           0,
                "breadth":         0.0,
                "ntc_reads":       ntc.get(genus, {}).get("align_ntc_reads", 0),
                "call":            "Negative",
                "cfr_reads":       cfr_r,
                "rescued":         "false",
                "align_confirmed": "false",
            })

    # --- Module 1: centrifuge-based calls for all non-Rickettsiales genera ---
    for row in cfr_rows:
        genus = row["genus"]
        if genus in RICK_GENERA:
            # Rickettsiales genera are handled via alignment above; skip here.
            # Their centrifuge counts are captured in cfr_reads_by_genus above.
            continue
        call, ntc_reads = call_centrifuge(
            genus, row["reads"], ntc,
            args.cfr_floor, args.cfr_fold,
        )
        results.append({
            "sample":          args.sample,
            "genus":           genus,
            "source":          "centrifuge",
            "reads":           row["reads"],
            "breadth":         "",   # not applicable for centrifuge
            "ntc_reads":       ntc_reads,
            "call":            call,
            "cfr_reads":       "",   # reads IS the cfr count; no separate column needed
            "rescued":         "",   # not applicable for centrifuge-source rows
            "align_confirmed": "",   # not applicable for centrifuge-source rows
        })

    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample", "genus", "source", "reads", "breadth",
                        "ntc_reads", "call", "cfr_reads", "rescued", "align_confirmed"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(results)

    n_pos = sum(1 for r in results if r["call"] in ("Confirmed", "Probable", "Detected"))
    print(
        f"Wrote {len(results)} rows ({n_pos} positive calls) to {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
