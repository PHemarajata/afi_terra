#!/usr/bin/env python3
import argparse
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd


# -----------------------------
# Utilities
# -----------------------------

def clean_expected_genus(expected_results: str) -> str:
    """
    expected_results examples:
      'Streptococcus pneumoniae'
      'Salmonella Enterica'
      'Lactobacillus fermentum*'
    Return genus only, cleaned.
    """
    if pd.isna(expected_results) or str(expected_results).strip() == "":
        return ""
    s = str(expected_results).strip()
    # remove trailing punctuation/asterisks
    s = re.sub(r"[^A-Za-z0-9\s_-]+$", "", s)
    genus = s.split()[0]
    genus = re.sub(r"[^A-Za-z0-9_-]", "", genus)
    return genus


def infer_sample_type(sample_name: str) -> str:
    """
    Deterministic type inference from sample_name:
      - NC / NTC
      - PC (mixed Zymo) like PC_S10, PC_S12, PC-20251016
      - MIXED4 (special) = 5_Mixed
      - PC_SINGLE = BS-20251016_S7 etc (single-organism controls)
      - CLINICAL = everything else
    """
    s = str(sample_name)

    if s.startswith("NTC"):
        return "NTC"
    if s.startswith("NC"):
        return "NC"
    if s == "5_Mixed":
        return "MIXED4"
    if s.startswith("PC"):
        return "PC_MIX8"
    # two-letter code with date-like tag
    if re.match(r"^[A-Z]{2}-\d{8}_S\d+$", s):
        return "PC_SINGLE"
    return "CLINICAL"


def normalize_sample_type(value: str) -> str:
    s = str(value).strip().upper()
    if "NTC" in s or s == "NC":
        return "NTC"
    if "PC_MIX8" in s:
        return "PC_MIX8"
    if "PC_MIX4" in s or "MIXED4" in s:
        return "MIXED4"
    if "PC_SINGLE" in s:
        return "PC_SINGLE"
    if "CLINICAL" in s:
        return "CLINICAL"
    return str(value).strip()

def resolve_report_path(report_dir: Path, report_name: str) -> Path:
    """
    Accept mapping names like:
      00126_S6_L001.report.txt
      00126_S6_L001_kraken2_report.txt
    and resolve to an existing path in report_dir by trying common alternates.
    """
    report_name = str(report_name).strip()
    p = report_dir / report_name
    if p.exists():
        return p

    # Try converting between naming conventions
    alts = []

    # .report.txt -> _kraken2_report.txt
    if report_name.endswith(".report.txt"):
        stem = report_name[:-len(".report.txt")]
        alts.append(f"{stem}_kraken2_report.txt")

    # _kraken2_report.txt -> .report.txt
    if report_name.endswith("_kraken2_report.txt"):
        stem = report_name[:-len("_kraken2_report.txt")]
        alts.append(f"{stem}.report.txt")

    # try also without any suffix (rare)
    alts.append(report_name.replace(".report.txt","").replace("_kraken2_report.txt","").replace("_kreport.txt",""))

    # _kreport.txt <-> _kraken2_report.txt / .report.txt
    if report_name.endswith("_kreport.txt"):
        stem = report_name[:-len("_kreport.txt")]
        alts.append(f"{stem}_kraken2_report.txt")
        alts.append(f"{stem}.report.txt")

    for a in alts:
        p2 = report_dir / a
        if p2.exists():
            return p2

    raise FileNotFoundError(f"Could not resolve report '{report_name}' in {report_dir}")

def parse_kraken_report_genus(path: Path) -> dict:
    """
    Parse Kraken2 --report output.
    Return genus -> clade_reads (col 2).
    """
    genus_counts = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                clade_reads = int(parts[1].strip())
            except ValueError:
                continue
            rank = parts[3].strip()
            name = parts[5].strip()
            if rank == "G":
                genus_counts[name] = clade_reads
    return genus_counts


def tier_from_fold(reads: int, ncmax: int, floor: int, fold: float) -> bool:
    """
    A simple boolean pass for validation:
      pass if reads >= floor AND reads >= fold * ncmax
    If ncmax == 0, pass requires reads >= floor only.
    """
    if reads < floor:
        return False
    if ncmax <= 0:
        return True
    return reads >= fold * ncmax


# -----------------------------
# Module 1: Kraken 16GB
# -----------------------------

def module_kk16g(meta: pd.DataFrame, kk16_dir: Path,
                floor: int = 500, fold: float = 5.0) -> pd.DataFrame:
    """
    For each sample, compute per-genus reads from 16GB report.
    Compute NCmax per run_id, genus using only NC/NTC within that run.
    Return long table: sample_name, run_id, genus, reads, ncmax, pass16g
    """
    # parse reports once per report filename
    cache = {}
    rows = []

    # Build run-wise NCmax
    ncmax = defaultdict(lambda: defaultdict(int))

    for _, r in meta.iterrows():
        run_id = int(r["run_id"])
        sname = r["sample_name"]
        stype = r["sample_type"]
        rep = kk16_dir / r["kk_report_name_16GB"]

        if rep.name not in cache:
            if not rep.exists():
                raise FileNotFoundError(f"Missing 16GB report: {rep}")
            cache[rep.name] = parse_kraken_report_genus(rep)

        if stype in ("NC", "NTC"):
            for g, c in cache[rep.name].items():
                if c > ncmax[run_id][g]:
                    ncmax[run_id][g] = c

    # Emit calls for all samples
    for _, r in meta.iterrows():
        run_id = int(r["run_id"])
        sname = r["sample_name"]
        rep = kk16_dir / r["kk_report_name_16GB"]
        gcounts = cache[rep.name]

        # For validation: we mainly care about expected genera, but keep all genera for later expansion.
        for g, reads in gcounts.items():
            base = ncmax[run_id].get(g, 0)
            passed = tier_from_fold(reads, base, floor=floor, fold=fold)
            rows.append({
                "run_id": run_id,
                "sample_name": sname,
                "genus": g,
                "kk16_reads": reads,
                "kk16_ncmax": base,
                "kk16_pass": passed,
            })

    return pd.DataFrame(rows)


# -----------------------------
# Module 2: Kraken rickettsiales DB
# -----------------------------

def module_kk_rick(meta: pd.DataFrame, kkr_dir: Path,
                   floor: int = 100, fold: float = 3.0) -> pd.DataFrame:
    """
    Parse rickettsiales-custom Kraken2 reports.

    IMPORTANT BEHAVIOR (by request):
    - Compute an ORDER-LEVEL aggregate signal representing "any taxon in the rickettsiales custom DB".
      This is implemented as the total genus-level reads in the report (sum across all genera in that report),
      with a per-run NCmax computed as the maximum of that total among NTCs from the same run.
    - Still emit genus-level rows for Orientia and Rickettsia for transparency/debugging.

    Output rows include:
      genus == "Rickettsiales"  -> order-level aggregate screening flag (kkr_pass)
      genus == "Orientia"       -> genus-level counts (kkr_pass also computed, but not required for rescue)
      genus == "Rickettsia"     -> genus-level counts (kkr_pass also computed, but not required for rescue)
    """
    target_genera = {"Orientia", "Rickettsia"}
    cache: dict[str, dict[str, int]] = {}
    ncmax_genus = defaultdict(lambda: defaultdict(int))   # run_id -> genus -> max reads in NTCs
    ncmax_order = defaultdict(int)                        # run_id -> max total reads in NTCs

    # --- cache reports and compute per-run NCmax (genus + order) ---
    for _, r in meta.iterrows():
        run_id = int(r["run_id"])
        stype = str(r["sample_type"])
        rep = resolve_report_path(kkr_dir, r["kk_report_name_rick"])

        if rep.name not in cache:
            if not rep.exists():
                raise FileNotFoundError(f"Missing rickettsiales report: {rep}")
            cache[rep.name] = parse_kraken_report_genus(rep)

        if stype in ("NC", "NTC"):
            gcounts = cache[rep.name]
            # genus NCmax for transparency
            for g in target_genera:
                c = int(gcounts.get(g, 0))
                if c > ncmax_genus[run_id][g]:
                    ncmax_genus[run_id][g] = c
            # ORDER NCmax = max total genus-level reads in the custom DB
            tot = int(sum(gcounts.values()))
            if tot > ncmax_order[run_id]:
                ncmax_order[run_id] = tot

    # --- emit per-sample rows ---
    rows = []
    for _, r in meta.iterrows():
        run_id = int(r["run_id"])
        sname = str(r["sample_name"])
        rep = resolve_report_path(kkr_dir, r["kk_report_name_rick"])
        gcounts = cache[rep.name]

        # ORDER row (any rickettsiales DB taxon)
        order_reads = int(sum(gcounts.values()))
        order_base = int(ncmax_order.get(run_id, 0))
        order_pass = (order_reads >= floor)  # order-level screen ignores NC because DB is rickettsiales-only
        rows.append({
            "run_id": run_id,
            "sample_name": sname,
            "genus": "Rickettsiales",
            "kkr_reads": order_reads,
            "kkr_ncmax": order_base,
            "kkr_pass": bool(order_pass),
        })

        # Genus rows (Orientia/Rickettsia)
        for g in sorted(target_genera):
            reads = int(gcounts.get(g, 0))
            base = int(ncmax_genus[run_id].get(g, 0))
            passed = tier_from_fold(reads, base, floor=floor, fold=fold)
            rows.append({
                "run_id": run_id,
                "sample_name": sname,
                "genus": g,
                "kkr_reads": reads,
                "kkr_ncmax": base,
                "kkr_pass": bool(passed),
            })

    return pd.DataFrame(rows)


# -----------------------------
# Module 3: rickettsiales 16S alignment metrics
# -----------------------------

def alignment_genus(ref: str) -> str:
    if ref.startswith("AM494475") or ref.startswith("AP008981"):
        return "Orientia"
    if ref.startswith("NC_006142") or ref.startswith("CP004888") or ref.startswith("NC_009882"):
        return "Rickettsia"
    # decoys / others ignored
    return "Other"


def module_align_rick16S(meta: pd.DataFrame, minimap_summary_tsv: Path,
                         confirmed_reads: int = 100, confirmed_breadth: float = 0.25,
                         fold_over_nc: float = 5.0,
                         debug_out_tsv: Path | None = None,
                         debug_fallback_out_tsv: Path | None = None) -> pd.DataFrame:
    """
    Use metrics/summary.tsv:
      sample ref mapped_reads ref_len breadth
    Aggregate to sample + genus:
      mapped_reads = sum
      max_breadth = max
    Call Confirmed if reads>=100 and breadth>=0.25 (V1–V3 aware)
    Then final positive if Confirmed and reads >= fold_over_nc * NCmax_align(run, genus)
    """
    aln = pd.read_csv(minimap_summary_tsv, sep="\t")
    # Drop undetermined and diluted technical samples at alignment stage too
    aln["sample"] = aln["sample"].astype(str).str.strip()
    aln = aln[~aln["sample"].str.startswith("Undetermined", na=False)].copy()
    aln = aln[~aln["sample"].isin({"003691_S11_L001","019371_S12_L001","027571_S13_L001"})].copy()
    aln["genus"] = aln["ref"].map(alignment_genus)
    aln = aln[aln["genus"].isin(["Orientia", "Rickettsia"])].copy()
    aln["breadth"] = aln["breadth"].fillna(0)

    agg = aln.groupby(["sample", "genus"]).agg(
        align_mapped_reads=("mapped_reads", "sum"),
        align_max_breadth=("breadth", "max"),
    ).reset_index()

    # Map alignment sample -> run_id/sample_name with robust fallbacks.
    # Priority:
    #   1) exact sample_name match (e.g. 00126_S6_L001)
    #   2) exact kk_report_name_16GB stem match
    #   3) stripped-lane fallback (_S\d+_L\d+ removed) ONLY when mapping is unambiguous
    meta_join = meta[["run_id", "sample_name", "kk_report_name_16GB"]].drop_duplicates().copy()
    meta_join["report_stem"] = meta_join["kk_report_name_16GB"].astype(str).str.strip().str.replace(
        r"(_kraken2_report\.txt(\.gz)?|\.report\.txt(\.gz)?|_kreport\.txt(\.gz)?)$",
        "",
        regex=True,
    )

    meta_strip_all = meta_join[["run_id", "sample_name", "report_stem"]].copy()
    meta_strip_all["sample_stripped"] = meta_strip_all["report_stem"].str.replace(r"_S\d+_L\d+$", "", regex=True)
    meta_strip_all = meta_strip_all.drop_duplicates(subset=["sample_stripped", "run_id", "sample_name", "report_stem"])

    ambiguous_map = (
        meta_strip_all.groupby("sample_stripped", as_index=False)
        .agg(
            n_candidate_run_ids=("run_id", "nunique"),
            n_candidate_sample_names=("sample_name", "nunique"),
            candidate_run_ids=("run_id", lambda s: ",".join(map(str, sorted(set(s))))),
            candidate_sample_names=("sample_name", lambda s: ",".join(sorted(set(map(str, s))))),
            candidate_report_stems=("report_stem", lambda s: ",".join(sorted(set(map(str, s))))),
        )
    )
    ambiguous_map = ambiguous_map[
        (ambiguous_map["n_candidate_run_ids"] > 1) |
        (ambiguous_map["n_candidate_sample_names"] > 1)
    ].copy()

    agg_merged = agg.copy()
    agg_merged["run_id"] = pd.NA
    agg_merged["sample_name"] = pd.NA
    agg_merged["map_stage"] = "unmapped"

    # 1) exact sample_name match
    by_sample_name = meta_join.drop_duplicates(subset=["sample_name"]).set_index("sample_name")
    by_sample_name_run = by_sample_name["run_id"]
    sample_name_identity = pd.Series(by_sample_name_run.index, index=by_sample_name_run.index)
    agg_merged["run_id"] = agg_merged["sample"].map(by_sample_name_run)
    agg_merged["sample_name"] = agg_merged["sample"].map(sample_name_identity)
    matched = agg_merged["run_id"].notna()
    agg_merged.loc[matched, "map_stage"] = "exact_sample_name"

    # 2) exact report stem match
    unmatched = agg_merged["run_id"].isna()
    if unmatched.any():
        by_report_stem = meta_join.drop_duplicates(subset=["report_stem"]).set_index("report_stem")
        agg_merged.loc[unmatched, "run_id"] = agg_merged.loc[unmatched, "sample"].map(by_report_stem["run_id"])
        agg_merged.loc[unmatched, "sample_name"] = agg_merged.loc[unmatched, "sample"].map(by_report_stem["sample_name"])
        newly_matched = unmatched & agg_merged["run_id"].notna()
        agg_merged.loc[newly_matched, "map_stage"] = "exact_report_stem"

    # 2b) strip lane suffix only (e.g. NC_S14_L001 -> NC_S14)
    unmatched = agg_merged["run_id"].isna()
    if unmatched.any():
        lane_trimmed = agg_merged.loc[unmatched, "sample"].str.replace(r"_L\d+$", "", regex=True)
        agg_merged.loc[unmatched, "run_id"] = lane_trimmed.map(by_sample_name_run)
        agg_merged.loc[unmatched, "sample_name"] = lane_trimmed.map(sample_name_identity)
        newly_matched = unmatched & agg_merged["run_id"].notna()
        agg_merged.loc[newly_matched, "map_stage"] = "lane_trim_sample_name"

        still_unmatched = agg_merged["run_id"].isna()
        if still_unmatched.any():
            by_report_stem = meta_join.drop_duplicates(subset=["report_stem"]).set_index("report_stem")
            lane_trimmed = agg_merged.loc[still_unmatched, "sample"].str.replace(r"_L\d+$", "", regex=True)
            agg_merged.loc[still_unmatched, "run_id"] = lane_trimmed.map(by_report_stem["run_id"])
            agg_merged.loc[still_unmatched, "sample_name"] = lane_trimmed.map(by_report_stem["sample_name"])
            newly_matched = still_unmatched & agg_merged["run_id"].notna()
            agg_merged.loc[newly_matched, "map_stage"] = "lane_trim_report_stem"

    # 3) lane-stripped fallback (only unique keys to avoid ambiguous one-to-many matches)
    unmatched = agg_merged["run_id"].isna()
    if unmatched.any():
        meta_strip = meta_join[["run_id", "sample_name", "report_stem"]].copy()
        meta_strip["sample_stripped"] = meta_strip["report_stem"].str.replace(r"_S\d+_L\d+$", "", regex=True)

        meta_strip_unique = meta_strip.drop_duplicates(subset=["sample_stripped", "run_id", "sample_name"])
        unique_counts = meta_strip_unique.groupby("sample_stripped")["run_id"].nunique()
        unique_keys = unique_counts[unique_counts == 1].index
        by_sample_stripped = (
            meta_strip_unique[meta_strip_unique["sample_stripped"].isin(unique_keys)]
            .drop_duplicates(subset=["sample_stripped"])
            .set_index("sample_stripped")
        )

        stripped_samples = agg_merged.loc[unmatched, "sample"].str.replace(r"_S\d+_L\d+$", "", regex=True)
        agg_merged.loc[unmatched, "run_id"] = stripped_samples.map(by_sample_stripped["run_id"])
        agg_merged.loc[unmatched, "sample_name"] = stripped_samples.map(by_sample_stripped["sample_name"])
        newly_matched = unmatched & agg_merged["run_id"].notna()
        agg_merged.loc[newly_matched, "map_stage"] = "stripped_unique"

    if debug_out_tsv is not None:
        unresolved = agg_merged[agg_merged["run_id"].isna()].copy()
        if unresolved.empty:
            unresolved_counts = pd.DataFrame(columns=["sample_stripped", "observed_unmatched_rows"])
        else:
            unresolved["sample_stripped"] = unresolved["sample"].str.replace(r"_S\d+_L\d+$", "", regex=True)
            unresolved_counts = (
                unresolved.groupby("sample_stripped", as_index=False)
                .size()
                .rename(columns={"size": "observed_unmatched_rows"})
            )

        debug_df = ambiguous_map.merge(unresolved_counts, on="sample_stripped", how="left")
        if "observed_unmatched_rows" not in debug_df.columns:
            debug_df["observed_unmatched_rows"] = 0
        debug_df["observed_unmatched_rows"] = debug_df["observed_unmatched_rows"].fillna(0).astype(int)
        debug_df = debug_df.sort_values(["observed_unmatched_rows", "sample_stripped"], ascending=[False, True])
        debug_out_tsv.parent.mkdir(parents=True, exist_ok=True)
        debug_df.to_csv(debug_out_tsv, sep="\t", index=False)

    if debug_fallback_out_tsv is not None:
        fallback_stages = ["lane_trim_sample_name", "lane_trim_report_stem", "stripped_unique"]
        fallback_debug = agg_merged[agg_merged["map_stage"].isin(fallback_stages)][
            ["sample", "map_stage", "run_id", "sample_name"]
        ].drop_duplicates().copy()
        fallback_debug = fallback_debug.sort_values(["map_stage", "sample"], ascending=[True, True])
        debug_fallback_out_tsv.parent.mkdir(parents=True, exist_ok=True)
        fallback_debug.to_csv(debug_fallback_out_tsv, sep="\t", index=False)

    agg = agg_merged

    if agg["run_id"].isna().any():
        missing = agg[agg["run_id"].isna()][["sample"]].drop_duplicates().head(10)
        raise ValueError(
            "Some alignment samples could not be mapped to run_id via kk_report_name_16GB stem. "
            f"Examples:\n{missing.to_string(index=False)}"
        )

    # alignment tier (V1–V3 aware)
    agg["align_confirmed"] = (agg["align_mapped_reads"] >= confirmed_reads) & (agg["align_max_breadth"] >= confirmed_breadth)

    # NCmax per run/genus
    meta_types = meta[["run_id", "sample_name", "sample_type"]].drop_duplicates()
    agg = agg.merge(meta_types, on=["run_id", "sample_name"], how="left")

    ncmax = defaultdict(lambda: defaultdict(int))
    for _, r in agg[agg["sample_type"].isin(["NC", "NTC"])].iterrows():
        run_id = int(r["run_id"])
        g = r["genus"]
        ncmax[run_id][g] = max(ncmax[run_id][g], int(r["align_mapped_reads"]))

    # also expose the NC baseline used (for debug)
    agg_keys = list(zip(agg["run_id"].astype(int), agg["genus"].astype(str)))
    ncmax_map = {(int(run_id), str(genus)): int(v)
                 for run_id, genus_dict in ncmax.items()
                 for genus, v in genus_dict.items()}
    agg["align_ncmax"] = pd.Series(agg_keys, index=agg.index).map(ncmax_map).fillna(0).astype(int)

    agg["align_final_pos"] = (
        agg["align_confirmed"] &
        (
            (agg["align_ncmax"] <= 0) |
            (agg["align_mapped_reads"] >= (fold_over_nc * agg["align_ncmax"]))
        )
    )

    return agg[[
        "run_id", "sample_name", "genus",
        "align_mapped_reads", "align_max_breadth",
        "align_confirmed", "align_ncmax", "align_final_pos"
    ]]


# -----------------------------
# Module 4: Final interpretation (validation mode)
# -----------------------------

def module_final_interpretation(meta: pd.DataFrame,
                                kk16: pd.DataFrame,
                                kkr: pd.DataFrame,
                                aln: pd.DataFrame,
                                out_prefix: Path) -> None:
    meta_exp = meta.copy()
    meta_exp["expected_genus"] = meta_exp["expected_results"].map(clean_expected_genus)

    kk16_key = kk16[["run_id", "sample_name", "genus", "kk16_pass", "kk16_reads", "kk16_ncmax"]].copy()
    kkr_cols = ["run_id", "sample_name", "genus", "kkr_pass", "kkr_reads", "kkr_ncmax"]
    if kkr.empty:
        kkr_key = pd.DataFrame(columns=kkr_cols)
    else:
        kkr_key = kkr[kkr_cols].copy()
    aln_key = aln[["run_id", "sample_name", "genus", "align_final_pos", "align_mapped_reads", "align_ncmax", "align_max_breadth"]].copy()

    # ---------- expected detection ----------
    # Alignment-first rescue for rickettsial expectations:
    # - Primary: any confirmed alignment evidence for Orientia/Rickettsia (order-level rescue by alignment)
    # - Genus-strong: expected genus passes alignment final rule
    # NOTE: rickettsiales-only Kraken is retained for flagging/debug, not core rescue.
    aln_confirm_map = {
        (int(r["run_id"]), str(r["sample_name"]), str(r["genus"])): bool(r["align_final_pos"])
        for _, r in aln_key.iterrows()
    }
    aln_any_confirmed_map = {
        (int(run_id), str(sample_name)): bool(group["align_confirmed"].any())
        for (run_id, sample_name), group in aln.groupby(["run_id", "sample_name"])
    }
    kk16_pass_map = {
        (int(r["run_id"]), str(r["sample_name"]), str(r["genus"])): bool(r["kk16_pass"])
        for _, r in kk16_key.iterrows()
    }
    kkr_order_flag_map = {
        (int(r["run_id"]), str(r["sample_name"])): bool(r["kkr_pass"])
        for _, r in kkr_key[kkr_key["genus"] == "Rickettsiales"].iterrows()
    }

    # Keep a single column name for downstream (validation scoring uses this)
    meta_exp["expected_genus"] = meta_exp["expected_genus"].fillna("")
    sample_keys = list(zip(meta_exp["run_id"].astype(int), meta_exp["sample_name"].astype(str)))
    genus_keys = list(zip(meta_exp["run_id"].astype(int), meta_exp["sample_name"].astype(str), meta_exp["expected_genus"].astype(str)))
    any_align = pd.Series(sample_keys, index=meta_exp.index).map(aln_any_confirmed_map).fillna(False)
    genus_align = pd.Series(genus_keys, index=meta_exp.index).map(aln_confirm_map).fillna(False)
    kk16_detect = pd.Series(genus_keys, index=meta_exp.index).map(kk16_pass_map).fillna(False)

    is_rick_expected = meta_exp["expected_genus"].isin(["Orientia", "Rickettsia"])
    has_expected = meta_exp["expected_genus"].ne("")
    meta_exp["expected_detected"] = False
    meta_exp.loc[has_expected & is_rick_expected, "expected_detected"] = (
        genus_align[has_expected & is_rick_expected] |
        any_align[has_expected & is_rick_expected]
    )
    meta_exp.loc[has_expected & ~is_rick_expected, "expected_detected"] = kk16_detect[has_expected & ~is_rick_expected]

    # ---------- observed set per sample ----------
    # Observed non-rickettsiales = kk16_pass True (Tier1-equivalent by your defaults)
    obs_non_rick = kk16_key[kk16_key["kk16_pass"] == True][["run_id", "sample_name", "genus"]].copy()
    obs_non_rick["source"] = "kk16"

    # Observed rickettsiales genus = alignment final positive only
    obs_rick = aln_key[aln_key["align_final_pos"] == True][["run_id", "sample_name", "genus"]].copy()
    obs_rick["source"] = "align"

    observed_long = pd.concat([obs_non_rick, obs_rick], ignore_index=True).drop_duplicates()

    # remove duplicates if both sources call same genus (keep one)
    observed_sets = (observed_long.groupby(["run_id", "sample_name"])["genus"]
                     .apply(lambda x: sorted(set(x)))
                     .reset_index(name="observed_genera_list"))

    # ---------- expected set per sample ----------
    expected_sets = (meta_exp.groupby(["run_id", "sample_name", "sample_type"])["expected_genus"]
                     .apply(lambda x: sorted([g for g in set(x) if g]))
                     .reset_index(name="expected_genera_list"))

    # ---------- per-sample counts ----------
    merged = expected_sets.merge(observed_sets, on=["run_id", "sample_name"], how="left")
    merged["observed_genera_list"] = [x if isinstance(x, list) else [] for x in merged["observed_genera_list"]]

        # ---------- per-sample row-level counts (validation scoring) ----------
    # TP_rows counts how many expected rows are detected under the current rules (including rickettsial order-level rescue).
    row_level = meta_exp[["run_id", "sample_name", "sample_type", "expected_genus", "expected_detected"]].copy()
    row_level = row_level[row_level["expected_genus"].notna() & (row_level["expected_genus"] != "")]
    per_sample = (row_level.groupby(["run_id", "sample_name"], as_index=False)
                  .agg(expected_rows=("expected_genus", "count"),
                       TP_rows=("expected_detected", "sum")))
    per_sample["TP_rows"] = per_sample["TP_rows"].astype(int)
    per_sample["expected_rows"] = per_sample["expected_rows"].astype(int)
    per_sample["FN_rows"] = per_sample["expected_rows"] - per_sample["TP_rows"]

    merged = merged.merge(per_sample, on=["run_id", "sample_name"], how="left")
    merged["expected_rows"] = merged["expected_rows"].fillna(0).astype(int)
    merged["TP_rows"] = merged["TP_rows"].fillna(0).astype(int)
    merged["FN_rows"] = merged["FN_rows"].fillna(0).astype(int)

    # Keep overlap-based counts too (helpful for debugging and for future unknown-mode), but do not use them for pass/fail.
    tp_overlap = []
    fn_overlap = []
    fp_overlap = []
    for exp_list, obs_list in zip(merged["expected_genera_list"], merged["observed_genera_list"]):
        exp = set(exp_list)
        obs = set(obs_list)
        tp_overlap.append(len(exp & obs))
        fn_overlap.append(len(exp - obs))
        fp_overlap.append(len(obs - exp))
    merged["TP_overlap"] = tp_overlap
    merged["FN_overlap"] = fn_overlap
    merged["FP_overlap"] = fp_overlap

    # Stringify for TSV readability
    merged["expected_genera"] = [",".join(x) for x in merged["expected_genera_list"]]
    merged["observed_genera"] = [",".join(x) for x in merged["observed_genera_list"]]

    # ---------- pass/fail rules ----------
    obs_len = merged["observed_genera_list"].str.len()
    merged["pass"] = False
    merged.loc[merged["sample_type"] == "PC_MIX8", "pass"] = merged["TP_rows"] >= 6
    merged.loc[merged["sample_type"] == "MIXED4", "pass"] = merged["TP_rows"] >= 3
    merged.loc[merged["sample_type"] == "PC_SINGLE", "pass"] = merged["TP_rows"] >= 1
    merged.loc[merged["sample_type"] == "CLINICAL", "pass"] = merged["TP_rows"] >= 1  # validation-only
    merged.loc[merged["sample_type"].isin(["NC", "NTC"]), "pass"] = obs_len == 0

    # ---------- flags ----------
    # rickettsiales detected in PCs = flag, not fail
    merged["flag_rick_in_pc"] = False
    is_pc = merged["sample_type"].isin(["PC_MIX8", "PC_SINGLE", "MIXED4"])
    merged.loc[is_pc, "flag_rick_in_pc"] = [
        bool(set(gs) & {"Orientia", "Rickettsia"})
        for gs in merged.loc[is_pc, "observed_genera_list"]
    ]

    # ---------- dataset-level metrics ----------
    # Metrics computed over expected rows (row-level), not over all possible genera.
    # This yields sensitivity and (panel-defined) precision.
    total_expected = int(meta_exp["expected_genus"].astype(bool).sum())
    total_tp = int(meta_exp["expected_detected"].sum())
    total_fn = total_expected - total_tp

    # For "FP" at row-level, we define a universe = union of expected genera across all samples.
    universe = sorted(set([g for g in meta_exp["expected_genus"].unique().tolist() if g]))
    U = set(universe)

    # Count FP by treating any observed genus not expected for that sample but within the universe as FP.
    fp_panel = 0
    tn_panel = 0
    for _, r in merged.iterrows():
        exp = set(r["expected_genera_list"])
        obs = set(r["observed_genera_list"])
        panel_neg = U - exp
        fp_panel += len((obs & U) - exp)
        tn_panel += len(panel_neg - (obs & U))

    sensitivity = total_tp / total_expected if total_expected else 0.0
    specificity_panel = tn_panel / (tn_panel + fp_panel) if (tn_panel + fp_panel) else 0.0
    precision_panel = total_tp / (total_tp + fp_panel) if (total_tp + fp_panel) else 0.0
    accuracy_panel = (total_tp + tn_panel) / (total_expected + tn_panel + fp_panel) if (total_expected + tn_panel + fp_panel) else 0.0

    metrics = pd.DataFrame([{
        "expected_rows_total": total_expected,
        "TP_rows": total_tp,
        "FN_rows": total_fn,
        "FP_panel": fp_panel,
        "TN_panel": tn_panel,
        "sensitivity": sensitivity,
        "specificity_panel": specificity_panel,
        "precision_panel": precision_panel,
        "accuracy_panel": accuracy_panel,
        "panel_universe_size": len(U),
    }])

    # ---------- control precision (reproducibility-style) ----------
    # What you actually want for precision in this validation context is:
    #   - binary reproducibility: % of controls that PASS within each control category
    #   - analyte detection rate: mean(TP_rows/expected_rows) for positive controls
    # NTC/NC are treated as the same negative class; each NTC instance counts separately.
    control = merged.copy()
    control["is_control"] = control["sample_type"].isin(["NTC", "NC", "PC_MIX8", "PC_SINGLE", "MIXED4"])
    control = control[control["is_control"]].copy()

    control["det_rate"] = pd.NA
    has_exp_rows = control["expected_rows"] > 0
    control.loc[has_exp_rows, "det_rate"] = (
        control.loc[has_exp_rows, "TP_rows"] / control.loc[has_exp_rows, "expected_rows"]
    ).astype(float)
    control["control_group"] = control["sample_type"].replace({"NC": "NEG_CONTROL", "NTC": "NEG_CONTROL"})

    ctrl_summary = (control.groupby("control_group", as_index=False)
                    .agg(n_samples=("sample_name", "nunique"),
                         n_instances=("sample_name", "size"),
                         pass_rate=("pass", "mean"),
                         det_rate_mean=("det_rate", "mean")))
    ctrl_summary = ctrl_summary.rename(columns={"control_group": "sample_type"})
    ctrl_summary["pass_rate"] = ctrl_summary["pass_rate"].fillna(0.0)
    ctrl_summary["det_rate_mean"] = ctrl_summary["det_rate_mean"].astype("Float64")

    # Add rescue transparency flags to per-sample output
    merged_sample_keys = list(zip(merged["run_id"].astype(int), merged["sample_name"].astype(str)))
    merged["flag_rick_align_rescue"] = pd.Series(merged_sample_keys, index=merged.index).map(aln_any_confirmed_map).fillna(False)
    merged["flag_kkr_order_screen"] = pd.Series(merged_sample_keys, index=merged.index).map(kkr_order_flag_map).fillna(False)

    # ---------- reviewer-facing primary metrics ----------
    def pass_rate_for(group_name: str):
        subset = merged[merged["sample_type"] == group_name]
        if subset.empty:
            return pd.NA
        return float(subset["pass"].mean())

    def det_rate_for(group_name: str):
        subset = control[control["sample_type"] == group_name]
        if subset.empty:
            return pd.NA
        if not subset["det_rate"].notna().any():
            return pd.NA
        return float(subset["det_rate"].dropna().mean())

    primary_metrics = pd.DataFrame([{
        "expected_rows_total": total_expected,
        "TP_rows": total_tp,
        "FN_rows": total_fn,
        "sensitivity_row_level": sensitivity,
        "clinical_n": int((merged["sample_type"] == "CLINICAL").sum()),
        "clinical_pass_rate": pass_rate_for("CLINICAL"),
        "controls_n_instances": int(control.shape[0]),
        "controls_pass_rate": float(control["pass"].mean()) if len(control) else pd.NA,
        "neg_control_n": int((merged["sample_type"].isin(["NC", "NTC"])).sum()),
        "neg_control_pass_rate": pass_rate_for("NTC"),
        "pc_single_n": int((merged["sample_type"] == "PC_SINGLE").sum()),
        "pc_single_pass_rate": pass_rate_for("PC_SINGLE"),
        "pc_single_det_rate_mean": det_rate_for("PC_SINGLE"),
        "pc_mix8_n": int((merged["sample_type"] == "PC_MIX8").sum()),
        "pc_mix8_pass_rate": pass_rate_for("PC_MIX8"),
        "pc_mix8_det_rate_mean": det_rate_for("PC_MIX8"),
        "mixed4_n": int((merged["sample_type"] == "MIXED4").sum()),
        "mixed4_pass_rate": pass_rate_for("MIXED4"),
        "mixed4_det_rate_mean": det_rate_for("MIXED4"),
    }])


# ---------- write outputs ----------
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    meta_exp.to_csv(f"{out_prefix}.expected_row_level.tsv", sep="\t", index=False)
    merged.drop(columns=["expected_genera_list", "observed_genera_list"]).to_csv(
        f"{out_prefix}.validation_summary.tsv", sep="\t", index=False
    )
    metrics.to_csv(f"{out_prefix}.dataset_metrics.tsv", sep="	", index=False)
    primary_metrics.to_csv(f"{out_prefix}.primary_metrics.tsv", sep="\t", index=False)
    ctrl_summary.to_csv(f"{out_prefix}.control_precision.tsv", sep="	", index=False)

    kk16.to_csv(f"{out_prefix}.kk16g_calls.tsv", sep="\t", index=False)
    kkr.to_csv(f"{out_prefix}.kkrick_calls.tsv", sep="\t", index=False)
    aln.to_csv(f"{out_prefix}.align_calls.tsv", sep="\t", index=False)


def module_final_interpretation_routine(meta: pd.DataFrame,
                                        kk16: pd.DataFrame,
                                        kkr: pd.DataFrame,
                                        aln: pd.DataFrame,
                                        out_prefix: Path) -> None:
    """
    Routine (unknown-sample) mode:
    - Determine observed taxa without using expected_results fields.
    - Keep rickettsial rescue/support flags for interpretation.
    - Do NOT compute validation metrics (TP/FN/sensitivity/etc.).
    """
    kk16_key = kk16[["run_id", "sample_name", "genus", "kk16_pass", "kk16_reads", "kk16_ncmax"]].copy()
    kkr_cols = ["run_id", "sample_name", "genus", "kkr_pass", "kkr_reads", "kkr_ncmax"]
    if kkr.empty:
        kkr_key = pd.DataFrame(columns=kkr_cols)
    else:
        kkr_key = kkr[kkr_cols].copy()
    aln_key = aln[["run_id", "sample_name", "genus", "align_confirmed", "align_final_pos", "align_mapped_reads", "align_ncmax", "align_max_breadth"]].copy()

    # Observed non-rickettsiales from broad-db screening
    obs_non_rick = kk16_key[kk16_key["kk16_pass"] == True][["run_id", "sample_name", "genus"]].copy()
    obs_non_rick["source"] = "kk16"

    # Observed rickettsiales genus from alignment final positive
    obs_rick = aln_key[aln_key["align_final_pos"] == True][["run_id", "sample_name", "genus"]].copy()
    obs_rick["source"] = "align"

    observed_long = pd.concat([obs_non_rick, obs_rick], ignore_index=True).drop_duplicates()
    observed_sets = (observed_long.groupby(["run_id", "sample_name"])["genus"]
                     .apply(lambda x: sorted(set(x)))
                     .reset_index(name="observed_genera_list"))

    out = meta[["run_id", "sample_name", "sample_type"]].drop_duplicates().copy()
    out = out.merge(observed_sets, on=["run_id", "sample_name"], how="left")
    out["observed_genera_list"] = [x if isinstance(x, list) else [] for x in out["observed_genera_list"]]
    out["observed_genera"] = [",".join(x) for x in out["observed_genera_list"]]

    # Routine call status (no expected comparator)
    out["detected_any"] = out["observed_genera_list"].str.len() > 0
    out["routine_call"] = out["detected_any"].map({True: "Detected", False: "Not detected"})

    # Rickettsial transparency/support flags
    aln_any_confirmed_map = {
        (int(run_id), str(sample_name)): bool(group["align_confirmed"].any())
        for (run_id, sample_name), group in aln_key.groupby(["run_id", "sample_name"])
    }
    kkr_order_flag_map = {
        (int(r["run_id"]), str(r["sample_name"])): bool(r["kkr_pass"])
        for _, r in kkr_key[kkr_key["genus"] == "Rickettsiales"].iterrows()
    }
    out_keys = list(zip(out["run_id"].astype(int), out["sample_name"].astype(str)))
    out["flag_rick_align_rescue"] = pd.Series(out_keys, index=out.index).map(aln_any_confirmed_map).fillna(False)
    out["flag_kkr_order_screen"] = pd.Series(out_keys, index=out.index).map(kkr_order_flag_map).fillna(False)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    out.drop(columns=["observed_genera_list"]).to_csv(
        f"{out_prefix}.routine_summary.tsv", sep="\t", index=False
    )

    # Keep call-level module outputs for traceability in both modes
    kk16.to_csv(f"{out_prefix}.kk16g_calls.tsv", sep="\t", index=False)
    kkr.to_csv(f"{out_prefix}.kkrick_calls.tsv", sep="\t", index=False)
    aln.to_csv(f"{out_prefix}.align_calls.tsv", sep="\t", index=False)

# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mapping_tsv", required=True, help="AFI_optimizeProtocol.tsv")
    ap.add_argument("--kk16_dir", required=True, help="Dir containing *_kraken2_report.txt")
    ap.add_argument("--kkrick_dir", required=False,
                    help="Dir containing *.c0.report.txt (required for validation/routine, optional for single-kraken modes)")
    ap.add_argument("--minimap_summary", required=True, help="metrics/summary.tsv from align step")
    ap.add_argument("--out_prefix", required=True, help="Output prefix, e.g. results/afi_validate")
    ap.add_argument("--mode",
                    choices=["validation", "routine", "validation_single_kraken", "routine_single_kraken"],
                    default="validation",
                    help=(
                        "validation: score against expected results using separate rickettsiales Kraken reports; "
                        "routine: unknown-sample observed calls using separate rickettsiales Kraken reports; "
                        "validation_single_kraken: validation using one custom Kraken report per sample (no kk_report_name_rick required); "
                        "routine_single_kraken: routine mode using one custom Kraken report per sample (no kk_report_name_rick required)"
                    ))

    # parameters (keep defaults aligned to what worked)
    ap.add_argument("--kk16_floor", type=int, default=500)
    ap.add_argument("--kk16_fold", type=float, default=5.0)
    ap.add_argument("--kkr_floor", type=int, default=200)
    ap.add_argument("--kkr_fold", type=float, default=3.0)
    ap.add_argument("--align_confirm_reads", type=int, default=100)
    ap.add_argument("--align_confirm_breadth", type=float, default=0.25)
    ap.add_argument("--align_fold_over_nc", type=float, default=5.0)

    args = ap.parse_args()

    mapping_tsv = Path(args.mapping_tsv)
    kk16_dir = Path(args.kk16_dir)
    kkr_dir = Path(args.kkrick_dir) if args.kkrick_dir else None
    minimap = Path(args.minimap_summary)
    out_prefix = Path(args.out_prefix)
    single_kraken_mode = args.mode in {"validation_single_kraken", "routine_single_kraken"}
    validation_mode = args.mode in {"validation", "validation_single_kraken"}

    meta = pd.read_csv(mapping_tsv, sep="\t")
    # Exclude diluted/technical samples and undetermined
    EXCLUDE_SAMPLE_NAMES = {
        "003691_S11_L001",
        "019371_S12_L001",
        "027571_S13_L001",
    }
    meta["sample_name"] = meta["sample_name"].astype(str).str.strip()
    meta = meta[~meta["sample_name"].isin(EXCLUDE_SAMPLE_NAMES)].copy()
    meta = meta[~meta["sample_name"].str.startswith("Undetermined", na=False)].copy()    
    req = {"run_id", "sample_name", "kk_report_name_16GB"}
    if validation_mode:
        req.add("expected_results")
    if not single_kraken_mode:
        req.add("kk_report_name_rick")

    missing = req - set(meta.columns)
    if missing:
        raise SystemExit(f"Missing required columns in mapping TSV: {missing}")

    if not single_kraken_mode and kkr_dir is None:
        raise SystemExit("--kkrick_dir is required for mode 'validation' and 'routine'.")

    if "sample_type" in meta.columns and meta["sample_type"].notna().any():
        meta["sample_type"] = meta["sample_type"].map(normalize_sample_type)
    else:
        meta["sample_type"] = meta["sample_name"].map(infer_sample_type)

    kk16 = module_kk16g(meta, kk16_dir, floor=args.kk16_floor, fold=args.kk16_fold)
    if single_kraken_mode:
        kkr = pd.DataFrame(columns=["run_id", "sample_name", "genus", "kkr_reads", "kkr_ncmax", "kkr_pass"])
    else:
        kkr = module_kk_rick(meta, kkr_dir, floor=args.kkr_floor, fold=args.kkr_fold)

    aln = module_align_rick16S(
        meta, minimap,
        confirmed_reads=args.align_confirm_reads,
        confirmed_breadth=args.align_confirm_breadth,
        fold_over_nc=args.align_fold_over_nc,
        debug_out_tsv=Path(f"{out_prefix}.align_mapping_debug.tsv"),
        debug_fallback_out_tsv=Path(f"{out_prefix}.align_mapping_fallbacks.tsv")
    )

    if validation_mode:
        module_final_interpretation(meta, kk16, kkr, aln, out_prefix)

        print("Wrote:")
        print(f"  {out_prefix}.validation_summary.tsv")
        print(f"  {out_prefix}.expected_row_level.tsv")
        print(f"  {out_prefix}.dataset_metrics.tsv")
        print(f"  {out_prefix}.primary_metrics.tsv")
        print(f"  {out_prefix}.control_precision.tsv")
        print(f"  {out_prefix}.kk16g_calls.tsv")
        print(f"  {out_prefix}.kkrick_calls.tsv")
        print(f"  {out_prefix}.align_calls.tsv")
        print(f"  {out_prefix}.align_mapping_debug.tsv")
        print(f"  {out_prefix}.align_mapping_fallbacks.tsv")
    else:
        module_final_interpretation_routine(meta, kk16, kkr, aln, out_prefix)

        print("Wrote:")
        print(f"  {out_prefix}.routine_summary.tsv")
        print(f"  {out_prefix}.kk16g_calls.tsv")
        print(f"  {out_prefix}.kkrick_calls.tsv")
        print(f"  {out_prefix}.align_calls.tsv")
        print(f"  {out_prefix}.align_mapping_debug.tsv")
        print(f"  {out_prefix}.align_mapping_fallbacks.tsv")


if __name__ == "__main__":
    main()
