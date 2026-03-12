#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from typing import Iterable


def normalize_text(value: str | None) -> str:
    return (value or "").strip()


def normalize_taxa_csv(value: str | None) -> str:
    items = [x.strip() for x in (value or "").split(",") if x.strip()]
    return ",".join(sorted(set(items)))


def discover_tsvs(path_text: str, suffix: str) -> list[Path]:
    path = Path(path_text)
    if not path.exists():
        raise SystemExit(f"Path does not exist: {path}")

    if path.is_file():
        if path.name.endswith(suffix):
            return [path]
        raise SystemExit(f"Expected file ending with '{suffix}', got: {path}")

    files = sorted(path.rglob(f"*{suffix}"))
    if not files:
        raise SystemExit(f"No files ending with '{suffix}' found under: {path}")
    return files


def read_rows(files: Iterable[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for file_path in files:
        with file_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                continue
            for row in reader:
                normalized = {k: normalize_text(v) for k, v in row.items()}
                normalized["__source_file"] = str(file_path)
                rows.append(normalized)
    return rows


def index_by_sample(rows: list[dict[str, str]], label: str) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for i, row in enumerate(rows, start=1):
        sample_id = normalize_text(row.get("sample_id") or row.get("sample"))
        if not sample_id:
            raise SystemExit(f"{label}: missing sample_id/sample at row {i} from {row.get('__source_file', 'unknown')}.")
        if sample_id in indexed:
            raise SystemExit(f"{label}: duplicate sample_id detected: {sample_id}")
        indexed[sample_id] = row
    return indexed


def compare_validation(single: dict[str, dict[str, str]], double: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    all_samples = sorted(set(single) | set(double))
    diffs: list[dict[str, str]] = []

    for sample_id in all_samples:
        s = single.get(sample_id)
        d = double.get(sample_id)

        if s is None:
            diffs.append({"sample_id": sample_id, "status": "only_in_double"})
            continue
        if d is None:
            diffs.append({"sample_id": sample_id, "status": "only_in_single"})
            continue

        s_taxa = normalize_taxa_csv(s.get("detected_taxa"))
        d_taxa = normalize_taxa_csv(d.get("detected_taxa"))
        s_result = normalize_text(s.get("validation_result"))
        d_result = normalize_text(d.get("validation_result"))

        if s_taxa != d_taxa or s_result != d_result:
            diffs.append(
                {
                    "sample_id": sample_id,
                    "status": "changed",
                    "single_validation_result": s_result,
                    "double_validation_result": d_result,
                    "single_detected_taxa": s_taxa,
                    "double_detected_taxa": d_taxa,
                    "single_sample_type": normalize_text(s.get("sample_type")),
                    "double_sample_type": normalize_text(d.get("sample_type")),
                    "expected_taxon": normalize_text(s.get("expected_taxon") or d.get("expected_taxon")),
                }
            )

    return diffs


def compare_routine(single: dict[str, dict[str, str]], double: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    all_samples = sorted(set(single) | set(double))
    diffs: list[dict[str, str]] = []

    for sample_id in all_samples:
        s = single.get(sample_id)
        d = double.get(sample_id)

        if s is None:
            diffs.append({"sample_id": sample_id, "status": "only_in_double"})
            continue
        if d is None:
            diffs.append({"sample_id": sample_id, "status": "only_in_single"})
            continue

        s_taxa = normalize_taxa_csv(s.get("taxa_present"))
        d_taxa = normalize_taxa_csv(d.get("taxa_present"))
        s_n = normalize_text(s.get("n_taxa_present"))
        d_n = normalize_text(d.get("n_taxa_present"))

        if s_taxa != d_taxa or s_n != d_n:
            diffs.append(
                {
                    "sample_id": sample_id,
                    "status": "changed",
                    "single_n_taxa_present": s_n,
                    "double_n_taxa_present": d_n,
                    "single_taxa_present": s_taxa,
                    "double_taxa_present": d_taxa,
                    "single_sample_type": normalize_text(s.get("sample_type")),
                    "double_sample_type": normalize_text(d.get("sample_type")),
                }
            )

    return diffs


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    if not rows:
        path.write_text("sample_id\tstatus\n", encoding="utf-8")
        return

    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, str]]) -> tuple[int, int, int]:
    changed = sum(1 for r in rows if r.get("status") == "changed")
    only_single = sum(1 for r in rows if r.get("status") == "only_in_single")
    only_double = sum(1 for r in rows if r.get("status") == "only_in_double")
    return changed, only_single, only_double


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare AFI single-classifier vs double-classifier outputs for validation and routine summaries."
    )
    parser.add_argument("--single-validation", help="Validation summary TSV file or directory for single-classifier run")
    parser.add_argument("--double-validation", help="Validation summary TSV file or directory for double-classifier run")
    parser.add_argument("--single-routine", help="Routine summary TSV file or directory for single-classifier run")
    parser.add_argument("--double-routine", help="Routine summary TSV file or directory for double-classifier run")
    parser.add_argument("--out-prefix", required=True, help="Output prefix (e.g., comparison/single_vs_double)")
    args = parser.parse_args()

    has_validation = bool(args.single_validation and args.double_validation)
    has_routine = bool(args.single_routine and args.double_routine)

    if not has_validation and not has_routine:
        raise SystemExit(
            "Provide at least one pair: (--single-validation and --double-validation) "
            "or (--single-routine and --double-routine)."
        )

    out_prefix = Path(args.out_prefix)
    report_lines: list[str] = []

    if has_validation:
        single_val_files = discover_tsvs(args.single_validation, "validation_summary.tsv")
        double_val_files = discover_tsvs(args.double_validation, "validation_summary.tsv")

        single_val_idx = index_by_sample(read_rows(single_val_files), "single-validation")
        double_val_idx = index_by_sample(read_rows(double_val_files), "double-validation")

        val_diffs = compare_validation(single_val_idx, double_val_idx)
        val_out = Path(f"{out_prefix}.validation_differences.tsv")
        write_tsv(val_out, val_diffs)

        changed, only_single, only_double = summarize(val_diffs)
        report_lines.append(f"validation_total_single\t{len(single_val_idx)}")
        report_lines.append(f"validation_total_double\t{len(double_val_idx)}")
        report_lines.append(f"validation_changed\t{changed}")
        report_lines.append(f"validation_only_in_single\t{only_single}")
        report_lines.append(f"validation_only_in_double\t{only_double}")
        report_lines.append(f"validation_diff_file\t{val_out}")

    if has_routine:
        single_rtn_files = discover_tsvs(args.single_routine, "routine_summary.tsv")
        double_rtn_files = discover_tsvs(args.double_routine, "routine_summary.tsv")

        single_rtn_idx = index_by_sample(read_rows(single_rtn_files), "single-routine")
        double_rtn_idx = index_by_sample(read_rows(double_rtn_files), "double-routine")

        rtn_diffs = compare_routine(single_rtn_idx, double_rtn_idx)
        rtn_out = Path(f"{out_prefix}.routine_differences.tsv")
        write_tsv(rtn_out, rtn_diffs)

        changed, only_single, only_double = summarize(rtn_diffs)
        report_lines.append(f"routine_total_single\t{len(single_rtn_idx)}")
        report_lines.append(f"routine_total_double\t{len(double_rtn_idx)}")
        report_lines.append(f"routine_changed\t{changed}")
        report_lines.append(f"routine_only_in_single\t{only_single}")
        report_lines.append(f"routine_only_in_double\t{only_double}")
        report_lines.append(f"routine_diff_file\t{rtn_out}")

    summary_out = Path(f"{out_prefix}.summary.tsv")
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    print(f"Wrote {summary_out}")
    if has_validation:
        print(f"Wrote {out_prefix}.validation_differences.tsv")
    if has_routine:
        print(f"Wrote {out_prefix}.routine_differences.tsv")


if __name__ == "__main__":
    main()
