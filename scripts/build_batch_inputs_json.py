#!/usr/bin/env python3
import argparse
import csv
import json
import re
from pathlib import Path

VALID_MODES = {"validation", "routine"}
VALID_CLASSIFIERS = {"single", "double"}


def parse_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    text = value.strip().lower()
    if text == "":
        return None
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    raise ValueError(f"Invalid boolean value: {value}")


def normalize_row(row: dict[str, str], row_index: int) -> dict:
    required = [
        "sample_id",
        "sample_type",
        "mode",
        "classifier_mode",
        "r1_fastq",
        "r2_fastq",
    ]

    missing = [key for key in required if key not in row or row[key].strip() == ""]
    if missing:
        raise ValueError(f"Row {row_index}: missing required columns/values: {', '.join(missing)}")

    mode = row["mode"].strip().lower()
    if mode not in VALID_MODES:
        raise ValueError(f"Row {row_index}: invalid mode '{row['mode']}', expected one of {sorted(VALID_MODES)}")

    classifier_mode = row["classifier_mode"].strip().lower()
    if classifier_mode not in VALID_CLASSIFIERS:
        raise ValueError(
            f"Row {row_index}: invalid classifier_mode '{row['classifier_mode']}', expected one of {sorted(VALID_CLASSIFIERS)}"
        )

    sample = {
        "sample_id": row["sample_id"].strip(),
        "sample_type": row["sample_type"].strip(),
        "mode": mode,
        "classifier_mode": classifier_mode,
        "r1_fastq": row["r1_fastq"].strip(),
        "r2_fastq": row["r2_fastq"].strip(),
    }

    expected_taxon = row.get("expected_taxon", "").strip()
    if expected_taxon:
        sample["expected_taxon"] = expected_taxon

    use_human_scrub = parse_bool(row.get("use_human_scrub"))
    if use_human_scrub is not None:
        sample["use_human_scrub"] = use_human_scrub

    return sample


def build_inputs(args: argparse.Namespace, samples: list[dict]) -> dict:
    inputs = {
        "AFI_Rickettsiales_Batch.samples": samples,
    }

    if args.rickettsiales_panel:
        inputs["AFI_Rickettsiales_Batch.rickettsiales_panel"] = args.rickettsiales_panel
    if args.ntc_background:
        inputs["AFI_Rickettsiales_Batch.ntc_background"] = args.ntc_background

    if args.default_use_human_scrub is not True:
        inputs["AFI_Rickettsiales_Batch.default_use_human_scrub"] = args.default_use_human_scrub
    if args.classify_threads != 16:
        inputs["AFI_Rickettsiales_Batch.classify_threads"] = args.classify_threads

    if args.kraken_db_16g:
        inputs["AFI_Rickettsiales_Batch.kraken_db_16g"] = args.kraken_db_16g
    if args.kraken_db_rick:
        inputs["AFI_Rickettsiales_Batch.kraken_db_rick"] = args.kraken_db_rick
    if args.centrifuger_db:
        inputs["AFI_Rickettsiales_Batch.centrifuger_db"] = args.centrifuger_db

    return inputs


def normalize_sample_type_from_mapping(value: str) -> str:
    text = (value or "").strip().upper()
    if text in {"PC", "POSITIVE_CONTROL", "PC_MIX8", "PC_SINGLE", "MIXED4"}:
        return "PC"
    if text in {"NTC", "NC", "NEGATIVE_CONTROL"}:
        return "NTC"
    return "clinical"


def normalize_expected_taxon(value: str) -> str:
    text = (value or "").strip()
    if text.upper() in {"", "NA", "N/A", "NONE", "NULL"}:
        return ""
    return text


def derive_fastq_sample_id(report_name: str, fallback_sample: str) -> str:
    text = (report_name or "").strip()
    if text:
        stem = re.sub(r"(_kraken2_report\.txt(\.gz)?|\.report\.txt(\.gz)?|_kreport\.txt(\.gz)?)$", "", text)
        if stem:
            return stem
    return fallback_sample.strip()


def infer_mode(expected_taxon: str, mode_policy: str) -> str:
    if mode_policy == "all_validation":
        return "validation"
    if mode_policy == "all_routine":
        return "routine"
    return "validation" if expected_taxon else "routine"


def load_samples_from_mapping(args: argparse.Namespace) -> list[dict]:
    mapping_path = Path(args.mapping_tsv)
    if not mapping_path.exists():
        raise SystemExit(f"Mapping TSV not found: {mapping_path}")

    samples: list[dict] = []
    with mapping_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("Mapping TSV appears empty")

        required_cols = {"sample_name", "sample_type", "expected_results", "kk_report_name_16GB"}
        missing = required_cols - set(reader.fieldnames)
        if missing:
            raise SystemExit(f"Mapping TSV missing required columns: {', '.join(sorted(missing))}")

        for row_index, row in enumerate(reader, start=2):
            sample_name = (row.get("sample_name") or "").strip()
            if not sample_name or sample_name.startswith("Undetermined"):
                continue

            expected_taxon = normalize_expected_taxon(row.get("expected_results", ""))
            mode = infer_mode(expected_taxon, args.mode_policy)

            sample_id = derive_fastq_sample_id(row.get("kk_report_name_16GB", ""), sample_name)
            r1 = f"{args.fastq_uri_prefix.rstrip('/')}/{sample_id}_R1.fastq.gz"
            r2 = f"{args.fastq_uri_prefix.rstrip('/')}/{sample_id}_R2.fastq.gz"

            item = {
                "sample_id": sample_id,
                "sample_type": normalize_sample_type_from_mapping(row.get("sample_type", "")),
                "mode": mode,
                "classifier_mode": args.default_classifier_mode,
                "r1_fastq": r1,
                "r2_fastq": r2,
            }

            if expected_taxon:
                item["expected_taxon"] = expected_taxon

            if args.per_sample_use_human_scrub is not None:
                item["use_human_scrub"] = args.per_sample_use_human_scrub

            samples.append(item)

    if not samples:
        raise SystemExit("No usable sample rows found in mapping TSV")

    return samples


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Terra input JSON for AFI_Rickettsiales_Batch from a sample sheet TSV."
    )
    parser.add_argument("--sample-sheet", help="TSV with one sample per row")
    parser.add_argument("--mapping-tsv", help="AFI mapping TSV (e.g., AFI_optimizeProtocol.tsv)")
    parser.add_argument("--out-json", required=True, help="Output JSON path")

    parser.add_argument("--rickettsiales-panel", help="GS path to rickettsiales_16S_panel.fasta")
    parser.add_argument("--ntc-background", help="GS path to ntc_background.tsv")

    parser.add_argument("--kraken-db-16g", help="Kraken2 16G database path")
    parser.add_argument("--kraken-db-rick", help="Kraken2 rickettsiales database path")
    parser.add_argument("--centrifuger-db", help="Centrifuger database path")

    parser.add_argument("--default-use-human-scrub", action="store_true", default=True)
    parser.add_argument("--no-default-use-human-scrub", dest="default_use_human_scrub", action="store_false")
    parser.add_argument("--classify-threads", type=int, default=16)

    parser.add_argument("--mode-policy", choices=["auto", "all_validation", "all_routine"], default="auto")
    parser.add_argument("--default-classifier-mode", choices=sorted(VALID_CLASSIFIERS), default="double")
    parser.add_argument("--fastq-uri-prefix", help="Required with --mapping-tsv, e.g. gs://bucket/fastq")
    parser.add_argument("--per-sample-use-human-scrub", choices=["true", "false"])

    args = parser.parse_args()

    out_path = Path(args.out_json)

    if bool(args.sample_sheet) == bool(args.mapping_tsv):
        raise SystemExit("Provide exactly one of --sample-sheet or --mapping-tsv")

    args.per_sample_use_human_scrub = parse_bool(args.per_sample_use_human_scrub)

    samples: list[dict] = []
    if args.sample_sheet:
        sheet_path = Path(args.sample_sheet)
        if not sheet_path.exists():
            raise SystemExit(f"Sample sheet not found: {sheet_path}")

        with sheet_path.open("r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise SystemExit("Sample sheet appears empty")

            for row_index, row in enumerate(reader, start=2):
                samples.append(normalize_row(row, row_index))
    else:
        if not args.fastq_uri_prefix:
            raise SystemExit("--fastq-uri-prefix is required with --mapping-tsv")
        samples = load_samples_from_mapping(args)

    if not samples:
        raise SystemExit("Sample sheet contains no data rows")

    payload = build_inputs(args, samples)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {out_path} with {len(samples)} samples")


if __name__ == "__main__":
    main()
