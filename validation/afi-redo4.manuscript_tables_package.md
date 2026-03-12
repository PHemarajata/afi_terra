# AFI Validation Manuscript Tables Package (afi-redo4)

This package contains publication-ready versions of the primary validation performance table and control reproducibility table.

## Table 1. Primary validation performance summary

| Metric | Value |
|---|---:|
| Expected rows total | 81 |
| True positive rows (TP) | 69 |
| False negative rows (FN) | 12 |
| Row-level sensitivity | 85.2% |
| Clinical samples (n) | 33 |
| Clinical pass rate | 63.6% |
| Control instances (n) | 15 |
| Overall controls pass rate | 100.0% |
| Negative controls (n) | 5 |
| Negative control pass rate | 100.0% |
| PC_SINGLE controls (n) | 4 |
| PC_SINGLE pass rate | 100.0% |
| PC_SINGLE mean detection rate | 100.0% |
| PC_MIX8 controls (n) | 5 |
| PC_MIX8 pass rate | 100.0% |
| PC_MIX8 mean detection rate | 100.0% |
| MIXED4 controls (n) | 1 |
| MIXED4 pass rate | 100.0% |
| MIXED4 mean detection rate | 100.0% |

## Table 2. Control reproducibility detail

| Control Type | Unique Samples (n) | Control Instances (n) | Pass Rate | Mean Detection Rate |
|---|---:|---:|---:|---:|
| NEG_CONTROL | 5 | 5 | 100.0% | N/A (no expected analytes) |
| PC_SINGLE | 4 | 4 | 100.0% | 100.0% |
| PC_MIX8 | 5 | 5 | 100.0% | 100.0% |
| MIXED4 | 1 | 1 | 100.0% | 100.0% |

## Footnotes (recommended)

1. **Primary endpoint definition:** Row-level sensitivity = TP expected rows / total expected rows.
2. **Control precision definition:** Pass-rate reproducibility is summarized by control type; mean detection rate is reported for controls with expected analytes.
3. **Negative controls:** NC and NTC are merged as `NEG_CONTROL`; mean detection rate is not applicable.
4. **Rickettsiales interpretation:** For expected rickettsial targets, alignment-supported order-level rescue is permitted due to limited V1–V2 genus resolution; rickettsiales-only Kraken order signal is used as supportive context.

## Source files

- `validation/afi-redo4.Table1_draft_pub.tsv`
- `validation/afi-redo4.Table2_control_precision_pub.tsv`
- `validation/afi-redo4.primary_metrics.tsv`
- `validation/afi-redo4.control_precision.tsv`
