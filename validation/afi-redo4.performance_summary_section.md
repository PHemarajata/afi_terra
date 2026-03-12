# Performance Summary (AFI 16S Validation)

This section mirrors the reviewer-style summary layout (category table, 2x2 table, and statistics with 95% CI).

## 1) Summary of assay performance

| Sample Category | Number Tested | Number Concordant | Notes |
|---|---:|---:|---|
| Positive Samples | 43 | 31 | Includes CLINICAL + positive controls with expected targets (n=43) |
| Negative Samples | 5 | 5 | Negative controls only (NTC/NC), expected none (n=5) |
| Accuracy |  | 36/48 = 75.0% | 95% CI 61.22% to 85.08% |

## 2) 2×2 concordance matrix (Expected vs Sequencing call)

| Expected Result | Sequencing Detected | Sequencing Not Detected |
|---|---:|---:|
| Detected | 31 | 12 |
| Not Detected | 0 | 5 |

Sensitivity = 72.09%    |    Specificity = 100.00%

## 3) Statistics

| Statistic | Value | 95% CI |
|---|---:|---|
| Analytical Sensitivity | 72.09% | 57.31% to 83.25% |
| Analytical Specificity | 100.00% | 56.55% to 100.00% |
| Overall Accuracy | 75.00% | 61.22% to 85.08% |

### Notes on definitions
- Positive sample = sample with at least one expected target row (`expected_rows > 0`).
- Negative sample = sample with no expected target row (`expected_rows = 0`; NTC/NC).
- Sequencing detected for positives = at least one expected target detected (`TP_rows > 0`).
- Sequencing detected for negatives = any observed genus called (non-empty observed list).