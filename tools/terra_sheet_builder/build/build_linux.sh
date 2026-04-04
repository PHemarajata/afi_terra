#!/usr/bin/env bash
# Run on Linux (AMD64 or ARM64 — build on target arch).
# Produces: dist/AFI_SheetBuilder
set -euo pipefail
pip install -r ../requirements.txt pyinstaller
pyinstaller \
  --onefile \
  --name AFI_SheetBuilder \
  ../terra_sheet_builder.py
echo "Built: dist/AFI_SheetBuilder"
