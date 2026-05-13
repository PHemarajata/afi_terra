#!/usr/bin/env bash
# Build AFI_SheetBuilder for Linux (AMD64 or ARM64 — build on target arch).
# Run from the build/ directory.
# Produces: dist/AFI_SheetBuilder
set -euo pipefail

pip install -r ../requirements.txt pyinstaller

pyinstaller \
  --onefile \
  --clean \
  --name AFI_SheetBuilder \
  --add-data "../aphl_style.qss:." \
  --add-data "../assets:assets" \
  ../terra_sheet_builder.py

echo "Built: dist/AFI_SheetBuilder"
