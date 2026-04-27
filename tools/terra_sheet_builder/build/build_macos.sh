#!/usr/bin/env bash
# Run on macOS (Apple Silicon or Intel).
# Produces: dist/AFI_SheetBuilder.app
set -euo pipefail
pip install -r ../requirements.txt pyinstaller
pyinstaller \
  --onefile \
  --windowed \
  --name AFI_SheetBuilder \
  ../terra_sheet_builder.py
echo "Built: dist/AFI_SheetBuilder.app"
