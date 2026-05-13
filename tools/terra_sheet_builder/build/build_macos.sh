#!/usr/bin/env bash
# Build AFI_SheetBuilder for macOS (Apple Silicon or Intel).
# Run from the build/ directory.
# Produces: dist/AFI_SheetBuilder.app
set -euo pipefail

pip install -r ../requirements.txt pyinstaller

pyinstaller \
  --onedir \
  --windowed \
  --clean \
  --name AFI_SheetBuilder \
  --add-data "../aphl_style.qss:." \
  --add-data "../assets:assets" \
  --icon "../assets/aphl-icon.icns" \
  ../terra_sheet_builder.py

echo "Built: dist/AFI_SheetBuilder.app"
