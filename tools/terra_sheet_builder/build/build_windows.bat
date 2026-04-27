@echo off
REM Run on Windows (x64).
REM Produces: dist\AFI_SheetBuilder.exe
pip install -r ..\requirements.txt pyinstaller
pyinstaller --onefile --windowed --name AFI_SheetBuilder ..\terra_sheet_builder.py
echo Built: dist\AFI_SheetBuilder.exe
