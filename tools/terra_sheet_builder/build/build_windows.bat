@echo off
REM Build AFI_SheetBuilder for Windows (x64)
REM Run from the build\ directory.
REM Produces: dist\AFI_SheetBuilder.exe

pip install -r ..\requirements.txt pyinstaller
if errorlevel 1 exit /b 1

pyinstaller ^
  --onefile ^
  --windowed ^
  --name AFI_SheetBuilder ^
  --add-data "..\aphl_style.qss;." ^
  --add-data "..\assets;assets" ^
  --icon "..\assets\aphl-icon.ico" ^
  ..\terra_sheet_builder.py

echo.
echo Built: dist\AFI_SheetBuilder.exe
