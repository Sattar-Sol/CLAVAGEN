@echo off
REM CLAVAGEN - double-click to start (Windows)
title CLAVAGEN
cd /d "%~dp0"

if exist ".venv\Scripts\pythonw.exe" goto launch

set "PY="
where py >nul 2>nul
if not errorlevel 1 set "PY=py -3"
if defined PY goto checkver
where python >nul 2>nul
if not errorlevel 1 set "PY=python"
if not defined PY goto nopython
:checkver
%PY% -c "import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)" >nul 2>nul
if errorlevel 1 goto nopython

echo First start: setting up CLAVAGEN (one time only, about 1 minute)...
%PY% -m venv .venv
if errorlevel 1 goto noenv
".venv\Scripts\python.exe" -m pip install --quiet --disable-pip-version-check certifi
".venv\Scripts\python.exe" -m pip install --quiet --disable-pip-version-check pywebview
if errorlevel 1 echo (Native window component could not be installed - an Edge app window will be used instead.)
".venv\Scripts\python.exe" -m pip install --quiet --disable-pip-version-check biopython >nul 2>nul
echo checked> ".venv\window_checked.txt"

:launch
REM window component missing (e.g. no internet at first start): try again, at most once a day
".venv\Scripts\python.exe" -c "import webview" >nul 2>nul
if not errorlevel 1 goto certs
forfiles /p ".venv" /m window_checked.txt /d 0 >nul 2>nul
if not errorlevel 1 goto certs
echo Installing the CLAVAGEN window component...
".venv\Scripts\python.exe" -m pip install --quiet --disable-pip-version-check --timeout 20 --retries 1 pywebview >nul 2>nul
echo checked> ".venv\window_checked.txt"
:certs
REM make sure secure-connection certificates are available (needed for ClinVar/NCBI/UniProt/gnomAD/Ensembl)
".venv\Scripts\python.exe" -c "import certifi" >nul 2>nul || ".venv\Scripts\python.exe" -m pip install --quiet --disable-pip-version-check certifi >nul 2>nul
start "" ".venv\Scripts\pythonw.exe" -m clavagen app
exit /b 0

:noenv
echo Could not create the local environment; starting with system Python.
start "" %PY% -m clavagen app
exit /b 0

:nopython
echo Python 3.8 or newer was not found.
echo Install it from https://www.python.org/downloads/
echo IMPORTANT: tick "Add Python to PATH" during installation, then double-click this file again.
start "" "https://www.python.org/downloads/"
pause
exit /b 1
