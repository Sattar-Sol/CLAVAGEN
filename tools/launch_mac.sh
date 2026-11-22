#!/bin/bash
# Shared macOS launcher (used by Start_Mac.command and CLAVAGEN.app)
cd "$(dirname "$0")/.." || exit 1
LOG="$HOME/Library/Logs/CLAVAGEN.log"
mkdir -p "$(dirname "$LOG")"
note() { osascript -e "display notification \"$1\" with title \"CLAVAGEN\"" >/dev/null 2>&1; echo "$1"; }

if [ ! -x .venv/bin/python ]; then
  PY=""
  for c in python3 /usr/local/bin/python3 /opt/homebrew/bin/python3 /Library/Frameworks/Python.framework/Versions/Current/bin/python3 /usr/bin/python3; do
    if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(0 if sys.version_info >= (3,8) else 1)' >/dev/null 2>&1; then
      PY="$c"; break
    fi
  done
  if [ -z "$PY" ]; then
    osascript -e 'display dialog "Python 3.8 or newer is required.\n\nInstall it from python.org, then open CLAVAGEN again." buttons {"Open python.org"} default button 1 with title "CLAVAGEN"' >/dev/null 2>&1
    open "https://www.python.org/downloads/"
    exit 1
  fi
  note "First start: setting up CLAVAGEN (one time, about 1 minute)..."
  if "$PY" -m venv .venv >>"$LOG" 2>&1; then
    .venv/bin/python -m pip install --quiet --disable-pip-version-check pywebview certifi biopython >>"$LOG" 2>&1 \
      || note "Native window component not installed - a Chrome/Edge app window or Safari will be used."
  else
    note "Could not create local environment - using system Python."
    exec "$PY" -m clavagen app >>"$LOG" 2>&1
  fi
fi
# secure-connection certificates (python.org macOS builds ship without them)
.venv/bin/python -c "import certifi" >/dev/null 2>&1 || .venv/bin/python -m pip install --quiet --disable-pip-version-check certifi >>"$LOG" 2>&1
exec .venv/bin/python -m clavagen app >>"$LOG" 2>&1
