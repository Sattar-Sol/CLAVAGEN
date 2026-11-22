#!/bin/bash
# CLAVAGEN - start (Linux)
cd "$(dirname "$0")" || exit 1
PY=$(command -v python3 || command -v python)
[ -z "$PY" ] && { echo "Python 3.8+ not found (install python3)."; exit 1; }
if [ ! -x .venv/bin/python ] && "$PY" -m venv .venv 2>/dev/null; then
  echo "First start: installing the native window component (optional)..."
  .venv/bin/python -m pip install --quiet pywebview certifi biopython || echo "(pywebview not installed - Chrome/Chromium app window will be used)"
fi
[ -x .venv/bin/python ] && PY=.venv/bin/python
exec "$PY" -m clavagen app
