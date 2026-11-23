#!/bin/bash
# CLAVAGEN - start (Linux)
cd "$(dirname "$0")" || exit 1
PY=$(command -v python3 || command -v python)
[ -z "$PY" ] && { echo "Python 3.8+ not found (install python3)."; exit 1; }
if [ ! -x .venv/bin/python ] && "$PY" -m venv .venv 2>/dev/null; then
  echo "First start: setting up CLAVAGEN (one time only)..."
  .venv/bin/python -m pip install --quiet certifi
  .venv/bin/python -m pip install --quiet pywebview || echo "(pywebview not installed - Chrome/Chromium app window will be used)"
  .venv/bin/python -m pip install --quiet biopython || echo "(Biopython not installed - optional, the built-in aligner is used)"
  touch .venv/.window_checked
fi
# window component missing: try again, at most once a day
if [ -x .venv/bin/python ] && ! .venv/bin/python -c "import webview" >/dev/null 2>&1; then
  if [ ! -f .venv/.window_checked ] || [ -n "$(find .venv/.window_checked -mtime +0 2>/dev/null)" ]; then
    .venv/bin/python -m pip install --quiet --timeout 20 --retries 1 pywebview || true
    touch .venv/.window_checked
  fi
fi
[ -x .venv/bin/python ] && PY=.venv/bin/python
exec "$PY" -m clavagen app
