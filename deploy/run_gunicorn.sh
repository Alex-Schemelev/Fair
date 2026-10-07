#!/usr/bin/env bash
# Production start (Linux). Run from anywhere; the script cds to the project root.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export HOST="${HOST:-0.0.0.0}"
export PORT="${PORT:-8000}"
export OPEN_BROWSER=0
# BLAST human stays off unless explicitly enabled in the UI / env for a run.
ulimit -n "${NOFILE_LIMIT:-65535}"
GUNICORN="${GUNICORN:-$ROOT/.venv/bin/gunicorn}"
if [[ ! -x "$GUNICORN" ]]; then
  GUNICORN="$(command -v gunicorn)"
fi
exec "$GUNICORN" --config "$ROOT/deploy/gunicorn.conf.py" web_app:app
