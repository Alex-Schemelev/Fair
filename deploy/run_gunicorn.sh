#!/usr/bin/env bash
# Production start (Linux). Run from project root with venv activated.
set -euo pipefail
cd "$(dirname "$0")/.."
export HOST="${HOST:-0.0.0.0}"
export PORT="${PORT:-8000}"
export OPEN_BROWSER=0
# BLAST human stays off unless explicitly enabled in the UI / env for a run.
exec gunicorn \
  --bind "${HOST}:${PORT}" \
  --workers 1 \
  --threads 4 \
  --timeout 1800 \
  --access-logfile - \
  --error-logfile - \
  "web_app:app"
