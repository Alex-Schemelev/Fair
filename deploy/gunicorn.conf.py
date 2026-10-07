"""Gunicorn limits for the primer app.

One gthread worker with the default 1000 accepted connections filled the
process file-descriptor limit (1024) and stopped accepting new clients.
Sync workers plus the deadlines in ``deadline_worker`` close that hole.
"""

import os
import resource

host = os.environ.get("HOST", "0.0.0.0")
port = os.environ.get("PORT", "8000")
bind = f"{host}:{port}"

# Two processes so one occupied worker still leaves the port responsive.
# Threads stay at 1: with threads > 1 gunicorn switches a sync worker to
# gthread, and gthread does not abort a stuck request.
workers = 2
threads = 1
worker_class = "deploy.deadline_worker.DeadlineWorker"

# Long enough for a PCR / BLAST run. Idle clients are cut off earlier by
# DeadlineWorker, not by this value.
timeout = 1800
graceful_timeout = 30
keepalive = 2

# Do not park hundreds of unaccepted connections in the listen queue.
backlog = 64

# Recycle workers so a slow leak cannot accumulate for weeks.
max_requests = 500
max_requests_jitter = 50

limit_request_line = 4094
limit_request_fields = 40
limit_request_field_size = 8190

accesslog = "-"
errorlog = "-"
capture_output = True

NOFILE = 65535


def on_starting(server):
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    target = min(NOFILE, hard)
    if soft < target:
        resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))
        server.log.info(
            "Raised RLIMIT_NOFILE soft limit from %s to %s (hard %s)",
            soft,
            target,
            hard,
        )
