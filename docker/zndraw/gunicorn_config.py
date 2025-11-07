"""Gunicorn configuration for production deployment.

This configuration is optimized for Flask-SocketIO with Redis message queue.
Uses threaded workers for WebSocket support without event loop conflicts.

Key points:
- Sync worker with threads for concurrent request handling
- Single worker process required for Flask-SocketIO with threaded mode
- Threads parameter controls concurrent connections
- No event loop, compatible with Zarr v3 synchronous operations
- For horizontal scaling, deploy multiple containers instead of multiple workers
"""

import os

# Server socket
bind = "0.0.0.0:5000"
backlog = 2048

# Worker processes
# Flask-SocketIO with threaded workers requires SINGLE worker process
# For horizontal scaling, deploy multiple containers instead of multiple workers
workers = int(os.getenv("GUNICORN_WORKERS", "1"))

worker_class = "eventlet"

# # Number of threads per worker
# # Each thread handles one concurrent connection
# # Default: 100 threads for testing, increase to 500-1000 for production
# threads = int(os.getenv("GUNICORN_THREADS", "100"))

# Timeouts
# Increased timeout for large file uploads (5 minutes)
timeout = 300
graceful_timeout = 30
keepalive = 5

# Logging
accesslog = "-"  # Log to stdout
errorlog = "-"  # Log to stderr
loglevel = os.getenv("GUNICORN_LOG_LEVEL", "info")
access_log_format = '%(h)s %(l)s %(u)s %(t)s "%(r)s" %(s)s %(b)s "%(f)s" "%(a)s" %(D)s'

# Process naming
proc_name = "zndraw"

# Server mechanics
daemon = False
pidfile = None
umask = 0
user = None
group = None
tmp_upload_dir = None

# SSL (if needed in future)
# keyfile = None
# certfile = None

# Preload app for better performance and memory sharing
preload_app = True

# Max requests per worker (helps prevent memory leaks)
# Workers restart after handling this many requests
max_requests = 10000
max_requests_jitter = 1000


def on_starting(server):
    """Called just before the master process is initialized."""
    server.log.info("=" * 80)
    server.log.info("ZnDraw Gunicorn server starting")
    server.log.info(f"Workers: {workers}")
    server.log.info(f"Worker class: {worker_class}")
    # server.log.info(f"Threads per worker: {threads}")
    # server.log.info(f"Max concurrent connections: {threads}")
    server.log.info(f"Timeout: {timeout}s")
    server.log.info("Note: Threaded workers with Redis for SocketIO coordination")
    server.log.info("Note: Redis coordinates SocketIO across multiple containers")
    server.log.info("=" * 80)


def on_reload(server):
    """Called to recycle workers during a reload via SIGHUP."""
    server.log.info("Reloading workers...")


def worker_int(worker):
    """Called when a worker receives the SIGINT or SIGQUIT signal."""
    worker.log.info(f"Worker received INT or QUIT signal (PID: {worker.pid})")


def worker_abort(worker):
    """Called when a worker receives the SIGABRT signal."""
    worker.log.info(f"Worker aborted (PID: {worker.pid})")
