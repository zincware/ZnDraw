"""Gunicorn configuration for production deployment.

This configuration is optimized for Flask-SocketIO with Redis message queue.
Uses threaded workers with simple-websocket for WebSocket support.
"""

import os

# Server socket
bind = "0.0.0.0:5000"
backlog = 2048

# Worker processes
# Start with conservative number of workers (1 per 16 cores as baseline)
# With 128 cores, we can scale up to 32+ workers if needed
workers = int(os.getenv("GUNICORN_WORKERS", "8"))

# Worker class - use sync (threaded) for Flask-SocketIO without gevent/eventlet
# Flask-SocketIO will use simple-websocket for WebSocket support
worker_class = "sync"

# Threads per worker
# Each worker can handle up to 'threads' concurrent requests
# 100 threads per worker = 800 concurrent connections with 8 workers
threads = int(os.getenv("GUNICORN_THREADS", "100"))

# Worker connections (for async workers like gevent/eventlet)
# Not used with sync workers, but including for future reference
# worker_connections = 1000

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
    server.log.info(f"Workers: {workers}, Threads per worker: {threads}")
    server.log.info(f"Total concurrent connections: {workers * threads}")
    server.log.info(f"Timeout: {timeout}s")
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
