"""Gunicorn configuration for production deployment.

This configuration is optimized for Flask-SocketIO with Redis message queue.
Uses threaded workers with simple-websocket for WebSocket support.
"""

import os

# Server socket
bind = "0.0.0.0:5000"
backlog = 2048

# Worker processes
# Flask-SocketIO with threaded workers requires SINGLE worker process
# Multiple workers are only supported with eventlet/gevent async frameworks
# For horizontal scaling, deploy multiple containers instead of multiple workers
workers = int(os.getenv("GUNICORN_WORKERS", "1"))

# Worker class - use sync (threaded) for Flask-SocketIO without gevent/eventlet
# Flask-SocketIO will use simple-websocket for WebSocket support
worker_class = "sync"

# Threads per worker
# Each worker can handle up to 'threads' concurrent requests
# With 1 worker, increase threads for higher concurrency
# 200 threads = 200 concurrent connections
# For 128 cores, you can increase this to 500-1000 threads if needed
threads = int(os.getenv("GUNICORN_THREADS", "200"))

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
    server.log.info(f"Max concurrent connections: {threads}")
    server.log.info(f"Worker class: {worker_class}")
    server.log.info(f"Timeout: {timeout}s")
    server.log.info("Note: Flask-SocketIO uses 1 worker + Redis for scaling")
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
