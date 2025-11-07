"""CLI entry point with eventlet monkey patching.

This module MUST be the entry point for the zndraw CLI command.
It applies eventlet monkey patching BEFORE any other imports to ensure
proper async behavior for Flask-SocketIO and Celery workers.

The monkey patching is only applied when running as a server (CLI mode),
not when using the Python API (ZnDraw class) in scripts or Jupyter notebooks.
"""
import eventlet

eventlet.monkey_patch()

# Now safe to import CLI and other zndraw modules
from zndraw.cli import app

if __name__ == "__main__":
    app()
