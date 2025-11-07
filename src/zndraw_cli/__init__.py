"""CLI entry point package with eventlet monkey patching.

This package provides CLI entry points that apply eventlet monkey patching
BEFORE importing any zndraw modules. This prevents conflicts between eventlet
and Flask/SocketIO/Celery.

The monkey patching is only applied when running CLI commands, not when using
the Python API (ZnDraw class) in scripts or Jupyter notebooks.
"""
import eventlet

# Apply eventlet monkey patching FIRST, before any other imports
eventlet.monkey_patch()
