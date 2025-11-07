"""CLI command entry point.

This module imports the main CLI after eventlet has been monkey patched
by the zndraw_cli package __init__.py.
"""
from zndraw.cli import app

if __name__ == "__main__":
    app()
