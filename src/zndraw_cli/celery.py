"""Celery worker entry point.

This module provides the Celery app after eventlet has been monkey patched
by the zndraw_cli package __init__.py.
"""
from zndraw.app.make_celery import celery_app

# Export celery_app for celery command line
__all__ = ["celery_app"]
