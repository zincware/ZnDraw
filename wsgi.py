"""WSGI entry point for production deployment with Gunicorn.

This module creates the Flask application instance for Gunicorn.
Configuration is read from environment variables:
- ZNDRAW_STORAGE_PATH: Path to Zarr storage (default: ./zndraw-data.zarr)
- ZNDRAW_REDIS_URL: Redis URL for state management (default: None)
"""

import os

from zndraw.server import create_app, socketio

# Create Flask app with configuration from environment variables
storage_path = os.getenv("ZNDRAW_STORAGE_PATH", "./zndraw-data.zarr")
redis_url = os.getenv("ZNDRAW_REDIS_URL")

app = create_app(storage_path=storage_path, redis_url=redis_url)

# Export both app and socketio for Gunicorn
# Gunicorn will use the 'app' object
__all__ = ["app", "socketio"]
