"""WSGI entry point for production deployment with Gunicorn.

This module creates the Flask application instance for Gunicorn.
Configuration is read from environment variables.

Uses threaded workers (sync worker class) to avoid event loop conflicts
with Zarr v3 synchronous operations. Redis message queue coordinates
SocketIO communication across containers.

Environment Variables:
- ZNDRAW_STORAGE_PATH: Path to Zarr storage (default: ./zndraw-data.zarr)
- ZNDRAW_REDIS_URL: Redis URL for state management (default: None)
- ZNDRAW_SIMGEN_ENABLED: Enable SiMGen molecular generation (default: false)
- FLASK_SECRET_KEY: Flask secret key for sessions (default: dev key)
- ZNDRAW_FILE_BROWSER_ENABLED: Enable file browser (default: false)
- ZNDRAW_FILE_BROWSER_ROOT: File browser root directory (default: cwd)
"""
import eventlet
eventlet.monkey_patch()

import os

from zndraw.server import create_app, socketio

# Create Flask app with configuration from environment variables
storage_path = os.getenv("ZNDRAW_STORAGE_PATH", "./zndraw-data.zarr")
redis_url = os.getenv("ZNDRAW_REDIS_URL")

app = create_app(storage_path=storage_path, redis_url=redis_url)

# Apply additional configuration from environment variables
# These are normally set by the CLI but need to be applied in WSGI mode
app.config["SIMGEN_ENABLED"] = os.getenv("ZNDRAW_SIMGEN_ENABLED", "false").lower() in ("true", "1", "yes")
app.config["SERVER_URL"] = os.getenv("ZNDRAW_SERVER_URL", "http://localhost:5000")
app.config["FILE_BROWSER_ENABLED"] = os.getenv("ZNDRAW_FILE_BROWSER_ENABLED", "false").lower() in ("true", "1", "yes")
app.config["FILE_BROWSER_ROOT"] = os.getenv("ZNDRAW_FILE_BROWSER_ROOT", os.getcwd())

# Export both app and socketio for Gunicorn
# Gunicorn will use the 'app' object
__all__ = ["app", "socketio"]
