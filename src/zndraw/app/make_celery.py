"""Celery application factory.

This module is loaded by celery workers using --pool=eventlet.
Monkey patching MUST happen before any other imports.
"""

import eventlet

eventlet.monkey_patch()

# Now safe to import Flask app
from zndraw import create_app  # noqa: E402

flask_app = create_app()
celery_app = flask_app.extensions["celery"]
