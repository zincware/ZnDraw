from zndraw.app import create_app

with create_app(None, False, True, None, None) as flask_app:
    celery_app = flask_app.extensions["celery"]