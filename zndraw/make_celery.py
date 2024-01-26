from zndraw.app import ZnDrawApp

with ZnDrawApp(None, False, True, None, None) as flask_app:
    celery_app = flask_app.extensions["celery"]
