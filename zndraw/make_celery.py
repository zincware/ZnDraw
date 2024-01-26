from zndraw.app import ZnDrawApp

with ZnDrawApp(None, False, True, None, None) as instance:
    celery_app = instance.app.extensions["celery"]
