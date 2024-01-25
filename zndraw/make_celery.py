from zndraw.app import create_app

flask_app = create_app(None, False, True, None, None)
celery_app = flask_app.extensions["celery"]
