from zndraw import tasks  # noqa used for registering tasks at the moment
from zndraw.app import create_app

flask_app = create_app()
celery_app = flask_app.extensions["celery"]
