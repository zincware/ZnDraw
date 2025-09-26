import logging

import redis
from flask import Flask
from flask_socketio import SocketIO
from celery import Celery, Task

log = logging.getLogger(__name__)

socketio = SocketIO()


def upload_data():
    from zndraw import Client

    c = Client(room="testroom", url="http://localhost:5000")
    c.connect()
    c.append({"x": 10, "y": 10})




def celery_init_app(app: Flask) -> Celery:
    class FlaskTask(Task):
        def __call__(self, *args: object, **kwargs: object) -> object:
            with app.app_context():
                return self.run(*args, **kwargs)

    celery_app = Celery(app.name, task_cls=FlaskTask)
    celery_app.config_from_object(app.config["CELERY"])
    celery_app.set_default()
    app.extensions["celery"] = celery_app
    return celery_app

def create_app(main: bool = False) -> Flask:
    app = Flask(__name__)

    from zndraw.app import main as main_blueprint
    from zndraw.app import tasks

    app.register_blueprint(main_blueprint)

    # Production
    app.config.from_mapping(
        CELERY=dict(
            broker_url="redis://localhost",
            result_backend="redis://localhost",
            task_ignore_result=True,
        ),
    )

    # Standalone / not recommended!
    # app.config.from_mapping(
    #     CELERY=dict(
    #         task_always_eager=True,
    #     ),
    # )

    app.config.from_prefixed_env()
    celery_init_app(app)

    socketio.init_app(app, cors_allowed_origins="*")
    r = redis.Redis(decode_responses=True)

    app.config["SECRET_KEY"] = "your_secret_key"
    app.config["redis"] = r
    
    if main:
        tasks.read_file.delay()

    return app
