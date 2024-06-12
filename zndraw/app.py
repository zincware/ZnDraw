import eventlet

eventlet.monkey_patch()


import time

import redis
import znsocket.exceptions
from celery import Celery, Task
from flask import Flask
from flask_socketio import SocketIO

from zndraw.server import init_socketio_events, main_blueprint


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


def storage_init_app(app: Flask) -> None:
    if app.config["STORAGE"].startswith("redis"):
        app.extensions["redis"] = redis.Redis.from_url(
            app.config["STORAGE"], decode_responses=True
        )
    elif app.config["STORAGE"].startswith("znsocket"):
        for _ in range(100):  # try to connect to znsocket for 10 s
            # if we start znsocket via celery it will take some time to start
            try:
                app.extensions["redis"] = znsocket.Client.from_url(app.config["STORAGE"])
                break
            except ConnectionError:
                # wait for znsocket to start, if started together with the server
                time.sleep(0.1)
    else:
        raise ValueError(f"Unknown storage type: {app.config['STORAGE']}")


def create_app() -> Flask:
    app = Flask(__name__)
    app.config["SECRET_KEY"] = "secret!"
    # loads all FLASK_ prefixed environment variables into the app config
    app.config.from_prefixed_env()

    # TODO: this will not work without redis!!!!
    if app.config["STORAGE"].startswith("redis"):
        app.config.from_mapping(
            CELERY={
                "broker_url": app.config["STORAGE"],
                "result_backend": app.config["STORAGE"],
                "task_ignore_result": True,
            },
        )
    else:
        raise ValueError(f"Unknown storage type: {app.config['STORAGE']}")

    # Initialize SocketIO
    socketio = SocketIO(app, message_queue=app.config["CELERY"]["broker_url"], cors_allowed_origins="*")

    # Initialize Celery
    celery_init_app(app)

    # Initialize storage
    storage_init_app(app)

    # Register routes and socketio events
    app.register_blueprint(main_blueprint)
    init_socketio_events(socketio)

    # Add socketio to app extensions for easy access
    app.extensions["socketio"] = socketio

    return app
