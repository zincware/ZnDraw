import logging

import redis
from celery import Celery, Task
from flask import Flask
from flask_socketio import SocketIO

log = logging.getLogger(__name__)

socketio = SocketIO(cors_allowed_origins="*")


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


def redis_init_app(app: Flask, redis_url: str) -> redis.Redis:
    r = redis.Redis.from_url(redis_url, decode_responses=True)
    app.extensions["redis"] = r
    return r


def create_app(
    storage_path: str = "./zndraw-data.zarr",
    redis_url: str = "redis://localhost:6379",
) -> Flask:
    app = Flask(__name__)

    from zndraw.app import main as main_blueprint
    from zndraw.app import tasks  # noqa: F401

    app.register_blueprint(main_blueprint)

    # Store configuration
    app.config["STORAGE_PATH"] = storage_path
    app.config["REDIS_URL"] = redis_url

    # Production
    app.config.from_mapping(
        CELERY=dict(
            broker_url=redis_url,
            result_backend=redis_url,
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
    redis_init_app(app, redis_url)

    socketio.init_app(app, cors_allowed_origins="*")

    app.config["SECRET_KEY"] = "your_secret_key"

    return app
