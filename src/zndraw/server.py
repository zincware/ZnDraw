import logging
import os
from pathlib import Path

import redis
from celery import Celery, Task
from flask import Flask
from flask_socketio import SocketIO
from znsocket import MemoryStorage

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


def redis_init_app(app: Flask, redis_url: str | None) -> redis.Redis | MemoryStorage:
    if redis_url is None:
        r = MemoryStorage()
        app.extensions["redis"] = r
    else:
        r = redis.Redis.from_url(redis_url, decode_responses=True)
        app.extensions["redis"] = r
    return r


def create_app(
    storage_path: str = "./zndraw-data.zarr",
    redis_url: str | None = None,
) -> Flask:
    # Priority: explicit parameter > environment variable > None
    if redis_url is None:
        redis_url = os.getenv("ZNDRAW_REDIS_URL")

    app = Flask(__name__)

    from zndraw.app import main as main_blueprint
    from zndraw.app import tasks  # noqa: F401
    from zndraw.app.file_browser import file_browser as file_browser_blueprint

    app.register_blueprint(main_blueprint)
    app.register_blueprint(file_browser_blueprint)

    # Store configuration
    app.config["STORAGE_PATH"] = storage_path
    app.config["REDIS_URL"] = redis_url

    # Upload configuration
    app.config["UPLOAD_TEMP_DIR"] = os.getenv(
        "ZNDRAW_UPLOAD_TEMP", "/tmp/zndraw_uploads"
    )
    app.config["MAX_CONTENT_LENGTH"] = (
        int(os.getenv("ZNDRAW_MAX_UPLOAD_MB", "500")) * 1024 * 1024
    )

    if redis_url is None:
        data_folder = Path("~/.zincware/zndraw/celery/out").expanduser()
        data_folder_processed = Path("~/.zincware/zndraw/celery/processed").expanduser()
        control_folder = Path("~/.zincware/zndraw/celery/ctrl").expanduser()

        data_folder.mkdir(parents=True, exist_ok=True)
        data_folder_processed.mkdir(parents=True, exist_ok=True)
        control_folder.mkdir(parents=True, exist_ok=True)

        app.config.from_mapping(
            CELERY={
                "broker_url": "filesystem://",
                "result_backend": "cache",
                "cache_backend": "memory",
                "task_ignore_result": True,
                "broker_transport_options": {
                    "data_folder_in": data_folder.as_posix(),
                    "data_folder_out": data_folder.as_posix(),
                    "data_folder_processed": data_folder_processed.as_posix(),
                    "control_folder": control_folder.as_posix(),
                },
            },
        )
    else:
        app.config.from_mapping(
            CELERY=dict(
                broker_url=redis_url,
                result_backend=redis_url,
                task_ignore_result=True,
            ),
        )

    app.config.from_prefixed_env()
    celery_init_app(app)
    redis_init_app(app, redis_url)

    socketio.init_app(app, cors_allowed_origins="*")

    # Configure SECRET_KEY from environment variable or use development default
    app.config["SECRET_KEY"] = os.getenv(
        "FLASK_SECRET_KEY", "dev-secret-key-change-in-production"
    )

    return app
