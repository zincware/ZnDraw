import pathlib
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
    app.config.update(SESSION_COOKIE_SAMESITE="None", SESSION_COOKIE_SECURE=True)
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
        # nothing else supported, using filesystem storage
        data_folder = pathlib.Path("~/.zincware/zndraw/celery/out").expanduser()
        data_folder_processed = pathlib.Path(
            "~/.zincware/zndraw/celery/processed"
        ).expanduser()
        control_folder = pathlib.Path("~/.zincware/zndraw/celery/ctrl").expanduser()

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

    # Initialize SocketIO
    message_queue = (
        app.config["CELERY"]["broker_url"]
        if app.config["STORAGE"].startswith("redis")
        else None
    )
    max_http_buffer_size = app.config.get("MAX_HTTP_BUFFER_SIZE")
    kwargs = {
        "message_queue": message_queue,
        "cors_allowed_origins": "*",
    }
    if "SOCKETIO_PING_TIMEOUT" in app.config:
        kwargs["ping_timeout"] = app.config["SOCKETIO_PING_TIMEOUT"]
    if max_http_buffer_size is not None:
        kwargs["max_http_buffer_size"] = int(max_http_buffer_size)

    socketio = SocketIO(app, **kwargs, logger=False, engineio_logger=False)

    # Initialize Celery
    celery_init_app(app)

    # Initialize storage
    storage_init_app(app)
    # we only need this server if we are using redis
    # otherwise a znsocket server will run anyhow
    if app.config["STORAGE"].startswith("redis"):
        # TODO: if we run standalone with znsocket running, don't start its own server but attach as well!
        from redis import Redis

        znsocket.attach_events(
            socketio.server,
            storage=Redis.from_url(app.config["STORAGE"], decode_responses=True),
        )
    else:
        znsocket.attach_events(
            socketio.server,
            storage=znsocket.Client.from_url(
                app.config["STORAGE"], decode_responses=True
            ),
        )

    # Register routes and socketio events
    app.register_blueprint(main_blueprint)
    init_socketio_events(socketio)

    # Add socketio to app extensions for easy access
    app.extensions["socketio"] = socketio

    return app
