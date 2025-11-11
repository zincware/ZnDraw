import logging
from pathlib import Path
from typing import TYPE_CHECKING

import redis
from celery import Celery, Task
from flask import Flask
from flask_socketio import SocketIO
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from prometheus_client import make_wsgi_app
from znsocket import MemoryStorage

if TYPE_CHECKING:
    from zndraw.config import ZnDrawConfig

log = logging.getLogger(__name__)

socketio = SocketIO(cors_allowed_origins="*")


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


def services_init_app(app: Flask) -> None:
    """Initialize service layer with Redis client.

    Services provide domain logic abstraction over Redis operations.
    Must be called after redis_init_app.

    Raises
    ------
    ValueError
        If admin environment variables are misconfigured
    """
    from zndraw.services import (
        AdminService,
        ClientService,
        RoomService,
        SettingsService,
        UserService,
    )

    redis_client = app.extensions["redis"]

    # Initialize AdminService first (validates env vars)
    app.extensions["admin_service"] = AdminService(redis_client)

    # Initialize UserService for user management
    app.extensions["user_service"] = UserService(redis_client)

    app.extensions["client_service"] = ClientService(redis_client)
    app.extensions["room_service"] = RoomService(redis_client)
    app.extensions["settings_service"] = SettingsService(redis_client)


def create_app(
    config: "ZnDrawConfig | None" = None,
    storage_path: str | None = None,
    redis_url: str | None = None,
) -> Flask:
    """Create and configure Flask application.

    Parameters
    ----------
    config : ZnDrawConfig | None
        Configuration object. If None, loads from environment via get_config().
    storage_path : str | None
        Override storage path (for backwards compatibility).
    redis_url : str | None
        Override Redis URL (for backwards compatibility).

    Returns
    -------
    Flask
        Configured Flask application instance.
    """
    from zndraw.config import get_config as _get_config

    # Load config from environment if not provided
    if config is None:
        config = _get_config()

    # Apply overrides for backwards compatibility
    if storage_path is not None:
        config.storage_path = storage_path
    if redis_url is not None:
        config.redis_url = redis_url

    # Configure logging from config
    log_level = getattr(logging, config.log_level.upper(), logging.WARNING)
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    log.info(f"Logging configured at level: {config.log_level}")

    app = Flask(__name__)

    # Store config object in extensions for direct access
    app.extensions["config"] = config

    # Set Flask secret key
    app.config["SECRET_KEY"] = config.flask_secret_key

    from zndraw.app import (
        bookmarks,
        extensions,
        filesystem_bp,
        frames,
        geometries,
        jobs,
        locks,
        media,
        rooms,
        tasks,  # noqa: F401
        utility,
        workers,
    )
    from zndraw.app.file_browser import file_browser as file_browser_blueprint

    app.register_blueprint(utility)
    app.register_blueprint(frames)
    app.register_blueprint(rooms)
    app.register_blueprint(extensions)
    app.register_blueprint(jobs)
    app.register_blueprint(geometries)
    app.register_blueprint(bookmarks)
    app.register_blueprint(media)
    app.register_blueprint(file_browser_blueprint)
    app.register_blueprint(filesystem_bp)
    app.register_blueprint(locks)
    app.register_blueprint(workers)

    # Configure Celery based on Redis availability
    if config.redis_url is None:
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
                broker_url=config.redis_url,
                result_backend=config.redis_url,
                task_ignore_result=True,
            ),
        )

    app.config.from_prefixed_env()
    celery_init_app(app)
    redis_init_app(app, config.redis_url)
    services_init_app(app)

    # Configure SocketIO with Redis message queue for multi-worker support
    if config.redis_url:
        log.info(f"Configuring SocketIO with Redis message queue: {config.redis_url}")
        socketio.init_app(app, message_queue=config.redis_url, cors_allowed_origins="*")
    else:
        log.info("Configuring SocketIO without message queue (single worker mode)")
        socketio.init_app(app, cors_allowed_origins="*")

    app.wsgi_app = DispatcherMiddleware(app.wsgi_app, {
        '/metrics': make_wsgi_app()
    })
    log.info("Prometheus metrics endpoint enabled at /metrics")

    return app
