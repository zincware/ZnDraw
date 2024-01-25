import contextlib
import threading
import uuid

from celery import Celery, Task
from flask import Flask
from flask_caching import Cache
from flask_socketio import SocketIO

socketio = SocketIO()
cache = Cache(
    config={"CACHE_TYPE": "SimpleCache", "CACHE_DEFAULT_TIMEOUT": 60 * 60 * 24}
)


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


def setup_cache():
    """Setup the cache."""
    cache.set("ROOM_HOSTS", {})
    cache.set("DEFAULT_PYCLIENT", None)

    # dict of {uuid: sid} for each client
    cache.set("pyclients", {})

    # dict of {token: dict}
    cache.set("PER-TOKEN-DATA", {})
    cache.set("MODIFIER", {"default_schema": {}, "active": None, "queue": []})


@contextlib.contextmanager
def create_app(
    use_token, upgrade_insecure_requests, compute_bonds, tutorial: str, auth_token: str
) -> Flask:
    """Create the Flask app."""

    app = Flask(__name__)
    app.config["SECRET_KEY"] = str(uuid.uuid4())

    app.config["TUTORIAL"] = tutorial
    app.config["AUTH_TOKEN"] = auth_token

    if not use_token:  # TODO: handle this differently
        app.config["token"] = "notoken"
    app.config["upgrade_insecure_requests"] = upgrade_insecure_requests
    app.config["compute_bonds"] = compute_bonds

    from .server import main as main_blueprint

    app.register_blueprint(main_blueprint)

    cache.init_app(app)
    socketio.init_app(app, cors_allowed_origins="*")
    setup_cache()

    app.config.from_mapping(
        CELERY=dict(
            broker_url="memory://",
            result_backend="cache",
            cache_backend="memory",
            task_ignore_result=True,
        ),
    )
    app.config.from_prefixed_env()
    celery_app = celery_init_app(app)

    worker = threading.Thread(
        target=celery_app.worker_main, args=(["worker", "--loglevel=info"],)
    )
    worker.start()

    yield app
    celery_app.control.shutdown()
    worker.join()
