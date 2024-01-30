import dataclasses
import subprocess
import uuid
import webbrowser

from celery import Celery, Task
from flask import Flask
from flask_caching import Cache
from flask_socketio import SocketIO

from .settings import GlobalConfig

socketio = SocketIO()


def get_cache():
    # read config for cache from zndraw config
    config = GlobalConfig.load()
    cache = Cache(config=config.cache.to_dict())
    return cache


cache = get_cache()


@dataclasses.dataclass
class FileIO:
    name: str = None
    start: int = 0
    stop: int = None
    step: int = 1
    remote: str = None
    rev: str = None


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
    cache.set("TEST", "Hello World!")

    # dict of {uuid: sid} for each client
    cache.set("pyclients", {})

    # dict of {token: dict}
    cache.set("PER-TOKEN-DATA", {})
    cache.set("MODIFIER", {"default_schema": {}, "active": None, "queue": []})


def setup_worker() -> list:
    """Setup the worker."""
    fast_worker = subprocess.Popen(
        [
            "celery",
            "-A",
            "zndraw.make_celery",
            "worker",
            "--loglevel=info",
            "--concurrency=16",
            "--hostname=fast_worker",
            "--queues=fast,celery",
        ]
    )

    slow_worker = subprocess.Popen(
        [
            "celery",
            "-A",
            "zndraw.make_celery",
            "worker",
            "--loglevel=info",
            "--concurrency=1",
            "--hostname=slow_worker",
            "--queues=slow",
        ]
    )
    return [fast_worker, slow_worker]


def create_app() -> Flask:
    """Create the Flask app."""

    app = Flask(__name__)
    from .server import main as main_blueprint

    app.register_blueprint(main_blueprint)

    cache.init_app(app)
    socketio.init_app(app, cors_allowed_origins="*")

    app.config.from_mapping(
        CELERY=GlobalConfig.load().celery.to_dict(),
    )
    app.config.from_prefixed_env()
    celery_init_app(app)
    return app


@dataclasses.dataclass
class ZnDrawServer:
    use_token: bool
    upgrade_insecure_requests: bool
    compute_bonds: bool
    tutorial: str
    auth_token: str
    port: int = 1234
    fileio: FileIO = None

    _workers: list = None
    app: Flask = None

    def __enter__(self):
        self.app = create_app()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._workers is None:
            return
        for worker in self._workers:
            worker.kill()
        cache.clear()
        for worker in self._workers:
            worker.wait()

    def update_cache(self):
        self.app.config["SECRET_KEY"] = str(uuid.uuid4())

        self.app.config["TUTORIAL"] = self.tutorial
        self.app.config["AUTH_TOKEN"] = self.auth_token
        self.app.config["USE_TOKEN"] = self.use_token
        self.app.config["PORT"] = self.port

        self.app.config["upgrade_insecure_requests"] = self.upgrade_insecure_requests
        self.app.config["compute_bonds"] = self.compute_bonds

        setup_cache()
        cache.set("FILEIO", self.fileio)
        # Mapping of modifier name f"name" to [request.sid, ...]
        cache.set("ROOM_MODIFIER_HOSTS", {})
        cache.set("MODIFIER_HOSTS", {})
        # Mapping of modifier name f"name" to [dict, ...]
        cache.set("ROOM_MODIFIER_SCHEMA", {})
        cache.set("MODIFIER_SCHEMA", {})
        # Keep track of available modifiers (bugfix for not being able to call("active")) due
        # to connection timeouts. This is a mapping {sid: bool} where bool is True if the
        # modifier is available and False otherwise (running something at the moment).
        cache.set("MODIFIER_AVAILABLE", {})
        # timeout per registered SID (not per modifier) in seconds
        cache.set("MODIFIER_TIMEOUT", {})

    def run(self, browser=False):
        self.update_cache()
        self._workers = setup_worker()
        if browser:
            webbrowser.open(self.url_root)
        socketio.run(self.app, port=self.port, host="0.0.0.0")

    @property
    def url_root(self):
        return f"http://127.0.0.1:{self.port}"
