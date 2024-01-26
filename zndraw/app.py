import dataclasses
import pathlib
import subprocess
import uuid
import webbrowser

from celery import Celery, Task
from flask import Flask
from flask_caching import Cache
from flask_socketio import SocketIO

socketio = SocketIO()


def get_cache():
    # read config for cache from zndraw config
    cache = Cache(
        config={
            "CACHE_TYPE": "FileSystemCache",
            "CACHE_DEFAULT_TIMEOUT": 60 * 60 * 24,
            "CACHE_DIR": ".zndraw/cache",
        }
    )
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
    cache.set("DEFAULT_PYCLIENT", None)
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


@dataclasses.dataclass
class ZnDrawApp:
    use_token: bool
    upgrade_insecure_requests: bool
    compute_bonds: bool
    tutorial: str
    auth_token: str
    port: int = 1234
    fileio: FileIO = None

    _workers: list = None

    def __enter__(self):
        """Create the Flask app."""
        # with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path.cwd() / ".zndraw"
        celery_data_folder_in = tmpdir / "celery" / "out"
        celery_data_folder_in.mkdir(parents=True, exist_ok=True)

        celery_data_folder_out = tmpdir / "celery" / "out"
        celery_data_folder_out.mkdir(parents=True, exist_ok=True)

        celery_data_folder_processed = tmpdir / "celery" / "processed"
        celery_data_folder_processed.mkdir(parents=True, exist_ok=True)

        cachdir = tmpdir / "cache"
        cachdir.mkdir(parents=True, exist_ok=True)

        app = Flask(__name__)

        from .server import main as main_blueprint

        app.register_blueprint(main_blueprint)

        cache.init_app(app)
        socketio.init_app(app, cors_allowed_origins="*")

        app.config.from_mapping(
            CELERY=dict(
                broker_url="filesystem://",
                broker_transport_options=dict(
                    data_folder_in=celery_data_folder_in,
                    data_folder_out=celery_data_folder_out,
                    data_folder_processed=celery_data_folder_processed,
                ),
                result_backend="cache",
                cache_backend="memory",
                task_ignore_result=True,
            ),
        )
        app.config.from_prefixed_env()
        celery_app = celery_init_app(app)
        self.app = app

        print(f"file io from cache: {cache.get('FILEIO')}")

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._workers is None:
            return
        for worker in self._workers:
            worker.terminate()
        cache.clear()
        for worker in self._workers:
            worker.wait()

        # remove tmpdir, but only if this is the main thread
        # and not a worker thread of celery that e.g. has been restarted

    def update_cache(self):
        self.app.config["SECRET_KEY"] = str(uuid.uuid4())

        self.app.config["TUTORIAL"] = self.tutorial
        self.app.config["AUTH_TOKEN"] = self.auth_token

        if not self.use_token:  # TODO: handle this differently
            self.app.config["token"] = "notoken"
        self.app.config["upgrade_insecure_requests"] = self.upgrade_insecure_requests
        self.app.config["compute_bonds"] = self.compute_bonds

        setup_cache()
        self.fileio.name = "This is a test"
        cache.set("FILEIO", self.fileio)

    def run(self, browser=False):
        self.update_cache()
        self._workers = setup_worker()
        if browser:
            webbrowser.open(self.url_root)
        socketio.run(self.app, port=self.port, host="0.0.0.0")

    @property
    def url_root(self):
        return f"http://127.0.0.1:{self.port}"
