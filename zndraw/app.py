import contextlib
import pathlib
import subprocess
import uuid
import dataclasses

from celery import Celery, Task
from flask import Flask
from flask_caching import Cache
from flask_socketio import SocketIO

socketio = SocketIO()
cache = Cache(
    config={
        "CACHE_TYPE": "FileSystemCache",
        "CACHE_DEFAULT_TIMEOUT": 60 * 60 * 24,
        "CACHE_DIR": ".zndraw/cache",
    }
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

    workers: list = dataclasses.field(init=False, default_factory=setup_worker)

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
        app.config["SECRET_KEY"] = str(uuid.uuid4())

        app.config["TUTORIAL"] = self.tutorial
        app.config["AUTH_TOKEN"] = self.auth_token

        if not self.use_token:  # TODO: handle this differently
            app.config["token"] = "notoken"
        app.config["upgrade_insecure_requests"] = self.upgrade_insecure_requests
        app.config["compute_bonds"] = self.compute_bonds

        from .server import main as main_blueprint

        app.register_blueprint(main_blueprint)

        cache.init_app(app)
        socketio.init_app(app, cors_allowed_origins="*")
        setup_cache()

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

        return app

    def __exit__(self, exc_type, exc_value, traceback): 
        for worker in self.workers:
            worker.terminate()
        cache.clear()
        for worker in self.workers:
            worker.wait()

        # remove tmpdir, but only if this is the main thread
        # and not a worker thread of celery that e.g. has been restarted
