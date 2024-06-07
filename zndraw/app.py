import dataclasses
import os
import platform
import subprocess
import time
import uuid
import webbrowser

from celery import Celery, Task
from flask import Flask
from flask_socketio import SocketIO
from redis import Redis
from znsocket.exceptions import ConnectionError

from zndraw.settings import CeleryFileSystemConfig
from zndraw.utils import get_port

from .base import FileIO
from .tasks import read_file, run_znsocket_server

socketio = SocketIO()


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


def setup_worker() -> list:
    """Setup the worker."""
    my_env = os.environ.copy()
    if platform.system() == "Darwin" and platform.processor() == "arm":
        # fix celery worker issue on apple silicon
        my_env["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

    worker = subprocess.Popen(
        [
            "celery",
            "-A",
            "zndraw.make_celery",
            "worker",
            "--loglevel=info",
            # "--concurrency=6",
            # "--hostname=default_worker",
            # "--queues=celery",
            # "--prefetch-multiplier=2",
        ],
        env=my_env,
    )
    return [worker]


def create_app() -> Flask:
    """Create the Flask app."""
    app = Flask(__name__)
    from .server import main as main_blueprint

    app.config["SECRET_KEY"] = str(uuid.uuid4())

    app.register_blueprint(main_blueprint)

    socketio.init_app(app, cors_allowed_origins="*")
    if "ZNDRAW_STORAGE" in os.environ and os.environ["ZNDRAW_STORAGE"].startswith(
        "redis"
    ):
        app.config["CELERY"] = {"broker_url": os.environ["ZNDRAW_STORAGE"]}
    else:
        app.config["CELERY"] = CeleryFileSystemConfig().to_dict()

    app.config.from_prefixed_env()
    celery_init_app(app)
    return app


@dataclasses.dataclass
class ZnDrawServer:
    tutorial: str
    auth_token: str
    storage: str
    port: int = 1234
    fileio: FileIO = None
    simgen: bool = False
    celery_worker: bool = True

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
        for worker in self._workers:
            worker.wait()

        # I'm not certain this will really do anything to stop the workers

        self.app.config["redis"].flushall()

    def update_cache(self):
        self.app.config["SECRET_KEY"] = str(uuid.uuid4())
        self.app.config["AUTH_TOKEN"] = self.auth_token
        self.app.config["PORT"] = self.port
        self.app.config["TUTORIAL"] = self.tutorial
        self.app.config["SIMGEN"] = self.simgen

        if self.storage.startswith("redis"):
            self.app.config["redis"] = Redis.from_url(self.storage, decode_responses=True)
        elif self.storage.startswith("znsocket"):
            import znsocket

            for _ in range(100):  # try to connect to znsocket for 10 s
                # if we start znsocket via celery it will take some time to start
                try:
                    self.app.config["redis"] = znsocket.Client.from_url(self.storage)
                    break
                except ConnectionError:
                    time.sleep(0.1)

    def run(self, browser=False):
        if self.celery_worker:
            self._workers = setup_worker()

        if self.storage is None:
            port = get_port(default=3018)
            self.storage = f"znsocket://127.0.0.1:{port}"
            run_znsocket_server.delay(port)

        self.update_cache()

        read_file.delay(
            fileio=self.fileio.to_dict(), io_port=self.port, storage=self.storage
        )

        if browser:
            webbrowser.open(self.url_root)
        socketio.run(self.app, port=self.port, host="0.0.0.0")

    @property
    def url_root(self):
        return f"http://127.0.0.1:{self.port}"
