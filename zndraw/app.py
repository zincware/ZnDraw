import dataclasses
import pathlib
import subprocess
import uuid
import webbrowser

from celery import Celery, Task
from flask import Flask
from flask_socketio import SocketIO
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from zndraw.db.schema import Base

from .settings import GlobalConfig

socketio = SocketIO()


@dataclasses.dataclass
class FileIO:
    name: str = None
    start: int = 0
    stop: int = None
    step: int = 1
    remote: str = None
    rev: str = None

    def to_dict(self):
        return dataclasses.asdict(self)


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


def setup_worker(silence: bool) -> list:
    """Setup the worker."""
    import os
    import platform

    my_env = os.environ.copy()
    if platform.system() == "Darwin" and platform.processor() == "arm":
        # fix celery worker issue on apple silicon
        my_env["OBJC_DISABLE_INITIALIZE_FORK_SAFETY"] = "YES"

    fast_worker = subprocess.Popen(
        [
            "celery",
            "-A",
            "zndraw.make_celery",
            "worker",
            "--loglevel=info",
            "--concurrency=6",
            "--hostname=fast_worker",
            "--queues=fast",
            "--prefetch-multiplier=20",
        ],
        env=my_env,
        stdout=subprocess.PIPE if silence else None,
        stderr=subprocess.PIPE if silence else None,
    )

    default_worker = subprocess.Popen(
        [
            "celery",
            "-A",
            "zndraw.make_celery",
            "worker",
            "--loglevel=info",
            "--concurrency=6",
            "--hostname=default_worker",
            "--queues=celery",
            "--prefetch-multiplier=2",
        ],
        env=my_env,
        stdout=subprocess.PIPE if silence else None,
        stderr=subprocess.PIPE if silence else None,
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
        ],
        env=my_env,
        stdout=subprocess.PIPE if silence else None,
        stderr=subprocess.PIPE if silence else None,
    )
    return [fast_worker, default_worker, slow_worker]


def create_app() -> Flask:
    """Create the Flask app."""

    app = Flask(__name__)
    from .server import main as main_blueprint

    # values to be overwritten by the server
    # these are the default
    app.config["SECRET_KEY"] = str(uuid.uuid4())

    app.config["TUTORIAL"] = ""
    app.config["AUTH_TOKEN"] = ""
    app.config["USE_TOKEN"] = True
    app.config["PORT"] = 1234
    # where is the port used?

    app.config["upgrade_insecure_requests"] = False
    app.config["compute_bonds"] = True

    app.register_blueprint(main_blueprint)

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
    simgen: bool = False

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

    def update_cache(self):
        self.app.config["SECRET_KEY"] = str(uuid.uuid4())

        self.app.config["TUTORIAL"] = self.tutorial
        self.app.config["AUTH_TOKEN"] = self.auth_token
        self.app.config["USE_TOKEN"] = self.use_token
        self.app.config["PORT"] = self.port
        self.app.config["SIMGEN"] = self.simgen

        self.app.config["upgrade_insecure_requests"] = self.upgrade_insecure_requests
        self.app.config["compute_bonds"] = self.compute_bonds

        self.app.config["FileIO"] = self.fileio.to_dict() if self.fileio else {}

    def run(self, browser=False):
        self.update_cache()
        self._workers = setup_worker(silence=not self.use_token)

        config = GlobalConfig.load()
        try:
            pathlib.Path(config.database.path).expanduser().unlink(missing_ok=True)
        except AttributeError:
            pass  # only for sqlite config
        engine = create_engine(config.database.get_path())
        Base.metadata.create_all(engine)
        session = sessionmaker(bind=engine)
        self._purge_old_modifier_clients(session)
        self._mark_old_queue_items_as_failed(session)
        if browser:
            webbrowser.open(self.url_root)
        socketio.run(self.app, port=self.port, host="0.0.0.0")

    @property
    def url_root(self):
        return f"http://127.0.0.1:{self.port}"

    @staticmethod
    def _purge_old_modifier_clients(session):
        from zndraw.db.schema import ModifierClient

        with session() as ses:
            ses.query(ModifierClient).delete()
            ses.commit()

    @staticmethod
    def _mark_old_queue_items_as_failed(session):
        from zndraw.db.schema import QueueItem

        with session() as ses:
            ses.query(QueueItem).filter(QueueItem.status == "queued").update(
                dict(status="failed:server_restart")
            )
            ses.commit()
