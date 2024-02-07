import threading
import time

import eventlet
import pytest
import socketio

from zndraw.app import create_app
from zndraw.db.schema import Base
from zndraw.utils import get_port
from zndraw.settings import GlobalConfig
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import subprocess


@pytest.fixture(scope="session")
def sio_server():
    port = get_port()

    def run_server(port):
        sio = socketio.Server(cors_allowed_origins="*")
        app = socketio.WSGIApp(sio)

        # react on every event
        @sio.on("*")
        def push_back(event, sid, data):
            sio.emit(event, data, to=sid)

        @sio.on("ping")
        def ping(sid):
            return "pong"

        eventlet.wsgi.server(eventlet.listen(("", port)), app)

    t = threading.Thread(target=run_server, args=(port,), daemon=True)
    t.start()
    time.sleep(1)
    yield f"http://localhost:{port}"


@pytest.fixture(scope="session")
def server():
    port = get_port()

    cmd = [
        "celery",
        "-A",
        "zndraw.make_celery",
        "worker",
        "--loglevel",
        "DEBUG",
        "--queues=io,fast,celery,slow",
    ]
    # proc = subprocess.Popen(cmd) # more verbose for testing
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    config = GlobalConfig.load()
    # TODO: use a temporary path for the database
    # has to be set in the workers somehow

    engine = create_engine(config.database.get_path())
    Base.metadata.create_all(engine)

    def run_server(port):
        app = create_app()

        eventlet.wsgi.server(eventlet.listen(("", port)), app)

    t = threading.Thread(target=run_server, args=(port,), daemon=True)
    t.start()
    time.sleep(1)
    yield f"http://localhost:{port}"
    config.database.unlink()
    proc.terminate()
    proc.wait()
