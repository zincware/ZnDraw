import os

import eventlet.wsgi

eventlet.monkey_patch()  # MUST BE THERE FOR THE TESTS TO WORK

import random
import subprocess

import ase.build
import ase.collections
import pytest
import socketio.exceptions

from zndraw.app import create_app


@pytest.fixture
def server():
    port = random.randint(10000, 20000)

    os.environ["FLASK_PORT"] = str(port)
    os.environ["FLASK_STORAGE"] = "redis://localhost:6379/0"
    os.environ["FLASK_SERVER_URL"] = f"http://localhost:{port}"

    proc = subprocess.Popen(
        [
            "celery",
            "-A",
            "zndraw.make_celery.celery_app",
            "worker",
            "--loglevel=info",
            "-P",
            "eventlet",
        ]
    )

    def start_server():
        os.environ["FLASK_PORT"] = str(port)
        os.environ["FLASK_STORAGE"] = "redis://localhost:6379/0"

        app = create_app()
        app.config["TESTING"] = True

        socketio = app.extensions["socketio"]
        try:
            socketio.run(
                app,
                host="0.0.0.0",
                port=app.config["PORT"],
            )
        finally:
            app.extensions["redis"].flushall()

    thread = eventlet.spawn(start_server)

    # wait for the server to be ready
    for _ in range(100):
        try:
            with socketio.SimpleClient() as client:
                client.connect(f"http://localhost:{port}")
                break
        except socketio.exceptions.ConnectionError:
            eventlet.sleep(0.1)
    else:
        raise TimeoutError("Server did not start in time")

    yield f"http://127.0.0.1:{port}"

    thread.kill()
    proc.kill()
    proc.wait()


@pytest.fixture
def s22() -> list[ase.Atoms]:
    """S22 dataset."""
    return list(ase.collections.s22)


@pytest.fixture
def water() -> ase.Atoms:
    """Water molecule."""
    return ase.build.molecule("H2O")
