import eventlet.wsgi
import os
import redis

eventlet.monkey_patch()  # MUST BE THERE FOR THE TESTS TO WORK

import random

import ase
import ase.collections
import pytest
import socketio.exceptions

from zndraw.app import create_app
from zndraw.standalone import run_celery_worker, run_znsocket


@pytest.fixture
def server():
    port = random.randint(10000, 20000)

    def start_server():
        os.environ["FLASK_PORT"] = str(port)
        os.environ["FLASK_STORAGE"] = "redis://localhost:6379/0"

        app = create_app()
        worker = run_celery_worker()

        socketio = app.extensions["socketio"]
        socketio.run(
            app,
            host="0.0.0.0",
            port=app.config["PORT"],
        )

        worker.terminate()


    thread = eventlet.spawn(start_server)

    # wait for the server to be ready
    for _ in range(100):
        try:
            with socketio.SimpleClient() as client:
                client.connect(f"http://localhost:{port}")
                client.disconnect()
                break
        except socketio.exceptions.ConnectionError:
            eventlet.sleep(0.1)
    else:
        raise TimeoutError("Server did not start in time")

    yield f"http://127.0.0.1:{port}"

    r = redis.Redis.from_url("redis://localhost:6379/0")
    r.flushall()

    thread.kill()


@pytest.fixture
def s22() -> list[ase.Atoms]:
    """S22 dataset."""
    return list(ase.collections.s22)
