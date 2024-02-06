import threading
import time

import eventlet
import pytest
import socketio

from zndraw.app import create_app
from zndraw.utils import get_port


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

    def run_server(port):
        app = create_app()

        eventlet.wsgi.server(eventlet.listen(("", port)), app)

    t = threading.Thread(target=run_server, args=(port,), daemon=True)
    t.start()
    time.sleep(1)
    yield f"http://localhost:{port}"
