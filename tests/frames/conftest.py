import eventlet

import random

import pytest
import socketio
import socketio.exceptions

from znsocket import Client, Server


@pytest.fixture
def eventlet_memory_server():
    port = random.randint(10000, 20000)

    def start_server():
        server = Server(port=port)
        server.run()

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

    yield f"znsocket://127.0.0.1:{port}"

    thread.kill()

@pytest.fixture
def znsclient(eventlet_memory_server):
    r = Client.from_url(eventlet_memory_server)
    yield r
    r.flushall()