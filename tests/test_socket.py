import pytest
from socketio import Client

from zndraw.app import create_app, socketio
from zndraw.zndraw_worker import ZnDrawWorker


@pytest.fixture()
def app():
    app = create_app()
    app.config.update(
        {
            "TESTING": True,
            "USE_TOKEN": False,
            "SECRET_KEY": "test",
            "upgrade_insecure_requests": False,
            "TUTORIAL": None,
        }
    )

    # other setup can go here

    yield app

    # clean up / reset resources here


@pytest.fixture()
def client(app):
    return app.test_client()


@pytest.fixture()
def sio_client(app):
    return socketio.test_client(app)


@pytest.fixture()
def runner(app):
    return app.test_cli_runner()


def test_socket_connection_test_client(sio_client):
    assert sio_client.emit("ping", callback=True) == "pong"


def test_socket_connection(sio_server):
    client = Client()
    client.connect(sio_server)
    assert client.call("ping") == "pong"


def test_zndraw_worker_log(sio_server):
    worker = ZnDrawWorker(token="test_token", url=sio_server)
    assert worker.socket.call("ping") == "pong"
    global answer
    answer = None

    def on_answer(data):
        global answer
        answer = data

    worker.socket.on("message:log", on_answer)
    worker.log("test")
    while answer is None:
        worker.socket.sleep(0.1)
    assert answer == {"message": "test", "token": "test_token"}
