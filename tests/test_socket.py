import pytest
from socketio import Client

from zndraw.app import create_app, socketio


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


@pytest.fixture()
def server():
    import threading

    import eventlet
    import socketio

    def run_server():
        sio = socketio.Server(cors_allowed_origins="*")
        app = socketio.WSGIApp(sio)

        @sio.on("ping")
        def ping(sid):
            return "pong"

        eventlet.wsgi.server(eventlet.listen(("", 8000)), app)

    t = threading.Thread(target=run_server, daemon=True)
    t.start()
    return


def test_socket_connection_test_client(sio_client):
    assert sio_client.emit("ping", callback=True) == "pong"


def test_socket_connection(server):
    client = Client()
    client.connect("http://localhost:8000")
    assert client.call("ping") == "pong"
