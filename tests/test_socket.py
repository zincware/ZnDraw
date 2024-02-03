import pytest
from zndraw.app import create_app, socketio
from socketio import Client
import functools

@pytest.fixture()
def app():
    app = create_app()
    app.config.update({
        "TESTING": True,
        "USE_TOKEN": False,
        "SECRET_KEY": "test",
        "upgrade_insecure_requests": False,
        "TUTORIAL": None,
    })

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

def test_socket_connection(sio_client):
    assert sio_client.emit("ping", callback=True) == "pong"
