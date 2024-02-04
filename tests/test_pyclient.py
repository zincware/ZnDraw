from unittest import mock
from zndraw.zndraw_worker import ZnDrawWorker
from zndraw import ZnDraw

import pytest
from socketio import Client
import eventlet

from zndraw.app import create_app, socketio
from zndraw.zndraw_worker import ZnDrawWorker
from zndraw.data import RoomGetData
from dataclasses import asdict


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

@pytest.fixture
def spawn_celery_worker():
    from zndraw.app import setup_worker
    workers = setup_worker()
    yield
    for worker in workers:
        worker.kill()
    for worker in workers:
        worker.wait()


def test_del_atoms(room_session, sio_client, spawn_celery_worker):
    with mock.patch("zndraw.zndraw_worker.Session", room_session):
        sio_client.emit("room:get", asdict(RoomGetData(step=True)))
        eventlet.sleep(5)
        
        answer = sio_client.get_received()
        assert len(answer) == 1
        assert answer[0] == "WTF"
        assert answer[0] == {'args': ['HelloWorld'], 'name': 'room:get', 'namespace': '/'}

