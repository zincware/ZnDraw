import subprocess
import time

import pytest
from ase.build import molecule
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from zndraw import ZnDraw
from zndraw.db.schema import Base
from zndraw.settings import GlobalConfig


@pytest.fixture
def session() -> sessionmaker:
    """pytest fixture to setup the database"""
    config = GlobalConfig.load()
    # TODO: use a temporary path for the database
    # has to be set in the workers somehow

    engine = create_engine(config.database.get_path())
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)


@pytest.fixture(scope="module")
def run_celery_worker():
    # Start the celery worker subprocess
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

    # Wait a bit for the subprocess to start up
    time.sleep(2)

    yield proc

    # Kill the subprocess after the test runs through
    proc.terminate()
    proc.wait()


def test_zndraw(server, run_celery_worker, session):
    # create a room by calling "server/token/test-room"
    # import requests
    # url = f"{server}/token/test-room"
    # response = requests.get(url)
    # assert response.status_code == 200
    # assert "ZnDraw" in response.text

    # vis = ZnDrawWorker(token="test-room", url=server)
    # vis[0] = molecule("CH4")

    # there is only a room once the webclient connects properly

    vis = ZnDraw(server, token="test-room")
    vis.socket.sleep(1)
    assert vis.socket.call("ping") == "pong"
    assert len(vis) == 0
    vis[0] = molecule("CH4")
    assert len(vis) == 1
    assert vis[0] == molecule("CH4")

    # with pytest.raises(exceptions.RoomNotFound):
    #     assert len(vis) == 0
