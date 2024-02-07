import subprocess
import time

import pytest
from ase.build import molecule
from ase.collections import s22
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import uuid

from zndraw import ZnDraw
from zndraw.db.schema import Base
from zndraw.settings import GlobalConfig



@pytest.fixture()
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

    # # Wait a bit for the subprocess to start up
    time.sleep(1)

    yield proc

    # Kill the subprocess after the test runs through
    proc.terminate()
    proc.wait()


def test_zndraw(server):
    vis = ZnDraw(server, token=str(uuid.uuid4()))
    assert vis.socket.call("ping") == "pong"
    assert len(vis) == 0
    vis[0] = molecule("CH4")
    assert len(vis) == 1
    assert vis[0] == molecule("CH4")

def test_zndraw_step(server):
    vis = ZnDraw(server, token=str(uuid.uuid4()))
    assert len(vis) == 0

    vis.extend(list(s22))
    assert len(vis) == len(s22)

    vis.step = 10
    assert vis.step == 10
    vis.step = 5
    assert vis.step == 5

    with pytest.raises(IndexError):
        vis.step = -1
    
    with pytest.raises(IndexError):
        vis.step = len(vis)

def test_zndraw_selection(server):
    vis = ZnDraw(server, token=str(uuid.uuid4()))
    assert len(vis) == 0

    vis.extend(list(s22))
    assert len(vis) == len(s22)

    vis.selection = [0, 1, 2]
    assert vis.selection == [0, 1, 2]
    
    vis.step = 0
    vis[0] = molecule("H2O")
    vis.selection = [0, 1]
    assert vis.selection == [0, 1]
    assert vis[0] == molecule("H2O")
    assert vis.step == 0
    # should we raise en error, if selection is bigger than the structure?
    with pytest.raises(ValueError):
        vis.selection = [2, 2]
    with pytest.raises(ValueError):
        vis.selection = [-1]
    with pytest.raises(ValueError):
        vis.selection = ["data"]
