import multiprocessing as mp
import threading
import time

import ase.build
import ase.collections
import eventlet
import pytest
import socketio
from selenium import webdriver
from selenium.webdriver.chrome.service import Service as ChromeService
from webdriver_manager.chrome import ChromeDriverManager

from zndraw.app import ZnDrawServer
from zndraw.utils import get_port
from zndraw.db.schema import Base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from ase.collections import s22
from zndraw.db import schema
import znframe
from zndraw.app import create_app, socketio

s22 = list(s22)

@pytest.fixture
def water() -> ase.Atoms:
    return ase.build.molecule("H2O")


@pytest.fixture
def ase_s22() -> list[ase.Atoms]:
    return list(ase.collections.s22)


@pytest.fixture()
def setup(request):
    options = webdriver.ChromeOptions()
    options.add_argument("--headless")
    options.add_argument("--disable-3d-apis")
    driver = webdriver.Chrome(
        service=ChromeService(ChromeDriverManager().install()), options=options
    )
    request.cls.driver = driver
    yield request.cls.driver
    request.cls.driver.close()


def run_server(port):
    with ZnDrawServer(None, False, True, None, None, port=port) as app:
        app.run(browser=False)


@pytest.fixture()
def server():
    port = get_port()

    server_proc = mp.Process(target=run_server, args=(port,))

    server_proc.start()
    time.sleep(2)
    try:
        yield f"http://127.0.0.1:{port}"
    finally:
        server_proc.terminate()
        server_proc.join()


@pytest.fixture()
def sio_server():
    port = get_port()

    def run_server(port):

        @socketio.on("*", namespace="/testing")
        def push_back(event, sid, data):
            socketio.emit(event, data, to=sid)

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

        eventlet.wsgi.server(eventlet.listen(("", port)), app)

    t = threading.Thread(target=run_server, args=(port,), daemon=True)
    t.start()
    eventlet.sleep(2)
    return f"http://localhost:{port}"


@pytest.fixture
def session() -> sessionmaker:
    """pytest fixture to setup the database"""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)


@pytest.fixture
def room_session(session) -> sessionmaker:
    with session() as s:
        room = schema.Room(
            token="test_token",
            currentStep=5,
            selection=[1, 2],
            points=[[0, 0, 0], [1, 1, 1]],
        )
        s.add(room)
        for idx, atoms in enumerate(s22):
            frame = schema.Frame(
                data=znframe.Frame.from_atoms(atoms).to_dict(built_in_types=False),
                index=idx,
                room=room,
            )
            s.add(frame)

        # add two bookmarks
        bookmark1 = schema.Bookmark(step=1, text="bm-1", room=room)
        s.add(bookmark1)
        bookmark2 = schema.Bookmark(step=2, text="bm-2", room=room)
        s.add(bookmark2)
        s.commit()

    return session
