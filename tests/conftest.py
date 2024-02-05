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
    t.daemon = True
    t.start()
    time.sleep(1)
    yield f"http://localhost:{port}"