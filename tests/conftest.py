import multiprocessing as mp
import time

import ase.build
import ase.collections
import pytest
from selenium import webdriver
from selenium.webdriver.chrome.service import Service as ChromeService
from webdriver_manager.chrome import ChromeDriverManager

from zndraw import ZnDraw
from zndraw.app import create_app, socketio
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


@pytest.fixture()
def server():
    port = get_port()

    def run_server():
        app = create_app(None, False, True)
        socketio.run(
            app, port=port, debug=False, host="0.0.0.0"
        )  # NEVER EVER USE  DEBUG=TRUE HERE!!!

    server_proc = mp.Process(
        target=run_server,
    )

    server_proc.start()
    time.sleep(5)
    yield f"http://localhost:{port}"
    server_proc.terminate()
    server_proc.join()


@pytest.fixture()
def vis(server) -> ZnDraw:
    vis = ZnDraw(url=server)
    # vis = ZnDraw(url="https://zndraw.icp.uni-stuttgart.de/", token="18f6d530e74246a9b96cca91f7fc55bc")
    yield vis
    # vis.close()
