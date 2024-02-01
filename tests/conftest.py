import multiprocessing as mp
import time

import ase.build
import ase.collections
import pytest
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
    time.sleep(1)
    try:
        yield f"http://127.0.0.1:{port}"
    finally:
        server_proc.terminate()
        server_proc.join()
