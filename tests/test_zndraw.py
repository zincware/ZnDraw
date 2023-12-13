# from zndraw.app import create_app


# def test_request_example(client):
#     response = client.app.get("/")
#     assert b"<h2>Hello, World!</h2>" in response.data

import numpy.testing as npt
import pytest
from selenium import webdriver
from selenium.webdriver.chrome.options import Options

from zndraw import ZnDraw

chrome_options = Options()
chrome_options.add_argument("--headless")  # for Chrome >= 109


@pytest.fixture
def vis() -> ZnDraw:
    visualizer = ZnDraw(token="test_token")
    driver = webdriver.Chrome(options=chrome_options)
    driver.get(f"{visualizer.url}/token/{visualizer.token}")
    yield visualizer
    visualizer.close()


@pytest.mark.chrome
def test_gui_running(water):
    vis = ZnDraw(token="test_token")

    driver = webdriver.Chrome(options=chrome_options)
    driver.get(f"{vis.url}/token/{vis.token}")
    assert "ZnDraw" in driver.title
    vis.close()


@pytest.mark.chrome
def test_vis_atoms(vis, water):
    vis[0] = water
    assert len(vis) == 1
    assert vis[0] == water

    assert vis.step == 0


@pytest.mark.chrome
def test_vis_selection(vis, water):
    vis[0] = water
    vis.selection = [1, 2]
    assert vis.selection == [1, 2]


@pytest.mark.parametrize("display_new", [True, False])
@pytest.mark.chrome
def test_vis_step(vis, ase_s22, display_new):
    vis.display_new = display_new
    vis.extend(ase_s22)

    assert len(vis) == 22
    if display_new:
        assert vis.step == 21
    else:
        assert vis.step == 0
    vis.step = 10
    assert vis.step == 10


@pytest.mark.chrome
def test_vis_points(vis, water):
    vis[0] = water
    npt.assert_array_equal(vis.points, [])
    npt.assert_array_equal(vis.segments, [])


@pytest.mark.chrome
def test_multiple_instances(vis, ase_s22):
    vis2 = ZnDraw(url=vis.url, token=vis.token)
    vis2.extend(ase_s22)
    assert len(vis2) == len(vis) == 22

    vis2.close()
