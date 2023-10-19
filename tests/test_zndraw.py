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

    assert len(vis) == 22

    vis2.close()

@pytest.mark.chrome
def test_delete_backwards(vis, ase_s22):
    vis.extend(ase_s22)
    
    assert len(vis) == 22
    for idx in range(22):
        assert vis[idx] == ase_s22[idx]
    

    for idx in range(22):
        del vis[len(vis) - 1]
        assert len(vis) == 22 - idx - 1
        for jdx in range(len(vis)):
            assert vis[jdx] == ase_s22[jdx]

@pytest.mark.chrome
def test_delete_forwards(vis, ase_s22):
    vis.extend(ase_s22)
    
    assert len(vis) == 22
    for idx in range(22):
        assert vis[idx] == ase_s22[idx]
    

    for idx in range(22):
        del vis[0]
        assert len(vis) == 22 - idx - 1
        for jdx in range(len(vis)):
            assert vis[jdx] == ase_s22[jdx + idx + 1]

@pytest.mark.chrome
def test_delete_middle(vis, ase_s22):
    vis.extend(ase_s22)
    
    assert len(vis) == 22
    for idx in range(22):
        assert vis[idx] == ase_s22[idx]
    
    for idx in range(0, 22):
        if len(vis) <= 10:
            with pytest.raises(IndexError):
                del vis[10]
        else:
            del vis[10]
            assert len(vis) == 22 - idx - 1
            for jdx in range(len(vis)):
                assert vis[jdx] == ase_s22[jdx] if jdx < 10 else ase_s22[jdx + 1]

@pytest.mark.chrome
def test_delete_slice_backwards(vis, ase_s22):
    vis.extend(ase_s22)
    
    assert len(vis) == 22
    for idx in range(22):
        assert vis[idx] == ase_s22[idx]
    

    for idx in range(22):
        del vis[len(vis) - 1:]
        assert len(vis) == 22 - idx - 1
        for jdx in range(len(vis)):
            assert vis[jdx] == ase_s22[jdx]

@pytest.mark.chrome
def test_delete_slice_forwards(vis, ase_s22):
    vis.extend(ase_s22)
    
    assert len(vis) == 22
    for idx in range(22):
        assert vis[idx] == ase_s22[idx]
    

    for idx in range(22):
        del vis[:1]
        assert len(vis) == 22 - idx - 1
        for jdx in range(len(vis)):
            assert vis[jdx] == ase_s22[jdx + idx + 1]

@pytest.mark.chrome
def test_delete_slice_middle(vis, ase_s22):
    vis.extend(ase_s22)
    
    assert len(vis) == 22
    for idx in range(22):
        assert vis[idx] == ase_s22[idx]
    
    for idx in range(0, 22):
        del vis[10:11]
        if len(vis) <= 10:
            assert len(vis) == 10
            for jdx in range(len(vis)):
                assert vis[jdx] == ase_s22[jdx] if jdx < 10 else ase_s22[jdx + 1]
        else:
            assert len(vis) == 22 - idx - 1
            for jdx in range(len(vis)):
                assert vis[jdx] == ase_s22[jdx] if jdx < 10 else ase_s22[jdx + 1]
