import time

import pytest

from zndraw.zndraw import ZnDraw


def test_vis_init_selection(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert vis.selection == ()


def test_vis_set_selection(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    vis.selection = [0, 1, 2]
    assert vis.selection == (0, 1, 2)
    vis.selection = []
    assert vis.selection == ()
    vis.selection = [0, 1, 2]
    assert vis.selection == (0, 1, 2)
    vis.selection = None
    assert vis.selection == ()


def test_vis_set_selection_invalid(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    with pytest.raises(ValueError):
        vis.selection = [100]
    with pytest.raises(ValueError):
        vis.selection = [-1]
    with pytest.raises(ValueError):
        vis.selection = "not a list"
    with pytest.raises(ValueError):
        vis.selection = [0, "a"]


def test_vis_selection_persistence(s22, server):
    vis1 = ZnDraw(url=server, room="testroom", user="testuser1")
    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")
    vis1.extend(s22)
    vis2.extend(s22)
    vis1.selection = [0, 1, 2]
    time.sleep(0.1)
    assert vis1.selection == (0, 1, 2)
    assert vis2.selection == (0, 1, 2)
    vis2.selection = [3, 4]
    time.sleep(0.1)
    assert vis1.selection == (3, 4)
    assert vis2.selection == (3, 4)


def test_vis_selection_different_rooms(s22, server):
    vis1 = ZnDraw(url=server, room="room1", user="testuser1")
    vis2 = ZnDraw(url=server, room="room2", user="testuser2")
    vis1.extend(s22)
    vis2.extend(s22)
    vis1.selection = [0, 1]
    vis2.selection = [2, 3]
    assert vis1.selection == (0, 1)
    assert vis2.selection == (2, 3)


def test_vis_selection_on_join(s22, server):
    vis1 = ZnDraw(url=server, room="testroom", user="testuser1")
    vis1.extend(s22)
    vis1.selection = [0, 1, 2]
    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")
    time.sleep(0.1)
    assert vis2.selection == (0, 1, 2)


# same test but for frame_selection


def test_vis_init_frame_selection(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert vis.frame_selection == ()


def test_vis_set_frame_selection(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    vis.frame_selection = [0, 1, 2]
    assert vis.frame_selection == (0, 1, 2)
    vis.frame_selection = []
    assert vis.frame_selection == ()
    vis.frame_selection = [0, 1, 2]
    assert vis.frame_selection == (0, 1, 2)
    vis.frame_selection = None
    assert vis.frame_selection == ()


def test_vis_set_frame_selection_invalid(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    with pytest.raises(ValueError):
        vis.frame_selection = [100]
    with pytest.raises(ValueError):
        vis.frame_selection = [-1]
    with pytest.raises(ValueError):
        vis.frame_selection = "not a list"
    with pytest.raises(ValueError):
        vis.frame_selection = [0, "a"]


def test_vis_frame_selection_persistence(s22, server):
    vis1 = ZnDraw(url=server, room="testroom", user="testuser1")
    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")
    vis1.extend(s22)
    vis2.extend(s22)
    vis1.frame_selection = [0, 1, 2]
    time.sleep(0.1)
    assert vis1.frame_selection == (0, 1, 2)
    assert vis2.frame_selection == (0, 1, 2)
    vis2.frame_selection = [3, 4]
    time.sleep(0.1)
    assert vis1.frame_selection == (3, 4)
    assert vis2.frame_selection == (3, 4)


def test_vis_frame_selection_different_rooms(s22, server):
    vis1 = ZnDraw(url=server, room="room1", user="testuser1")
    vis2 = ZnDraw(url=server, room="room2", user="testuser2")
    vis1.extend(s22)
    vis2.extend(s22)
    vis1.frame_selection = [0, 1]
    vis2.frame_selection = [2, 3]
    assert vis1.frame_selection == (0, 1)
    assert vis2.frame_selection == (2, 3)


def test_vis_frame_selection_on_join(s22, server):
    vis1 = ZnDraw(url=server, room="testroom", user="testuser1")
    vis1.extend(s22)
    vis1.frame_selection = [0, 1, 2]
    vis2 = ZnDraw(url=server, room="testroom", user="testuser2")
    time.sleep(0.1)
    assert vis2.frame_selection == (0, 1, 2)
