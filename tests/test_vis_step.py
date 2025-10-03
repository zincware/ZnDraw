import pytest

from zndraw.zndraw import ZnDraw


def test_vis_step(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert vis.bookmarks == {}
    vis.extend(s22)
    assert vis.step == 0

    vis.step = 5
    assert vis.step == 5

    with pytest.raises(ValueError):
        vis.step = -1

    with pytest.raises(ValueError):
        vis.step = 100


def test_vis_step_two_users(server, s22):
    vis1 = ZnDraw(url=server, room="testroom", user="user1")
    vis2 = ZnDraw(url=server, room="testroom", user="user2")

    vis1.extend(s22)

    assert vis1.step == 0
    assert vis2.step == 0

    vis1.step = 3
    vis1.socket.sio.sleep(0.1)
    assert vis1.step == 3
    assert vis2.step == 3

    vis2.step = 5
    vis1.socket.sio.sleep(0.1)
    assert vis1.step == 5
    assert vis2.step == 5

    with pytest.raises(ValueError):
        vis1.step = 100
    with pytest.raises(ValueError):
        vis2.step = -1

    assert vis1.step == 5
