import pytest

from zndraw.zndraw import ZnDraw


def test_vis_init_len(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert len(vis) == 0  # []


def test_vis_append(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.append(s22[0])
    assert len(vis) == 1
    vis.append(s22[1])
    assert len(vis) == 2


def test_vis_append_invalid(server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    with pytest.raises(TypeError):
        vis.append("not an atoms object")


def test_vis_extend(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    assert len(vis) == len(s22) == 22


def test_vis_extend_invalid(server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    with pytest.raises(TypeError):
        vis.extend(["not", "atoms", "objects"])


def test_vis_getitem(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    for i in range(len(s22)):
        assert vis[i] == s22[i]

    # TODO, needs tests in storage
    # assert vis[0:3] == s22[0:3]
    # assert vis[[0, 2, 4]] == [s22[0], s22[2], s22[4]]
    # assert vis[-1] == s22[-1]


def test_vis_getitem_invalid(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    with pytest.raises(TypeError):
        vis["not an index"]
    with pytest.raises(IndexError):
        vis[100]


def test_vis_setitem(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    vis[0] = s22[1]
    assert vis[0] == s22[1]
    assert s22[0] != s22[1]
    assert len(vis) == len(s22)

    # TODO, needs tests in storage
    # vis[1:3] = s22[3:5]
    # vis[[3,4]] = [s22[5], s22[6]]


def test_vis_setitem_invalid(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    with pytest.raises(TypeError):
        vis[0] = "not an atoms object"
    with pytest.raises(IndexError):
        vis[100] = s22[0]


def test_vis_delitem(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    del vis[0]
    assert len(vis) == len(s22) - 1
    assert vis[0] == s22[1]

    # TODO, needs tests in storage
    # del vis[1:3]
    # del vis[[3,4]]


def test_vis_delitem_invalid(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    with pytest.raises(IndexError):
        del vis[100]

    with pytest.raises(TypeError):
        del vis["not an index"]


def test_vis_insert(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    vis.insert(0, s22[0])
    assert len(vis) == len(s22) + 1
    assert vis[0] == s22[0]
    assert vis[1] == s22[0]
    vis.insert(100, s22[0])
    assert len(vis) == len(s22) + 2
    assert vis[-1] == s22[0]


def test_vis_insert_invalid(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)
    with pytest.raises(TypeError):
        vis.insert(0, "not an atoms object")
