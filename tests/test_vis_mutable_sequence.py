import time
from zndraw.zndraw import ZnDraw
import pytest

def test_vis_init_len(s22, server):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert len(vis) == 0 # []

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
    
    # assert vis[0:3] == s22[0:3]
    # assert vis[[0, 2, 4]] == [s22[0], s22[2], s22[4]]
