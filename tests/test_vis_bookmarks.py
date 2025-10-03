import pytest

from zndraw.zndraw import ZnDraw


def test_vis_bookmarks(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    assert vis.bookmarks == {}
    vis.extend(s22)

    vis.bookmarks = {1: "First Frame", 5: "Middle Frame", 9: "Last Frame"}
    del vis[2]
    assert vis.bookmarks == {1: "First Frame", 4: "Middle Frame", 8: "Last Frame"}
    del vis[1]
    assert vis.bookmarks == {3: "Middle Frame", 7: "Last Frame"}
    vis.insert(0, s22[0])
    assert vis.bookmarks == {4: "Middle Frame", 8: "Last Frame"}
    # replacing the frame removes the bookmark as well
    vis[4] = s22[3]
    assert vis.bookmarks == {8: "Last Frame"}
    vis.bookmarks = None
    assert vis.bookmarks == {}


def test_vis_bookmarks_errors(server, s22):
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    with pytest.raises(IndexError):
        vis.bookmarks = {999: "Out of bounds"}

    with pytest.raises(TypeError):
        vis.bookmarks = {0: 123}  # Non-string label
