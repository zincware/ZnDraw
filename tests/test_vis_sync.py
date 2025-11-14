from zndraw import ZnDraw
from zndraw.geometries import Box


def test_step_sync(server, s22):
    vis1 = ZnDraw(url=server, room="testroom", user="u1")
    vis1.extend(s22)
    assert vis1.step == 0

    vis1.step = 7
    assert vis1.step == 7
    
    vis2 = ZnDraw(url=server, room="testroom", user="u2")
    assert vis2.step == 7

def test_settings_sync(server, s22):
    """Test that settings are per-user (not shared across users).

    Each user has their own settings in the same room.
    """
    vis1 = ZnDraw(url=server, room="testroom", user="u1")
    vis1.extend(s22)
    vis1.settings.camera.far_plane = 100.0
    assert vis1.settings.camera.far_plane == 100.0

    # vis2 is a different user, should have default settings
    vis2 = ZnDraw(url=server, room="testroom", user="u2")
    assert vis2.settings.camera.far_plane == 300  # Default value, not synced from u1

def test_bookmarks_sync(server, s22):
    vis1 = ZnDraw(url=server, room="testroom", user="u1")
    vis1.extend(s22)
    assert vis1.bookmarks == {}

    vis1.bookmarks[5] = "My Bookmark"
    assert vis1.bookmarks == {5: "My Bookmark"}

    vis2 = ZnDraw(url=server, room="testroom", user="u2")
    assert vis2.bookmarks == {5: "My Bookmark"}

def test_geometries_sync(server, s22):
    vis1 = ZnDraw(url=server, room="testroom", user="u1")
    vis1.geometries["box"] = Box(position=[(0, 0, 0)], size=(1, 1, 1))
    assert "box" in vis1.geometries
    vis2 = ZnDraw(url=server, room="testroom", user="u2")
    assert "box" in vis2.geometries

def test_selection_sync(server, s22):
    vis1 = ZnDraw(url=server, room="testroom", user="u1")
    vis1.extend(s22)
    vis1.selection = [0, 1, 2]
    assert vis1.selection == (0, 1, 2)

    vis2 = ZnDraw(url=server, room="testroom", user="u2")
    assert vis2.selection == (0, 1, 2)

def test_frame_selection_sync(server, s22):
    vis1 = ZnDraw(url=server, room="testroom", user="u1")
    vis1.extend(s22)
    vis1.frame_selection = [0, 1, 2]
    assert vis1.frame_selection == (0, 1, 2)

    vis2 = ZnDraw(url=server, room="testroom", user="u2")
    assert vis2.frame_selection == (0, 1, 2)