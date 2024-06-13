import pytest

from zndraw import ZnDraw


def test_camera(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    assert vis.bookmarks == {}
    vis.bookmarks = {5: "Hey there!"}
    assert vis.bookmarks == {5: "Hey there!"}

    with pytest.raises(ValueError):
        vis.bookmarks = ["Bookmarks are not a list!"]

    with pytest.raises(ValueError):
        vis.bookmarks = {5: object()}

    with pytest.raises(ValueError):
        vis.bookmarks = {"string": "Hey there!"}
