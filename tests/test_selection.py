import pytest

from zndraw import ZnDraw


def test_selection(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    assert len(vis.atoms) == 8

    with pytest.raises(ValueError):
        vis.selection = ["Hello"]

    with pytest.raises(IndexError):
        vis.selection = [0, 1, 2, 25]

    with pytest.raises(IndexError):
        vis.selection = [0, 1, 2, -10]

    vis.selection = [0, 7, 6, 5, 4]
    assert vis.selection == [0, 7, 6, 5, 4]
