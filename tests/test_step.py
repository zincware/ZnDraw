import pytest

from zndraw import ZnDraw


def test_step(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    assert vis.step == 0
    assert isinstance(vis.step, int)

    assert len(vis) == 22

    with pytest.raises(ValueError):
        vis.step = -1

    with pytest.raises(IndexError):
        vis.step = 22

    with pytest.raises(ValueError):
        vis.step = "string"

    vis.step = 5
    assert vis.step == 5
    assert isinstance(vis.step, int)
