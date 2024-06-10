import ase.build

from zndraw import ZnDraw


def test_zndraw(server):
    """Test the server fixture."""


def test_set_get(server):
    """Test the server fixture."""
    zndraw = ZnDraw(url=server, token="test_token")
    zndraw[0] = ase.build.molecule("H2O")
    zndraw.socket.sleep(0.5)
    assert len(zndraw) == 1
    assert zndraw[0] == ase.build.molecule("H2O")
