from zndraw import ZnDraw


def test_zndraw(server):
    """Test the server fixture."""


def test_set_get(server, s22):
    """Test the server fixture."""
    zndraw = ZnDraw(url=server, token="test_token")
    zndraw[0] = s22[0]
    zndraw.socket.sleep(0.5)
    assert len(zndraw) == 1
    assert zndraw[0] == s22[0]


def test_extend(server, s22):
    """Test the server fixture."""
    zndraw = ZnDraw(url=server, token="test_token")
    zndraw.extend(s22)
    zndraw.socket.sleep(0.5)
    assert len(zndraw) == len(s22)
    # TODO: when using extend, it should really be 23, should it not?
