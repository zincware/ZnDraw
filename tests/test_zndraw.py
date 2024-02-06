from zndraw import ZnDraw


def test_zndraw(server):
    vis = ZnDraw(server)
    assert vis.socket.call("ping") == "pong"