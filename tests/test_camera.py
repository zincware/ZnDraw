from zndraw import ZnDraw


def test_camera(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    assert vis.camera["position"] == [-10, -10, -10]
    assert vis.camera["target"] == [0, 0, 0]
    assert len(vis.camera) == 2

    vis.camera = {"position": [1, 0, 0], "target": [0, 1, 0]}

    assert vis.camera["position"] == [1, 0, 0]
    assert vis.camera["target"] == [0, 1, 0]
