import pytest

from zndraw import Extension, ZnDraw


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


def test_run_selection(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)
    vis.selection = [0]

    vis.socket.emit("selection:run", {"method": {"discriminator": "ConnectedParticles"}})
    vis.socket.sleep(5)
    assert vis.selection == [0, 1, 2, 3]


def test_register_modifier(server, s22, water):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    class AppendWater(Extension):
        def run(self, vis, **kwargs) -> None:
            vis.append(water)
            vis.step = len(vis) - 1

    vis.register_modifier(AppendWater)

    vis.socket.emit("modifier:run", {"method": {"discriminator": "AppendWater"}})
    vis.socket.sleep(5)

    assert vis.atoms == water
