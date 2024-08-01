from ase.build import molecule

from zndraw import Extension, ZnDraw


def test_run_selection(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)
    vis.selection = [0]
    vis.socket.emit("selection:run", {"method": {"discriminator": "ConnectedParticles"}})
    vis.socket.sleep(7)
    assert vis.selection == [0, 1, 2, 3]


def test_run_modifier(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.selection = [0]
    assert len(vis) == 1
    assert len(vis[-1]) == 3
    vis.socket.emit("modifier:run", {"method": {"discriminator": "Delete"}})
    vis.socket.sleep(7)
    assert len(vis) == 2
    assert len(vis[-1]) == 2


def test_register_modifier(server, s22, water):
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    class AppendWater(Extension):
        def run(self, vis, **kwargs) -> None:
            vis.append(water)
            vis.step = len(vis) - 1

    vis.register_modifier(AppendWater)

    vis.socket.emit("modifier:run", {"method": {"discriminator": "AppendWater"}})
    vis.socket.sleep(7)

    assert vis.atoms == water


def test_locked(server):
    vis = ZnDraw(url=server, token="test_token")

    assert vis.locked is False
    vis.locked = True
    assert vis.locked is True
