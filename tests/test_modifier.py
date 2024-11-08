import numpy as np
import numpy.testing as npt
import znsocket
from ase.build import bulk, molecule

from zndraw import Extension, ZnDraw


def run_queue(vis, key, msg: dict):
    modifier_queue = znsocket.Dict(
        r=vis.r,
        socket=vis._refresh_client,
        key=f"queue:{vis.token}:{key}",
    )
    modifier_queue.update(msg)
    vis.socket.emit("room:worker:run")
    vis.socket.sleep(5)


def test_run_selection(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)
    vis.selection = [0]

    run_queue(vis, "selection", {"ConnectedParticles": {}})

    assert vis.selection == [0, 1, 2, 3]


def test_run_modifier(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.selection = [0]
    assert len(vis) == 1
    assert len(vis[-1]) == 3

    run_queue(vis, "modifier", {"Delete": {}})

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

    run_queue(vis, "modifier", {"AppendWater": {}})

    assert vis.atoms == water


def test_locked(server):
    vis = ZnDraw(url=server, token="test_token")

    assert vis.locked is False
    vis.locked = True
    assert vis.locked is True


##### Tests for each available modifier #####


def test_modify_delete(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.selection = [0]
    assert len(vis) == 1
    assert len(vis[-1]) == 3

    run_queue(vis, "modifier", {"Delete": {}})

    assert len(vis) == 2
    assert len(vis[-1]) == 2


def test_modify_rotate(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.selection = [0, 1]
    vis.points = [[0, 0, 0], [1, 0, 0]]

    run_queue(vis, "modifier", {"Rotate": {"steps": 10}})
    vis.socket.sleep(5)

    assert len(vis) == 11
    # TODO: test that the atoms rotated correctly


def test_modify_translate(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.selection = [0, 1, 2]
    vis.points = [[1, 0, 0], [0, 0, 0]]

    run_queue(vis, "modifier", {"Translate": {"steps": 10}})

    assert len(vis) == 11

    orig_pos = vis[0].positions
    npt.assert_allclose(vis[10].positions, orig_pos - np.array([1, 0, 0]))
    # spline interpolation is not an exact line
    npt.assert_allclose(vis[5].positions, orig_pos - np.array([0.5, 0, 0]), rtol=0.015)

    assert len(vis.points) == 2


def test_modify_duplicate(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.selection = [0]

    run_queue(vis, "modifier", {"Duplicate": {}})

    assert len(vis) == 2
    assert len(vis[0]) == 3
    assert len(vis[1]) == 4  # one duplicated atom


def test_modify_change_type(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.selection = [0]

    run_queue(vis, "modifier", {"ChangeType": {"symbol": "He"}})

    assert vis[1].symbols[0] == "He"


def test_modify_wrap(server):
    vis = ZnDraw(url=server, token="test_token")
    copper = bulk("Cu", cubic=True)
    copper.positions += 5  # shift, so wrapped is recognizable
    vis.extend([copper, copper])

    run_queue(vis, "modifier", {"Wrap": {"all": True}})

    # Wrap is an inplace modifier
    assert len(vis) == 2
    for idx in range(2):
        wrapped_atoms = vis[idx]
        wrapped_atoms.wrap()
        assert not np.allclose(vis[idx].positions, copper.positions)
        assert np.allclose(vis[idx].positions, wrapped_atoms.positions)


def test_modify_replicate(server):
    vis = ZnDraw(url=server, token="test_token")
    wurtzite = bulk("ZnO", "wurtzite", a=3.25, c=5.2)
    vis.extend([wurtzite, wurtzite])

    run_queue(vis, "modifier", {"Replicate": {"x": 2, "y": 2, "z": 2, "all": True}})
    vis.socket.sleep(5)

    # Replicate is an inplace modifier
    assert len(vis) == 2
    for idx in range(2):
        assert len(vis[idx]) == 8 * len(wurtzite)


def test_modify_AddLineParticles(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.points = [[0, 0, 0], [1, 0, 0]]

    run_queue(vis, "modifier", {"AddLineParticles": {"steps": 10, "symbol": "He"}})

    assert len(vis[0]) == 3
    assert len(vis[1]) == 5
    npt.assert_allclose(vis[1].positions[3], [0, 0, 0])
    npt.assert_allclose(vis[1].positions[4], [1, 0, 0])
    assert vis[1].symbols[4] == "He"
    assert vis[1].symbols[3] == "He"


# def test_modify_connect(server):
#     Camera is not available without a webclient
#     vis = ZnDraw(url=server, token="test_token")
#     vis.append(molecule("H2O"))
#     vis.selection = [0, 1]
#     vis.socket.emit("modifier:run", {"method": {"discriminator": "Connect"}})


def test_modify_center(server):
    vis = ZnDraw(url=server, token="test_token")
    copper = bulk("Cu", cubic=True)
    vis.append(copper)
    vis.selection = [0]

    run_queue(vis, "modifier", {"Center": {"all": True}})

    assert np.allclose(vis[0][0].position, np.diag(vis[0].cell) / 2)
    assert not np.allclose(vis[0].positions, copper.positions)


def test_modify_RemoveAtoms(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("H2O"))
    vis.append(molecule("H2O"))
    assert len(vis) == 2
    vis.step = 0

    run_queue(vis, "modifier", {"RemoveAtoms": {}})

    assert len(vis) == 1
