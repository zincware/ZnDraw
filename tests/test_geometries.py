import pytest
import znsocket
from ase.build import molecule

from zndraw import Box, Sphere, ZnDraw


def test_geometries(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    assert vis.geometries == []
    vis.geometries = [Box(position=[1, 2, 3]), Sphere()]
    vis.socket.sleep(0.1)

    assert len(vis.geometries) == 2
    assert isinstance(vis.geometries[0], Box)
    assert isinstance(vis.geometries[1], Sphere)

    assert vis.geometries[0].position == [1, 2, 3]
    assert vis.geometries[1].position == [0, 0, 0]

    with pytest.raises(ValueError):
        vis.geometries = ["Geometries are not string!"]

    with pytest.raises(ValueError):
        vis.geometries = "Geometries are not string!"


def test_geometry_selection_position(server):
    vis = ZnDraw(url=server, token="test_token")
    vis.append(molecule("CH4"))
    assert len(vis.geometries) == 0
    geometry_queue = znsocket.Dict(
        r=vis.r,
        socket=vis._refresh_client,
        key=f"queue:{vis.token}:geometry",
    )
    geometry_queue["Plane"] = {
        "material": {
            "color": "#62929e",
            "opacity": 0.2,
            "wireframe": False,
            "outlines": False,
        },
        "width": 10,
        "height": 10,
    }
    vis.socket.emit("room:worker:run")
    vis.socket.sleep(5)

    assert len(vis.geometries) == 1
    assert vis.geometries[0].position == [0, 0, 0]

    # now with a selection
    vis.selection = [1]
    geometry_queue["Plane"] = {
        "material": {
            "color": "#62929e",
            "opacity": 0.2,
            "wireframe": False,
            "outlines": False,
        },
        "width": 10,
        "height": 10,
    }
    vis.socket.emit("room:worker:run")
    vis.socket.sleep(5)

    assert len(vis.geometries) == 2
    assert vis.geometries[1].position == vis.atoms.positions[1].tolist()
    assert vis.geometries[1].position != [0, 0, 0]

    # now with a selection of multiple atoms
    vis.selection = [1, 2]
    geometry_queue["Plane"] = {
        "material": {
            "color": "#62929e",
            "opacity": 0.2,
            "wireframe": False,
            "outlines": False,
        },
        "width": 10,
        "height": 10,
    }
    vis.socket.emit("room:worker:run")
    vis.socket.sleep(5)

    assert len(vis.geometries) == 3
    assert (
        vis.geometries[2].position
        == vis.atoms.get_center_of_mass(indices=[1, 2]).tolist()
    )
