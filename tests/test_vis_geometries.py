import pytest
import requests

from zndraw import ZnDraw
from zndraw.geometries import Sphere, Bond


def test_rest_get_geometries(server):
    """Test listing geometry keys and getting individual geometries."""
    room = "test-room-geom"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    # Test listing geometry keys
    response = requests.get(f"{server}/api/rooms/{room}/geometries")
    assert response.status_code == 200
    data = response.json()
    assert "geometries" in data
    assert set(data["geometries"]) == {"particles", "bonds", "curve"}
    
    # Test getting individual geometry - particles
    response = requests.get(f"{server}/api/rooms/{room}/geometries/particles")
    assert response.status_code == 200
    data = response.json()
    assert data["key"] == "particles"
    assert data["geometry"]["type"] == "Sphere"
    assert data["geometry"]["data"]["color"] == "arrays.colors"
    assert data["geometry"]["data"]["scale"] == 0.7
    
    # Test getting individual geometry - bonds
    response = requests.get(f"{server}/api/rooms/{room}/geometries/bonds")
    assert response.status_code == 200
    data = response.json()
    assert data["key"] == "bonds"
    assert data["geometry"]["type"] == "Bond"
    assert data["geometry"]["data"]["connectivity"] == "info.connectivity"
    assert data["geometry"]["data"]["scale"] == 0.15
    
    # Test getting non-existent geometry
    response = requests.get(f"{server}/api/rooms/{room}/geometries/nonexistent")
    assert response.status_code == 404


def test_rest_update_geometries(server):
    """Test creating/updating geometries via POST endpoint."""
    room = "test-room-geom-update"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    new_geometry_data = {
        "color": [1.0, 0.0, 0.0],
        "position": [1.0, 1.0, 1.0],
        "radius": 1.0,
    }

    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "particles",
            "data": new_geometry_data,
            "type": "Sphere",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    # Verify the geometry was updated
    response = requests.get(f"{server}/api/rooms/{room}/geometries/particles")
    assert response.status_code == 200
    data = response.json()
    assert data["geometry"]["type"] == "Sphere"
    assert data["geometry"]["data"]["color"] == [1.0, 0.0, 0.0]
    assert data["geometry"]["data"]["position"] == [1.0, 1.0, 1.0]
    assert data["geometry"]["data"]["radius"] == 1.0 


def test_rest_partial_update_geometries(server):
    """Test partially updating a geometry without losing existing data."""
    room = "test-room-geom-partial-update"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    # First, create a geometry with some data
    initial_geometry_data = {
        "position": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        "color": "#FF0000",
    }
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "curve",
            "data": initial_geometry_data,
            "type": "Curve",
        },
    )
    assert response.status_code == 200

    # Now, update only the 'active' status
    partial_update_data = {"active": False}
    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "curve",
            "data": partial_update_data,
            "type": "Curve",
        },
    )
    assert response.status_code == 200

    # Verify the geometry was updated and old data was preserved
    response = requests.get(f"{server}/api/rooms/{room}/geometries/curve")
    assert response.status_code == 200
    data = response.json()
    assert data["geometry"]["type"] == "Curve"
    assert data["geometry"]["data"]["active"] is False
    assert data["geometry"]["data"]["position"] == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
    assert data["geometry"]["data"]["color"] == "#FF0000"


def test_rest_add_unknown_geometry(server):
    """Test that creating geometry with unknown type returns error."""
    room = "test-room-geom-unknown"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    response = requests.post(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "unknown",
            "data": {},
            "type": "UnknownType",
        },
    )
    assert response.status_code == 400
    data = response.json()
    assert data["type"] == "ValueError"
    assert "Unknown geometry type" in data["error"]


def test_rest_delete_geometry(server):
    """Test deleting geometries."""
    room = "test-room-geom-delete"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    response = requests.delete(f"{server}/api/rooms/{room}/geometries/particles")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    response = requests.delete(f"{server}/api/rooms/{room}/geometries/bonds")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    response = requests.delete(f"{server}/api/rooms/{room}/geometries/curve")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    # Verify no geometries remain
    response = requests.get(f"{server}/api/rooms/{room}/geometries")
    assert response.status_code == 200
    data = response.json()
    assert data == {"geometries": []}


def test_rest_delete_unknown_geometry(server):
    room = "test-room-geom-delete-unknown"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    response = requests.delete(f"{server}/api/rooms/{room}/geometries/unknown")
    assert response.status_code == 404
    data = response.json()
    assert data["type"] == "KeyError"
    assert "Geometry with key 'unknown' does not exist" in data["error"]


def test_vis_list_geometries(server):
    from zndraw.geometries import Curve
    
    vis = ZnDraw(url=server, room="test-room-vis-list-geom", user="tester")
    assert len(vis.geometries) == 3
    assert vis.geometries["particles"] == Sphere()
    assert vis.geometries["bonds"] == Bond()
    assert vis.geometries["curve"] == Curve()

def test_vis_add_update_delete_geometry(server):
    vis1 = ZnDraw(url=server, room="room1", user="tester")
    vis2 = ZnDraw(url=server, room="room1", user="tester2")

    new_sphere = Sphere(color=(0.0, 1.0, 0.0), position=(0.0, 0.0, 0.0), radius=2.0)
    vis1.geometries["new_sphere"] = new_sphere
    vis2.socket.sio.sleep(0.5)
    assert vis1.geometries["new_sphere"] == new_sphere
    assert vis2.geometries["new_sphere"] == new_sphere

    del vis2.geometries["new_sphere"]
    vis2.socket.sio.sleep(0.5)
    assert "new_sphere" not in vis1.geometries
    assert "new_sphere" not in vis2.geometries


def test_vis_add_unknown_geometry(server):
    class MyGeometry:
        pass

    vis = ZnDraw(url=server, room="room-unknown-geom", user="tester")
    with pytest.raises(ValueError, match="Unknown geometry type"):
        vis.geometries["my_geom"] = MyGeometry()


def test_vis_delete_unknown_geometry(server):
    vis = ZnDraw(url=server, room="room-del-unknown-geom", user="tester")
    with pytest.raises(KeyError, match="Geometry with key 'unknown' does not exist"):
        del vis.geometries["unknown"]


def test_vis_geometries_key_error(server):
    vis = ZnDraw(url=server, room="room-geom-key-error", user="tester")
    with pytest.raises(
        KeyError, match="Geometry with key 'nonexistent' does not exist"
    ):
        _ = vis.geometries["nonexistent"]