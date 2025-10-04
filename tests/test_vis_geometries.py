import pytest
import requests

from zndraw import ZnDraw
from zndraw.geometries import Sphere


def test_rest_get_geometries(server):
    room = "test-room-geom"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    response = requests.get(f"{server}/api/rooms/{room}/geometries")
    assert response.status_code == 200
    data = response.json()
    assert data == {
        "particles": {
            "type": "Sphere",
            "data": {
                "color": "arrays.colors",
                "position": "arrays.positions",
                "radius": "arrays.radii",
            },
        }
    }


def test_rest_update_geometries(server):
    room = "test-room-geom-update"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    new_geometries = {
        "particles": {
            "type": "Sphere",
            "data": {
                "color": [1.0, 0.0, 0.0],
                "position": [1.0, 1.0, 1.0],
                "radius": 1.0,
            },
        }
    }

    response = requests.put(
        f"{server}/api/rooms/{room}/geometries",
        json={
            "key": "particles",
            "data": new_geometries["particles"]["data"],
            "type": "Sphere",
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    response = requests.get(f"{server}/api/rooms/{room}/geometries")
    assert response.status_code == 200
    data = response.json()
    assert data == new_geometries


def test_rest_add_unknown_geometry(server):
    room = "test-room-geom-unknown"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    response = requests.put(
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
    room = "test-room-geom-delete"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    response = requests.delete(f"{server}/api/rooms/{room}/geometries/particles")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"

    response = requests.get(f"{server}/api/rooms/{room}/geometries")
    assert response.status_code == 200
    data = response.json()
    assert data == {}


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
    vis = ZnDraw(url=server, room="test-room-vis-list-geom", user="tester")
    assert len(vis.geometries) == 1
    assert vis.geometries["particles"] == Sphere()


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
