import pytest
import requests

from zndraw.zndraw import ZnDraw


def test_rest_join_new_room(server):
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    assert isinstance(rooms, list)
    assert len(rooms) == 0  # no rooms yet

    room = "test-room-1"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == room
    assert data["frameCount"] == 0
    assert data["selection"] is None
    assert data["frame_selection"] is None
    assert data["created"] is True
    assert data["presenter-lock"] is None
    assert data["step"] == 0
    assert data["bookmarks"] is None
    assert data["template"] == "empty"
    assert data["geometries"] == {
        "particles": {
            "type": "Sphere",
            "data": {
                "color": "arrays.colors",
                "position": "arrays.positions",
                "radius": "arrays.radii",
            },
        }
    }
    assert "settings" in data
    assert "joinToken" in data
    assert "clientId" in data

    # list all rooms again to see if the new room is there
    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200
    rooms = response.json()
    assert isinstance(rooms, list)
    assert len(rooms) == 1
    assert rooms[0]["id"] == room
    assert rooms[0]["template"] == "empty"

    # getting any frame will fail with index error
    for frame_idx in [0, 1, -1, 100]:
        response = requests.get(
            f"{server}/api/rooms/{room}/frames", params={"indices": str(frame_idx)}
        )
        assert response.status_code == 404
        data = response.json()
        assert data["type"] == "IndexError"
        assert (
            data["error"]
            == f"Index out of range for data with 0 frames in room '{room}'"
        )


def test_join_existing_room(server):
    # create a room first
    room = "test-room-1"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200

    # join the existing room with a different user
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["roomId"] == room
    assert data["frameCount"] == 0
    assert data["selection"] is None
    assert data["frame_selection"] is None
    assert data["created"] is False
    assert data["presenter-lock"] is None
    assert data["template"] == "empty"


def test_join_room_with_template(server):
    # create a room and promote to template
    room = "template-room"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 200
    response = requests.post(
        f"{server}/api/rooms/{room}/promote",
        json={"name": "Template Room", "description": "A custom template"},
    )
    assert response.status_code == 200

    # create a new room using the above template
    new_room = "test-room-1"
    response = requests.post(
        f"{server}/api/rooms/{new_room}/join", json={"template": room}
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["template"] == room

    # create a new room without template
    another_room = "test-room-2"
    response = requests.post(f"{server}/api/rooms/{another_room}/join", json={})
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["template"] == "empty"


def test_join_room_invalid_name(server):
    room = "invalid:room"
    response = requests.post(f"{server}/api/rooms/{room}/join", json={})
    assert response.status_code == 400
    data = response.json()
    assert "error" in data
    assert data["error"] == "Room ID cannot contain ':' character"

    with pytest.raises(RuntimeError):
        _ = ZnDraw(url=server, room=room, user="user1")


def test_vis_empty_room(server, s22):
    vis1 = ZnDraw(url=server, room="room-a", user="user1")
    assert len(vis1) == 0
    with pytest.raises(IndexError):
        _ = vis1[0]

    vis1.extend(s22)
    assert len(vis1) == len(s22)
    for i in range(len(s22)):
        assert vis1[i] == s22[i]

    with pytest.raises(IndexError):
        _ = vis1[9999]


def test_vis_len(server, s22):
    vis1 = ZnDraw(url=server, room="room-a", user="user1")
    vis1.extend(s22)
    assert len(vis1) == len(s22)

    vis2 = ZnDraw(url=server, room="room-a", user="user2")
    assert len(vis2) == len(s22)
