from zndraw import ZnDraw
import requests
from fsspec.implementations.local import LocalFileSystem
import pytest
import time


def test_register_room_fs_basic(server):
    """Test registering a per-room filesystem."""
    vis = ZnDraw(url=server, room="test_room", user="user1")
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="LocalFS", public=False)
    response = requests.get(
        f"{server}/api/rooms/test_room/filesystems",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 200
    filesystems = response.json()
    assert filesystems == [
        {
            "name": "LocalFS",
            "fsType": "LocalFileSystem",
            "public": False,
            "sessionId": vis.api.session_id,
        }
    ]


def test_register_global_fs_basic(server):
    """Test registering a global filesystem."""
    vis = ZnDraw(url=server, room="global_room", user="admin_user")
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="GlobalLocalFS", public=True)
    response = requests.get(
        f"{server}/api/rooms/global_room/filesystems",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 200
    filesystems = response.json()
    assert filesystems == [
        {
            "name": "GlobalLocalFS",
            "fsType": "LocalFileSystem",
            "public": True,
            "sessionId": vis.api.session_id,
        }
    ]


def test_register_room_fs_twice(server):
    """Test that registering the same room filesystem twice fails."""
    vis = ZnDraw(url=server, room="test_room", user="user1")
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="DuplicateFS", public=False)

    with pytest.raises(ValueError, match="already registered"):
        vis.register_filesystem(fs, name="DuplicateFS", public=False)


def test_register_global_fs_twice(server):
    """Test that registering the same global filesystem twice fails."""
    vis = ZnDraw(url=server, room="global_room", user="admin_user")
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="GlobalDuplicateFS", public=True)

    with pytest.raises(ValueError, match="already registered"):
        vis.register_filesystem(fs, name="GlobalDuplicateFS", public=True)


def test_room_and_global_fs_same_name(server):
    """Test that room and global filesystems can have the same name."""
    vis = ZnDraw(url=server, room="test-room", user="admin_user")

    fs_global = LocalFileSystem()
    vis.register_filesystem(fs_global, name="SharedFS", public=True)

    fs_room = LocalFileSystem()
    vis.register_filesystem(fs_room, name="SharedFS", public=False)

    response_global = requests.get(
        f"{server}/api/rooms/test-room/filesystems",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    response_json = response_global.json()
    assert len(response_json) == 2
    names = [fs["name"] for fs in response_json]
    assert names.count("SharedFS") == 2
    publics = [fs["public"] for fs in response_json if fs["name"] == "SharedFS"]
    assert sorted(publics) == [False, True]

@pytest.mark.parametrize(
    "register_global",
    [True, False],
)
def test_register_fs_disconnect(server, register_global):
    """Test that registering a filesystem that is disconnected fails."""
    vis = ZnDraw(url=server, room="test_room", user="user1")
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="ConnectedFS", public=register_global)
    response = requests.get(
        f"{server}/api/rooms/test_room/filesystems",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    filesystems = response.json()
    assert any(f["name"] == "ConnectedFS" for f in filesystems)

    vis.disconnect()
    time.sleep(0.1)

    response = requests.get(
        f"{server}/api/rooms/test_room/filesystems",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    filesystems = response.json()
    assert all(f["name"] != "ConnectedFS" for f in filesystems)
