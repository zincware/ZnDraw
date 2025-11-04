import uuid

import ase
import ase.io
import os
import requests
from fsspec.implementations.local import LocalFileSystem

from conftest import get_jwt_auth_headers
from zndraw import ZnDraw


def test_register_filesystem(server):
    """Test basic filesystem registration with a local filesystem."""

    room = uuid.uuid4().hex
    vis = ZnDraw(room=room, url=server)

    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="test-local")

    # Verify it was registered locally
    assert "test-local" in vis._filesystems
    assert vis._filesystems["test-local"]["fs"] == fs
    assert vis._filesystems["test-local"]["public"] is False

    # use request to list the filesystem via the server API
    response = requests.get(
        f"{server}/api/rooms/{room}/filesystems",
        headers=get_jwt_auth_headers(server, "test-user"),
    )
    assert response.status_code == 200
    data = response.json()
    assert "filesystems" in data
    assert len(data["filesystems"]) == 1
    assert data["filesystems"][0]["name"] == "test-local"
    assert data["filesystems"][0]["fsType"] == "LocalFileSystem"
    assert data["filesystems"][0]["public"] is False

    # Clean up
    vis.disconnect()
    # assert it was unregistered on disconnect
    response = requests.get(
        f"{server}/api/rooms/{room}/filesystems",
        headers=get_jwt_auth_headers(server, "test-user"),
    )
    assert response.status_code == 200
    data = response.json()
    assert "filesystems" in data
    assert len(data["filesystems"]) == 0

def test_register_filesystem_public(server):
    """Test filesystem registration with public=True and cross-room accessibility."""

    room1 = uuid.uuid4().hex
    vis1 = ZnDraw(room=room1, url=server)

    fs = LocalFileSystem()
    vis1.register_filesystem(fs, name="test-local-public", public=True)

    # Verify it was registered locally
    assert "test-local-public" in vis1._filesystems
    assert vis1._filesystems["test-local-public"]["fs"] == fs
    assert vis1._filesystems["test-local-public"]["public"] is True

    # Verify via server API that it was registered as public in room1
    response = requests.get(
        f"{server}/api/rooms/{room1}/filesystems",
        headers=get_jwt_auth_headers(server, "test-user"),
    )
    assert response.status_code == 200
    data = response.json()
    assert "filesystems" in data
    assert len(data["filesystems"]) == 1
    assert data["filesystems"][0]["name"] == "test-local-public"
    assert data["filesystems"][0]["fsType"] == "LocalFileSystem"
    assert data["filesystems"][0]["public"] is True

    # Create a second room and verify the public filesystem is accessible
    room2 = uuid.uuid4().hex
    vis2 = ZnDraw(room=room2, url=server)

    response = requests.get(
        f"{server}/api/rooms/{room2}/filesystems",
        headers=get_jwt_auth_headers(server, "test-user"),
    )
    assert response.status_code == 200
    data = response.json()
    assert "filesystems" in data
    # Should see the public filesystem from room1
    public_fs = next(
        (fs for fs in data["filesystems"] if fs["name"] == "test-local-public"),
        None
    )
    assert public_fs is not None, "Public filesystem should be visible from room2"
    assert public_fs["fsType"] == "LocalFileSystem"
    assert public_fs["public"] is True

    # Clean up
    vis1.disconnect()
    vis2.disconnect()

    # Verify it was unregistered on disconnect from room1
    response = requests.get(
        f"{server}/api/rooms/{room2}/filesystems",
        headers=get_jwt_auth_headers(server, "test-user"),
    )
    assert response.status_code == 200
    data = response.json()
    assert "filesystems" in data
    # Public filesystem should be gone after disconnect
    public_fs = next(
        (fs for fs in data["filesystems"] if fs["name"] == "test-local-public"),
        None
    )
    assert public_fs is None, "Public filesystem should be unregistered after disconnect"


def test_list_with_filesystem(server, tmp_path, s22):
    """Test listing files from a registered filesystem."""
    os.chdir(tmp_path)
    ase.io.write("s22.xyz", s22)

    # register local filesystem
    room = uuid.uuid4().hex
    vis = ZnDraw(room=room, url=server)
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="test-local")
    # list files via server API
    response = requests.get(
        f"{server}/api/rooms/{room}/filesystems/test-local/list",
        headers=get_jwt_auth_headers(server, "test-user"),
        params={"path": "."},
    )
    assert response.status_code == 200
    data = response.json()
    assert "files" in data
    filenames = [f["name"] for f in data["files"]]
    assert "s22.xyz" in filenames

    # Clean up
    vis.disconnect()


def test_load_file_with_filesystem(server, tmp_path, s22):
    """Test loading a file from a registered filesystem into a room."""
    os.chdir(tmp_path)

    # Write test file
    filepath = tmp_path / "s22.xyz"
    ase.io.write(str(filepath), s22)

    # Create room and register filesystem
    room = uuid.uuid4().hex
    vis = ZnDraw(room=room, url=server)
    assert len(vis) == 0  # no frames yet
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="test-local")

    # Load file via server API
    response = requests.post(
        f"{server}/api/rooms/{room}/filesystems/test-local/load",
        headers=get_jwt_auth_headers(server, "test-user"),
        json={"path": str(filepath)},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True
    assert data["frameCount"] == len(s22)

    # Verify frames were actually loaded into the room by checking the frame count
    # ZnDraw client queries the server, so this will see the loaded frames
    assert len(vis) == len(s22), f"Expected {len(s22)} frames, got {len(vis)}"

    frames = list(vis)
    for f_loaded, f_original in zip(frames, s22):
        assert f_loaded == f_original

    # Clean up
    vis.disconnect() 