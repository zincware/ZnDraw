"""Tests for filesystem endpoint namespace separation.

Tests verify that:
1. Room endpoints successfully list files from room-scoped filesystems
2. Global endpoints successfully list files from global filesystems
3. Room endpoints return 404 for global filesystems
4. Global endpoints return 404 for room-scoped filesystems
5. Same name in different namespaces works correctly
6. Load endpoints work correctly for both room and global filesystems
"""

import tempfile
from pathlib import Path

import ase.io
import requests
from fsspec.implementations.local import LocalFileSystem
from fsspec.implementations.dirfs import DirFileSystem
from zndraw import ZnDraw


def test_room_endpoint_lists_room_filesystem(server, s22):
    """Room endpoint should successfully list files from room-scoped filesystem."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test files
        test_file1 = Path(tmpdir) / "test1.xyz"
        test_file2 = Path(tmpdir) / "test2.xyz"
        ase.io.write(test_file1, s22[:2])
        ase.io.write(test_file2, s22[:3])

        vis = ZnDraw(url=server, room="test_room", user="admin")
        fs = DirFileSystem(path=tmpdir, fs=LocalFileSystem())
        vis.register_filesystem(fs, name="local", public=False)

        # List files using room endpoint - should succeed
        response = requests.get(
            f"{server}/api/rooms/test_room/filesystems/local/list",
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 200
        data = response.json()
        assert "files" in data
        assert isinstance(data["files"], list)
        # Verify we see our test files
        file_names = [f["name"] for f in data["files"]]
        assert "test1.xyz" in file_names
        assert "test2.xyz" in file_names


def test_global_endpoint_lists_global_filesystem(server, s22):
    """Global endpoint should successfully list files from global filesystem."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test files
        test_file1 = Path(tmpdir) / "global1.xyz"
        test_file2 = Path(tmpdir) / "global2.xyz"
        ase.io.write(test_file1, s22[:2])
        ase.io.write(test_file2, s22[:4])

        vis = ZnDraw(url=server, room="test_room", user="admin")
        fs = DirFileSystem(path=tmpdir, fs=LocalFileSystem())
        vis.register_filesystem(fs, name="local", public=True)

        # List files using global endpoint - should succeed
        response = requests.get(
            f"{server}/api/filesystems/local/list",
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 200
        data = response.json()
        assert "files" in data
        assert isinstance(data["files"], list)
        # Verify we see our test files
        file_names = [f["name"] for f in data["files"]]
        assert "global1.xyz" in file_names
        assert "global2.xyz" in file_names


def test_room_endpoint_rejects_global_filesystem(server):
    """Room endpoint should return 404 for global filesystems."""
    vis = ZnDraw(url=server, room="test_room", user="admin")
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="GlobalFS", public=True)

    # Try to list files using room endpoint - should fail
    response = requests.get(
        f"{server}/api/rooms/test_room/filesystems/GlobalFS/list",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 404
    assert "room-scoped" in response.json()["error"].lower()


def test_global_endpoint_rejects_room_filesystem(server):
    """Global endpoint should return 404 for room-scoped filesystems."""
    vis = ZnDraw(url=server, room="test_room", user="user1")
    fs = LocalFileSystem()
    vis.register_filesystem(fs, name="RoomFS", public=False)

    # Try to list files using global endpoint - should fail
    response = requests.get(
        f"{server}/api/filesystems/RoomFS/list",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 404
    assert "global" in response.json()["error"].lower()


def test_same_name_both_namespaces_list_correctly(server, s22):
    """Filesystems with same name in different namespaces list files correctly."""
    with tempfile.TemporaryDirectory() as tmpdir1, tempfile.TemporaryDirectory() as tmpdir2:
        # Create distinct files in each directory
        room_file = Path(tmpdir1) / "room_file.xyz"
        global_file = Path(tmpdir2) / "global_file.xyz"
        ase.io.write(room_file, s22[:2])
        ase.io.write(global_file, s22[:3])

        vis = ZnDraw(url=server, room="test_room", user="admin")

        # Register room-scoped filesystem
        fs_room = DirFileSystem(path=tmpdir1, fs=LocalFileSystem())
        vis.register_filesystem(fs_room, name="local", public=False)

        # Register global filesystem with same name
        fs_global = DirFileSystem(path=tmpdir2, fs=LocalFileSystem())
        vis.register_filesystem(fs_global, name="local", public=True)

        # List files from room-scoped using room endpoint - should succeed
        response = requests.get(
            f"{server}/api/rooms/test_room/filesystems/local/list",
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 200
        room_data = response.json()
        assert "files" in room_data
        assert isinstance(room_data["files"], list)
        # Verify we see the room file, not the global file
        room_file_names = [f["name"] for f in room_data["files"]]
        assert "room_file.xyz" in room_file_names
        assert "global_file.xyz" not in room_file_names

        # List files from global using global endpoint - should succeed
        response = requests.get(
            f"{server}/api/filesystems/local/list",
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 200
        global_data = response.json()
        assert "files" in global_data
        assert isinstance(global_data["files"], list)
        # Verify we see the global file, not the room file
        global_file_names = [f["name"] for f in global_data["files"]]
        assert "global_file.xyz" in global_file_names
        assert "room_file.xyz" not in global_file_names


def test_room_filesystem_isolated_between_rooms(server):
    """Room-scoped filesystems are isolated between rooms."""
    # Room 1 registers a filesystem
    vis1 = ZnDraw(url=server, room="room1", user="user1")
    fs1 = LocalFileSystem()
    vis1.register_filesystem(fs1, name="local", public=False)

    # List files from room1 - should succeed
    response = requests.get(
        f"{server}/api/rooms/room1/filesystems/local/list",
        headers={"Authorization": f"Bearer {vis1.api.jwt_token}"},
    )
    assert response.status_code == 200

    # Room 2 tries to access room1's filesystem - should fail
    vis2 = ZnDraw(url=server, room="room2", user="user2")
    response = requests.get(
        f"{server}/api/rooms/room2/filesystems/local/list",
        headers={"Authorization": f"Bearer {vis2.api.jwt_token}"},
    )
    assert response.status_code == 404
    assert "room-scoped" in response.json()["error"].lower()


def test_listing_filesystems_returns_both(server):
    """Listing filesystems should return both room and global filesystems."""
    vis = ZnDraw(url=server, room="test_room", user="admin")

    # Register one of each type
    fs_room = LocalFileSystem()
    vis.register_filesystem(fs_room, name="RoomFS", public=False)

    fs_global = LocalFileSystem()
    vis.register_filesystem(fs_global, name="GlobalFS", public=True)

    # Get filesystem list
    response = requests.get(
        f"{server}/api/rooms/test_room/filesystems",
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 200
    filesystems = response.json()

    # Should contain both
    assert len(filesystems) == 2
    names = [fs["name"] for fs in filesystems]
    assert "RoomFS" in names
    assert "GlobalFS" in names

    # Verify public flags
    room_fs = next(fs for fs in filesystems if fs["name"] == "RoomFS")
    global_fs = next(fs for fs in filesystems if fs["name"] == "GlobalFS")
    assert room_fs["public"] is False
    assert global_fs["public"] is True


def test_room_endpoint_loads_room_filesystem(server, s22):
    """Room endpoint should successfully load file from room-scoped filesystem."""
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        ase.io.write(test_file, s22[:5])  # Write first 5 frames

        # Register room-scoped filesystem
        vis = ZnDraw(url=server, room="test_room", user="admin")
        fs = DirFileSystem(path=tmpdir, fs=LocalFileSystem())
        vis.register_filesystem(fs, name="local", public=False)

        # Load file using room endpoint - should succeed (path relative to DirFileSystem root)
        response = requests.post(
            f"{server}/api/rooms/test_room/filesystems/local/load",
            json={"path": "test.xyz", "targetRoom": "test_room"},
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["frameCount"] == 5


def test_global_endpoint_loads_global_filesystem(server, s22):
    """Global endpoint should successfully load file from global filesystem."""
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        ase.io.write(test_file, s22[:3])  # Write first 3 frames

        # Register global filesystem
        vis = ZnDraw(url=server, room="test_room", user="admin")
        fs = DirFileSystem(path=tmpdir, fs=LocalFileSystem())
        vis.register_filesystem(fs, name="local", public=True)

        # Load file using global endpoint - should succeed (path relative to DirFileSystem root)
        response = requests.post(
            f"{server}/api/filesystems/local/load",
            json={"path": "test.xyz", "targetRoom": "test_room"},
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["success"] is True
        assert data["frameCount"] == 3


def test_room_endpoint_rejects_global_filesystem_load(server, s22):
    """Room endpoint should return 404 when trying to load from global filesystem."""
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        ase.io.write(test_file, s22[:2])

        # Register global filesystem
        vis = ZnDraw(url=server, room="test_room", user="admin")
        fs = LocalFileSystem()
        vis.register_filesystem(fs, name="GlobalFS", public=True)

        # Try to load using room endpoint - should fail
        response = requests.post(
            f"{server}/api/rooms/test_room/filesystems/GlobalFS/load",
            json={"path": str(test_file), "targetRoom": "test_room"},
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 404
        assert "room-scoped" in response.json()["error"].lower()


def test_global_endpoint_rejects_room_filesystem_load(server, s22):
    """Global endpoint should return 404 when trying to load from room filesystem."""
    with tempfile.TemporaryDirectory() as tmpdir:
        test_file = Path(tmpdir) / "test.xyz"
        ase.io.write(test_file, s22[:2])

        # Register room-scoped filesystem
        vis = ZnDraw(url=server, room="test_room", user="user1")
        fs = LocalFileSystem()
        vis.register_filesystem(fs, name="RoomFS", public=False)

        # Try to load using global endpoint - should fail
        response = requests.post(
            f"{server}/api/filesystems/RoomFS/load",
            json={"path": str(test_file), "targetRoom": "test_room"},
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 404
        assert "global" in response.json()["error"].lower()


def test_same_name_both_namespaces_load_correctly(server, s22):
    """Load endpoints work correctly with same filesystem name in different namespaces."""
    with tempfile.TemporaryDirectory() as tmpdir1, tempfile.TemporaryDirectory() as tmpdir2:
        room_file = Path(tmpdir1) / "room.xyz"
        global_file = Path(tmpdir2) / "global.xyz"
        ase.io.write(room_file, s22[:4])
        ase.io.write(global_file, s22[:6])

        vis = ZnDraw(url=server, room="test_room", user="admin")

        # Register both room and global with same name, different directories
        fs_room = DirFileSystem(path=tmpdir1, fs=LocalFileSystem())
        vis.register_filesystem(fs_room, name="local", public=False)

        fs_global = DirFileSystem(path=tmpdir2, fs=LocalFileSystem())
        vis.register_filesystem(fs_global, name="local", public=True)

        # Load from room-scoped using room endpoint (path relative to tmpdir1)
        response = requests.post(
            f"{server}/api/rooms/test_room/filesystems/local/load",
            json={"path": "room.xyz", "targetRoom": "test_room"},
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 200
        room_data = response.json()
        assert room_data["success"] is True
        assert room_data["frameCount"] == 4

        # Load from global using global endpoint (path relative to tmpdir2)
        response = requests.post(
            f"{server}/api/filesystems/local/load",
            json={"path": "global.xyz", "targetRoom": "test_room"},
            headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
        )
        assert response.status_code == 200
        global_data = response.json()
        assert global_data["success"] is True
        assert global_data["frameCount"] == 6
