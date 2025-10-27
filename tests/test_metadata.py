"""Tests for room metadata functionality."""

from datetime import datetime, timezone
from pathlib import Path

import pytest
import requests

import zndraw
from zndraw.app.metadata_manager import RoomMetadataManager


def test_metadata_manager_crud(redis_client):
    """Test basic CRUD operations."""
    manager = RoomMetadataManager(redis_client, "test_room")

    # Test set
    manager.set("key1", "value1")
    assert manager.get("key1") == "value1"

    # Test get non-existent
    assert manager.get("nonexistent") is None

    # Test get_all
    manager.set("key2", "value2")
    all_data = manager.get_all()
    assert all_data == {"key1": "value1", "key2": "value2"}

    # Test update
    manager.update({"key1": "updated", "key3": "value3"})
    assert manager.get("key1") == "updated"
    assert manager.get("key3") == "value3"

    # Test delete
    manager.delete("key2")
    assert manager.get("key2") is None
    assert "key2" not in manager.get_all()

    # Test exists
    assert manager.exists() is True

    # Test clear
    manager.clear()
    assert manager.get_all() == {}
    assert manager.exists() is False


def test_metadata_get_nonexistent_room(redis_client):
    """Test getting metadata for room that doesn't exist."""
    manager = RoomMetadataManager(redis_client, "nonexistent_room")

    assert manager.get("key") is None
    assert manager.get_all() == {}
    assert manager.exists() is False


def test_metadata_update_atomic(redis_client):
    """Test atomic updates."""
    manager = RoomMetadataManager(redis_client, "test_room")

    # Set initial data
    manager.set("key1", "value1")
    manager.set("key2", "value2")

    # Update multiple fields atomically
    manager.update({"key1": "new1", "key2": "new2", "key3": "new3"})

    all_data = manager.get_all()
    assert all_data["key1"] == "new1"
    assert all_data["key2"] == "new2"
    assert all_data["key3"] == "new3"


def test_metadata_clear(redis_client):
    """Test clearing all metadata."""
    manager = RoomMetadataManager(redis_client, "test_room")

    manager.set("key1", "value1")
    manager.set("key2", "value2")
    assert len(manager.get_all()) == 2

    manager.clear()
    assert manager.get_all() == {}
    assert manager.exists() is False


def test_metadata_rest_get(server):
    """Test REST API GET endpoint"""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")
    metadata = {"key1": "value1", "key2": "value2"}
    vis.metadata.update(metadata)

    # Get metadata via REST
    response = requests.get(f"{server}/api/rooms/test_room/metadata")
    assert response.status_code == 200
    data = response.json()
    assert data["metadata"]["key1"] == "value1"
    assert data["metadata"]["key2"] == "value2"


def test_post_room_metadata(server):
    """Test creating/updating metadata."""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")

    # Post metadata
    metadata = {"file_path": "test.xyz", "file_size": "1234"}
    response = requests.post(f"{server}/api/rooms/test_room/metadata", json=metadata)

    assert response.status_code == 200
    data = response.json()
    assert data["success"] is True

    # Verify it was stored
    response = requests.get(f"{server}/api/rooms/test_room/metadata")
    assert response.status_code == 200
    stored = response.json()["metadata"]
    assert stored["file_path"] == "test.xyz"
    assert stored["file_size"] == "1234"


def test_metadata_respects_lock(server, redis_client):
    """Test that metadata writes check permanent room lock (not trajectory lock)."""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")

    # Set permanent room lock
    redis_client.set(f"room:test_room:locked", "1")

    try:
        # Try to update metadata - should be rejected
        response = requests.post(
            f"{server}/api/rooms/test_room/metadata", json={"key": "value"}
        )
        assert response.status_code == 403  # Forbidden due to permanent lock
    finally:
        # Clean up lock
        redis_client.delete(f"room:test_room:locked")

    # Now without lock, should work
    response = requests.post(
        f"{server}/api/rooms/test_room/metadata", json={"key": "value"}
    )
    assert response.status_code == 200


def test_delete_metadata_field(server):
    """Test deleting specific field."""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")

    # Set metadata
    requests.post(
        f"{server}/api/rooms/test_room/metadata",
        json={"key1": "value1", "key2": "value2"},
    )

    # Delete one field
    response = requests.delete(f"{server}/api/rooms/test_room/metadata/key1")
    assert response.status_code == 200

    # Verify key1 is gone but key2 remains
    response = requests.get(f"{server}/api/rooms/test_room/metadata")
    metadata = response.json()["metadata"]
    assert "key1" not in metadata
    assert metadata["key2"] == "value2"


def test_client_metadata_get_set(server):
    """Test vis.metadata dict-like interface."""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")

    # Set via dict interface
    vis.metadata["file"] = "test.xyz"
    vis.metadata["size"] = "1234"

    # Get via dict interface
    assert vis.metadata["file"] == "test.xyz"
    assert vis.metadata["size"] == "1234"

    # Test __len__
    assert len(vis.metadata) == 2

    # Test __iter__
    keys = list(vis.metadata)
    assert "file" in keys
    assert "size" in keys

    # Test __delitem__
    del vis.metadata["size"]
    assert "size" not in vis.metadata
    assert len(vis.metadata) == 1


def test_client_metadata_lazy_load(server):
    """Test lazy loading behavior."""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")

    # Set metadata via API directly
    requests.post(f"{server}/api/rooms/test_room/metadata", json={"key": "value"})

    # Access should trigger lazy load
    assert vis.metadata["key"] == "value"

    # Test refresh
    requests.post(f"{server}/api/rooms/test_room/metadata", json={"key": "updated"})

    vis.metadata.refresh()
    assert vis.metadata["key"] == "updated"


def test_room_list_includes_metadata(server):
    """Test that room list includes metadata."""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")
    vis.metadata["file"] = "test.xyz"

    response = requests.get(f"{server}/api/rooms")
    assert response.status_code == 200

    rooms = response.json()
    test_room = next((r for r in rooms if r["id"] == "test_room"), None)
    assert test_room is not None
    assert "metadata" in test_room
    assert test_room["metadata"]["file"] == "test.xyz"


def test_search_rooms_by_metadata(server):
    """Test searching rooms by metadata values."""
    # Create rooms with different metadata
    vis1 = zndraw.ZnDraw(url=server, room="room1", user="tester")
    vis1.metadata["file"] = "experiment_A.xyz"

    vis2 = zndraw.ZnDraw(url=server, room="room2", user="tester")
    vis2.metadata["file"] = "experiment_B.xyz"

    vis3 = zndraw.ZnDraw(url=server, room="room3", user="tester")
    vis3.metadata["file"] = "control.xyz"

    # Search for "experiment"
    response = requests.get(f"{server}/api/rooms?search=experiment")
    assert response.status_code == 200

    rooms = response.json()
    room_ids = [r["id"] for r in rooms]

    assert "room1" in room_ids
    assert "room2" in room_ids
    assert "room3" not in room_ids


def test_metadata_string_validation(server):
    """Test that non-string values are rejected."""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")

    # Should raise TypeError for non-string
    with pytest.raises(TypeError):
        vis.metadata["number"] = 123

    with pytest.raises(TypeError):
        vis.metadata["list"] = ["a", "b"]


def test_metadata_requires_connection(server):
    """Test that writes require active connection."""
    vis = zndraw.ZnDraw(url=server, room="test_room", user="tester")

    # Disconnect
    vis.disconnect()

    # Should raise RuntimeError
    with pytest.raises(RuntimeError, match="not connected"):
        vis.metadata["key"] = "value"

    with pytest.raises(RuntimeError, match="not connected"):
        del vis.metadata["key"]


def test_file_browser_search(server):
    """Test file browser search functionality."""
    # Test search endpoint with pattern
    response = requests.get(f"{server}/api/files?path=/&search=experiment")
    assert response.status_code == 200
