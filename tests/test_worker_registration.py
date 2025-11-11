"""Tests for worker registration via REST endpoint.

Tests cover:
1. Basic registration with valid schema
2. Re-registration with same schema (allowed)
3. Schema mismatch error (different schema for same extension)
4. Room vs global extension conflict (allowed)
5. Global extension duplication with different schema (error)
6. Admin privilege enforcement
7. Reserved name protection
"""

import time

import pytest
import requests
from ase import Atoms
from zndraw import ZnDraw
from zndraw.extensions import Category, Extension


class TestExtension(Extension):
    """A test extension for registration tests."""

    category = Category.MODIFIER

    def run(self, *args, **kwargs):
        """Dummy run method."""
        pass

class RemoveAtoms(Extension):
    """Extension with different schema for testing."""

    category = Category.MODIFIER
    parameter: str = "does-not-exist-on-global-extension"

    def run(self, atoms: Atoms, indices: list[int], **kwargs):
        ...

# We can't easily create two extension classes with the same __name__ in Python
# So we'll dynamically create a class with a different schema for testing
def create_test_extension_with_different_schema():
    """Dynamically create TestExtension class with different schema."""

    class TestExtension(Extension):
        """A test extension with different schema (has param field)."""
        category = Category.MODIFIER
        param: int = 42  # This field makes schema different

        def run(self, *args, **kwargs):
            """Dummy run method."""
            pass

    return TestExtension


def test_register_room_extension_basic(server):
    """Test basic room-scoped extension registration."""
    vis = ZnDraw(url=server, room="test_room", user="test_user")

    # Register extension
    vis.register_extension(TestExtension, public=False)

    # Verify extension appears in schema
    response = requests.get(f"{server}/api/rooms/test_room/schema/modifiers")
    assert response.status_code == 200
    schemas = response.json()
    assert isinstance(schemas, list), "Schema should be a list"
    extension_names = [ext["name"] for ext in schemas]
    assert "TestExtension" in extension_names


def test_register_same_schema_twice_allowed(server):
    """Test that registering an extension with the same schema twice is allowed."""
    vis1 = ZnDraw(url=server, room="test_room", user="worker1")
    vis2 = ZnDraw(url=server, room="test_room", user="worker2")

    # Both workers register the same extension with same schema
    vis1.register_extension(TestExtension, public=False)
    vis2.register_extension(TestExtension, public=False)

    # Verify extension is registered (should succeed)
    response = requests.get(f"{server}/api/rooms/test_room/schema/modifiers")
    assert response.status_code == 200
    schemas = response.json()
    assert isinstance(schemas, list), "Schema should be a list"
    extension_names = [ext["name"] for ext in schemas]
    assert "TestExtension" in extension_names


def test_register_different_schema_error(server):
    """Test that registering an extension with different schema fails."""
    vis1 = ZnDraw(url=server, room="test_room", user="worker1")
    vis2 = ZnDraw(url=server, room="test_room", user="worker2")

    # First worker registers TestExtension
    vis1.register_extension(TestExtension, public=False)
    DifferentTestExtension = create_test_extension_with_different_schema()
    assert DifferentTestExtension.__name__ == TestExtension.__name__

    # This should fail with RuntimeError
    with pytest.raises(RuntimeError, match="different schema"):
        vis2.register_extension(DifferentTestExtension, public=False)


def test_room_extension_overrides_global_allowed(server):
    """Test that room-scoped extension can override global extension (allowed)."""
    vis_global = ZnDraw(url=server, room="global_room", user="admin_user")
    vis_room = ZnDraw(url=server, room="specific_room", user="room_user")

    # Register as global extension first
    vis_global.register_extension(TestExtension, public=True)
    # Register same extension in specific room (should be allowed)
    vis_room.register_extension(TestExtension, public=False)

    # Both should exist
    global_response = requests.get(f"{server}/api/rooms/global_room/schema/modifiers")
    global_schemas = global_response.json()
    assert isinstance(global_schemas, list), "Schema should be a list"
    global_names = [ext["name"] for ext in global_schemas]
    assert "TestExtension" in global_names

    room_response = requests.get(f"{server}/api/rooms/specific_room/schema/modifiers")
    room_schemas = room_response.json()
    assert isinstance(room_schemas, list), "Schema should be a list"
    room_names = [ext["name"] for ext in room_schemas]
    assert "TestExtension" in room_names


def test_global_extension_duplicate_schema_error(server):
    """Test that registering a global extension twice with different schema fails."""
    vis1 = ZnDraw(url=server, room="room1", user="admin1")
    vis2 = ZnDraw(url=server, room="room2", user="admin2")

    # First admin registers global extension
    vis1.register_extension(TestExtension, public=True)
    DifferentTestExtension = create_test_extension_with_different_schema()

    assert DifferentTestExtension.__name__ == TestExtension.__name__

    # Second admin tries to register global extension with different schema but same name
    with pytest.raises(RuntimeError, match="different schema"):
        vis2.register_extension(DifferentTestExtension, public=True)


def test_reserved_name_protection_per_room(server):
    """Test that extensions cannot override server-side extension names."""
    vis = ZnDraw(url=server, room="test_room", user="test_user")
    vis.register_extension(RemoveAtoms, public=False) # non-global is allowed

    response = requests.get(f"{server}/api/rooms/test_room/schema/modifiers")
    assert response.status_code == 200
    schemas = response.json()
    assert isinstance(schemas, list), "Schema should be a list"
    extension_names = [ext["name"] for ext in schemas]
    assert "RemoveAtoms" in extension_names

def test_reserved_name_protection_global(server):
    """Test that extensions cannot override server-side extension names globally."""
    vis = ZnDraw(url=server, room="any_room", user="admin_user")

    # Attempt to register RemoveAtoms as global extension should fail, because 
    # it is a server side extension
    with pytest.raises(RuntimeError, match="reserved"):
        vis.register_extension(RemoveAtoms, public=True)


def test_session_id_required(server):
    """Test that sessionId is required in the request."""
    vis = ZnDraw(url=server, room="test_room", user="test_user")

    # Make a direct REST call without sessionId
    response = requests.post(
        f"{server}/api/workers/register",
        json={
            # Missing sessionId
            "roomId": "test_room",
            "name": "TestExtension",
            "category": "modifiers",
            "schema": {},
            "public": False,
        },
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )

    assert response.status_code == 400
    data = response.json()
    assert data["success"] is False
    assert "sessionId" in data["error"].lower() or data.get("code") == "MISSING_FIELD"


def test_schema_validation_on_registration(server):
    """Test that schema is properly validated and stored."""
    vis = ZnDraw(url=server, room="test_room", user="test_user")

    # Register extension
    vis.register_extension(TestExtension, public=False)
    time.sleep(0.5)

    # Verify schema was stored correctly
    response = requests.get(f"{server}/api/rooms/test_room/schema/modifiers")
    assert response.status_code == 200
    schemas = response.json()

    assert isinstance(schemas, list), "Schema should be a list"
    extension_names = [ext["name"] for ext in schemas]
    assert "TestExtension" in extension_names

    # Find the TestExtension in the list
    test_ext = next(ext for ext in schemas if ext["name"] == "TestExtension")
    schema = test_ext["schema"]
    # Schema should be a dict with type definitions
    assert isinstance(schema, dict)


def test_disconnect_cleans_up_extension(server):
    """Test that disconnecting removes the extension if no other workers."""
    vis = ZnDraw(url=server, room="test_room", user="worker_cleanup")

    # Register extension
    vis.register_extension(TestExtension, public=False)
    # Verify extension exists
    response = requests.get(f"{server}/api/rooms/test_room/schema/modifiers")
    schemas = response.json()
    assert isinstance(schemas, list), "Schema should be a list"
    extension_names = [ext["name"] for ext in schemas]
    assert "TestExtension" in extension_names

    # Disconnect
    vis.socket.disconnect()
    vis.socket.sio.disconnect()
    time.sleep(1)

    # Extension should be removed (assuming no queued jobs)
    response = requests.get(f"{server}/api/rooms/test_room/schema/modifiers")
    schemas_after = response.json()
    assert isinstance(schemas_after, list), "Schema should be a list"
    extension_names_after = [ext["name"] for ext in schemas_after]
    assert "TestExtension" not in extension_names_after
