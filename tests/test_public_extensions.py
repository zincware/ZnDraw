"""Tests for public (global) extensions with admin authentication."""

import time
import requests
import pytest
from ase import Atoms
from zndraw import ZnDraw
from zndraw.extensions import Extension, Category


class GlobalTestExtension(Extension):
    """A test extension that will be registered globally."""

    category = Category.MODIFIER

    def run(self, *args, **kwargs):
        """Dummy run method."""
        pass


class RoomModifyingExtension(Extension):
    """Extension that modifies the room's atoms to verify correct room context."""

    category = Category.MODIFIER
    tag: int = 999

    def run(self, vis: ZnDraw, **kwargs):
        """Add a tag to all atoms in the room to mark which room was modified."""
        if len(vis) > 0:
            atoms = vis[0]
            # Add tag to all atoms
            atoms.set_tags([self.tag] * len(atoms))
            vis[0] = atoms


def test_public_extension_visible_across_rooms(server):
    """Test that public extensions registered in one room are visible in other rooms.

    This test verifies the complete flow:
    1. Admin user registers a public extension in room1
    2. Extension is stored with global Redis keys
    3. Extension is visible when querying schemas from room2
    4. Extension is visible when querying schemas from room3

    In LOCAL mode (default for tests), all users are admin.
    """
    # Step 1: Connect to room1 and register a public extension
    vis_room1 = ZnDraw(url=server, room="room1", user="admin_user")
    assert vis_room1.is_admin  # In LOCAL mode, all users are admin

    # Register the extension as public (global)
    vis_room1.register_extension(GlobalTestExtension, public=True)

    # Verify it was registered locally in public namespace
    assert "GlobalTestExtension" in vis_room1._public_extensions
    assert vis_room1._public_extensions["GlobalTestExtension"]["public"] is True

    # Step 2: Query the schema endpoint from room2 (different room)
    # The global extension should be visible here
    response = requests.get(f"{server}/api/rooms/room2/schema/modifiers")
    assert response.status_code == 200
    schemas_room2 = response.json()

    # Schema is now a list of extension objects
    assert isinstance(schemas_room2, list), "Schema should be a list"

    # The global extension should be in the list
    extension_names = [ext["name"] for ext in schemas_room2]
    assert "GlobalTestExtension" in extension_names, (
        f"Global extension not found in room2. Available: {extension_names}"
    )

    # Step 3: Query from room3 (yet another different room)
    response = requests.get(f"{server}/api/rooms/room3/schema/modifiers")
    assert response.status_code == 200
    schemas_room3 = response.json()

    assert isinstance(schemas_room3, list), "Schema should be a list"
    extension_names = [ext["name"] for ext in schemas_room3]
    assert "GlobalTestExtension" in extension_names, (
        f"Global extension not found in room3. Available: {extension_names}"
    )

    # Step 4: Verify the schema is the same across rooms
    # Find the extension object by name
    schema_room2 = next(ext for ext in schemas_room2 if ext["name"] == "GlobalTestExtension")
    schema_room3 = next(ext for ext in schemas_room3 if ext["name"] == "GlobalTestExtension")
    assert schema_room2 == schema_room3, "Schema should be identical across rooms"


def test_global_extension_modifies_correct_room(server):
    """Test that global extensions execute in the correct room context.

    This is the critical test that verifies the fix:
    1. Register a global extension in room1 (worker stays in room1)
    2. Trigger the extension from room2
    3. Verify that room2's data is modified, NOT room1's data

    This test catches the bug where extensions would modify the wrong room.
    """
    # Setup: Create atoms in both rooms
    atoms_room1 = Atoms("H2O", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])
    atoms_room2 = Atoms("CH4", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0)])

    # Connect to room1 and add atoms
    vis_room1 = ZnDraw(url=server, room="room1", user="worker_user", auto_pickup_jobs=True)
    vis_room1.extend([atoms_room1])
    assert len(vis_room1) == 1
    assert len(vis_room1[0]) == 3  # H2O has 3 atoms

    # Connect to room2 and add atoms
    vis_room2 = ZnDraw(url=server, room="room2", user="user2", auto_pickup_jobs=False)
    vis_room2.extend([atoms_room2])
    assert len(vis_room2) == 1
    assert len(vis_room2[0]) == 5  # CH4 has 5 atoms

    # Register a global extension in room1 with a specific tag value
    vis_room1.register_extension(RoomModifyingExtension, public=True)
    time.sleep(0.5)  # Give registration time to complete

    # Trigger the extension from room2 with tag=222
    response = requests.post(
        f"{server}/api/rooms/room2/extensions/public/modifiers/RoomModifyingExtension/submit",
        json={"data": {"tag": 222}},
        headers={"Authorization": f"Bearer {vis_room2.api.jwt_token}"},
    )
    assert response.status_code in [200, 201], f"Failed to submit job: {response.text}"

    # Wait for job to complete
    time.sleep(2)

    # Verify room2 was modified (should have tag 222)
    vis_room2_check = ZnDraw(url=server, room="room2", user="checker2")
    atoms_room2_after = vis_room2_check[0]
    assert all(tag == 222 for tag in atoms_room2_after.get_tags()), (
        f"Room2 atoms should have tag 222, but got {atoms_room2_after.get_tags()}. "
        "This means the extension executed in the wrong room!"
    )

    # Verify room1 was NOT modified (should still have no tags or default tags)
    vis_room1_check = ZnDraw(url=server, room="room1", user="checker1")
    atoms_room1_after = vis_room1_check[0]
    room1_tags = atoms_room1_after.get_tags()
    assert all(tag != 222 for tag in room1_tags), (
        f"Room1 atoms should NOT have tag 222, but got {room1_tags}. "
        "This means the extension incorrectly modified room1 instead of room2!"
    )


def test_global_extension_cleanup_on_disconnect(server):
    """Test that global extensions are removed from all rooms when worker disconnects.

    This verifies the cleanup behavior:
    1. Register a global extension in room1
    2. Verify extension is visible in room2 and room3
    3. Disconnect the worker (close the ZnDraw instance)
    4. Verify extension is no longer visible in any room
    """
    # Step 1: Connect to room1 and register a global extension
    vis_room1 = ZnDraw(url=server, room="room1", user="worker_disconnect")
    vis_room1.register_extension(GlobalTestExtension, public=True)
    time.sleep(0.5)  # Give registration time to complete

    # Step 2: Verify extension is visible in multiple rooms
    response_room2 = requests.get(f"{server}/api/rooms/room2/schema/modifiers")
    assert response_room2.status_code == 200
    schemas_room2_before = response_room2.json()

    # Schema is now a list of extension objects
    assert isinstance(schemas_room2_before, list), "Schema should be a list"
    extension_names_room2 = [ext["name"] for ext in schemas_room2_before]
    assert "GlobalTestExtension" in extension_names_room2, (
        "Global extension should be visible in room2 before disconnect"
    )

    response_room3 = requests.get(f"{server}/api/rooms/room3/schema/modifiers")
    assert response_room3.status_code == 200
    schemas_room3_before = response_room3.json()

    assert isinstance(schemas_room3_before, list), "Schema should be a list"
    extension_names_room3 = [ext["name"] for ext in schemas_room3_before]
    assert "GlobalTestExtension" in extension_names_room3, (
        "Global extension should be visible in room3 before disconnect"
    )

    # Step 3: Disconnect the worker
    vis_room1.socket.disconnect()
    vis_room1.socket.sio.disconnect()  # Ensure full disconnect
    time.sleep(1)  # Give server time to process disconnect

    # Step 4: Verify extension is no longer visible in any room
    response_room2_after = requests.get(f"{server}/api/rooms/room2/schema/modifiers")
    assert response_room2_after.status_code == 200
    schemas_room2_after = response_room2_after.json()

    assert isinstance(schemas_room2_after, list), "Schema should be a list"
    extension_names_room2_after = [ext["name"] for ext in schemas_room2_after]
    assert "GlobalTestExtension" not in extension_names_room2_after, (
        f"Global extension should be removed from room2 after disconnect. "
        f"Available extensions: {extension_names_room2_after}"
    )

    response_room3_after = requests.get(f"{server}/api/rooms/room3/schema/modifiers")
    assert response_room3_after.status_code == 200
    schemas_room3_after = response_room3_after.json()

    assert isinstance(schemas_room3_after, list), "Schema should be a list"
    extension_names_room3_after = [ext["name"] for ext in schemas_room3_after]
    assert "GlobalTestExtension" not in extension_names_room3_after, (
        f"Global extension should be removed from room3 after disconnect. "
        f"Available extensions: {extension_names_room3_after}"
    )

    # Also verify it's removed from room1
    response_room1_after = requests.get(f"{server}/api/rooms/room1/schema/modifiers")
    assert response_room1_after.status_code == 200
    schemas_room1_after = response_room1_after.json()

    assert isinstance(schemas_room1_after, list), "Schema should be a list"
    extension_names_room1_after = [ext["name"] for ext in schemas_room1_after]
    assert "GlobalTestExtension" not in extension_names_room1_after, (
        f"Global extension should be removed from room1 after disconnect. "
        f"Available extensions: {extension_names_room1_after}"
    )
