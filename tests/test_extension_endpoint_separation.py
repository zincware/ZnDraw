"""Tests for extension endpoint namespace separation.

Tests verify that:
1. Private endpoint successfully submits to room-scoped extensions
2. Public endpoint successfully submits to global extensions
3. Private endpoint returns 404 for global extensions
4. Public endpoint returns 404 for room-scoped extensions
5. Celery extensions only work on public endpoint
6. Same name in different namespaces works correctly
"""

import requests
from zndraw import ZnDraw
from zndraw.extensions import Extension, Category


class TestExtension(Extension):
    """Test extension for testing."""

    category = Category.MODIFIER
    value: int = 1

    def run(self, *args, **kwargs):
        """Dummy run method."""
        pass


class GlobalExt(Extension):
    """Global test extension."""

    category = Category.MODIFIER
    value: int = 1

    def run(self, *args, **kwargs):
        """Dummy run method."""
        pass


class RoomExt(Extension):
    """Room test extension."""

    category = Category.MODIFIER
    value: int = 1

    def run(self, *args, **kwargs):
        """Dummy run method."""
        pass


def test_private_endpoint_submits_to_room_extension(server):
    """Private endpoint should successfully submit to room-scoped extension."""
    vis = ZnDraw(url=server, room="test_room", user="user1")
    vis.register_extension(TestExtension, public=False)

    # Submit using private endpoint - should succeed
    response = requests.post(
        f"{server}/api/rooms/test_room/extensions/private/modifiers/TestExtension/submit",
        json={"data": {"value": 42}},
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"
    assert "jobId" in data


def test_public_endpoint_submits_to_global_extension(server):
    """Public endpoint should successfully submit to global extension."""
    vis = ZnDraw(url=server, room="test_room", user="admin")
    vis.register_extension(TestExtension, public=True)

    # Submit using public endpoint - should succeed
    response = requests.post(
        f"{server}/api/rooms/test_room/extensions/public/modifiers/TestExtension/submit",
        json={"data": {"value": 99}},
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"
    assert "jobId" in data


def test_private_endpoint_rejects_global_extension(server):
    """Private endpoint should return 404 for global extensions."""
    vis = ZnDraw(url=server, room="test_room", user="admin")
    vis.register_extension(GlobalExt, public=True)

    # Try to submit using private endpoint - should fail
    response = requests.post(
        f"{server}/api/rooms/test_room/extensions/private/modifiers/GlobalExt/submit",
        json={"data": {"value": 1}},
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 404
    response_data = response.json()
    assert "room-scoped" in response_data["error"].lower()
    assert response_data["code"] == "EXTENSION_NOT_FOUND"
    assert response_data["namespace"] == "room-scoped"


def test_public_endpoint_rejects_room_extension(server):
    """Public endpoint should return 404 for room-scoped extensions."""
    vis = ZnDraw(url=server, room="test_room", user="user1")
    vis.register_extension(RoomExt, public=False)

    # Try to submit using public endpoint - should fail
    response = requests.post(
        f"{server}/api/rooms/test_room/extensions/public/modifiers/RoomExt/submit",
        json={"data": {"value": 1}},
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 404
    response_data = response.json()
    assert "global" in response_data["error"].lower()
    assert response_data["code"] == "EXTENSION_NOT_FOUND"
    assert response_data["namespace"] == "global"


def test_celery_extension_requires_public_endpoint(server):
    """Server-side celery extensions must use public endpoint."""
    vis = ZnDraw(url=server, room="test_room", user="admin")

    # Try to submit celery extension using private endpoint - should fail
    # "Delete" is a server-side celery extension in modifiers
    response = requests.post(
        f"{server}/api/rooms/test_room/extensions/private/modifiers/Delete/submit",
        json={"data": {}},
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 400
    response_data = response.json()
    assert "public endpoint" in response_data["error"].lower()
    assert response_data["code"] == "WRONG_ENDPOINT"
    assert response_data["expectedEndpoint"] == "public"


def test_celery_extension_works_on_public_endpoint(server):
    """Server-side celery extensions work on public endpoint."""
    vis = ZnDraw(url=server, room="test_room", user="admin")

    # Submit celery extension using public endpoint - should succeed
    response = requests.post(
        f"{server}/api/rooms/test_room/extensions/public/modifiers/Delete/submit",
        json={"data": {}},
        headers={"Authorization": f"Bearer {vis.api.jwt_token}"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "success"


def test_room_extension_isolated_between_rooms(server):
    """Room-scoped extensions are isolated between rooms."""
    # Room 1 registers an extension
    vis1 = ZnDraw(url=server, room="room1", user="user1")
    vis1.register_extension(RoomExt, public=False)

    # Submit to room1 - should succeed
    response = requests.post(
        f"{server}/api/rooms/room1/extensions/private/modifiers/RoomExt/submit",
        json={"data": {"value": 1}},
        headers={"Authorization": f"Bearer {vis1.api.jwt_token}"},
    )
    assert response.status_code == 200

    # Room 2 tries to access room1's extension - should fail
    vis2 = ZnDraw(url=server, room="room2", user="user2")
    response = requests.post(
        f"{server}/api/rooms/room2/extensions/private/modifiers/RoomExt/submit",
        json={"data": {"value": 1}},
        headers={"Authorization": f"Bearer {vis2.api.jwt_token}"},
    )
    assert response.status_code == 404
    response_data = response.json()
    assert "room-scoped" in response_data["error"].lower()
    assert response_data["code"] == "EXTENSION_NOT_FOUND"


def test_public_extension_accessible_from_multiple_rooms(server):
    """Global extensions are accessible from any room."""
    # Register global extension from room1
    vis1 = ZnDraw(url=server, room="room1", user="admin")
    vis1.register_extension(GlobalExt, public=True)

    # Submit from room1 - should succeed
    response = requests.post(
        f"{server}/api/rooms/room1/extensions/public/modifiers/GlobalExt/submit",
        json={"data": {"value": 1}},
        headers={"Authorization": f"Bearer {vis1.api.jwt_token}"},
    )
    assert response.status_code == 200

    # Submit from room2 - should also succeed (global extension)
    vis2 = ZnDraw(url=server, room="room2", user="user2")
    response = requests.post(
        f"{server}/api/rooms/room2/extensions/public/modifiers/GlobalExt/submit",
        json={"data": {"value": 2}},
        headers={"Authorization": f"Bearer {vis2.api.jwt_token}"},
    )
    assert response.status_code == 200


def test_schema_endpoint_includes_public_flag(server):
    """Schema endpoint returns public flag for each extension."""
    vis = ZnDraw(url=server, room="test_room", user="admin")

    # Register one of each type
    vis.register_extension(RoomExt, public=False)
    vis.register_extension(GlobalExt, public=True)

    # Get schemas
    response = requests.get(
        f"{server}/api/rooms/test_room/schema/modifiers",
    )
    assert response.status_code == 200
    schemas = response.json()

    # Schemas is now a list of extension objects
    assert isinstance(schemas, list), "Schemas should be a list"

    # Helper to find extension by name
    def find_extension(name):
        return next((ext for ext in schemas if ext.get("name") == name), None)

    # Check room-scoped extension has public=False
    room_ext = find_extension("RoomExt")
    assert room_ext is not None, "RoomExt should be in schemas"
    assert room_ext["public"] is False

    # Check global extension has public=True
    global_ext = find_extension("GlobalExt")
    assert global_ext is not None, "GlobalExt should be in schemas"
    assert global_ext["public"] is True

    # Check server-side (celery) extensions have public=True
    # "Delete" is a server-side extension
    delete_ext = find_extension("Delete")
    assert delete_ext is not None, "Delete should be in schemas"
    assert delete_ext["public"] is True
    assert delete_ext["provider"] == "celery"
