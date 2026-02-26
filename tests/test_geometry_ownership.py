"""Integration tests for geometry ownership and admin lock enforcement.

Tests use the real server (server_auth fixture) with admin/guest role separation
and the Pydantic model API: vis.geometries["key"] = Sphere(owner=...).
"""

import uuid

import jwt
import pytest

from zndraw.client import RoomLockedError, ZnDraw
from zndraw.geometries import Sphere


def _get_user_id(vis: ZnDraw) -> str:
    """Extract user ID from the client's JWT token."""
    assert vis.api.token is not None
    payload = jwt.decode(vis.api.token, options={"verify_signature": False})
    return payload["sub"]


def _lock_room(vis: ZnDraw) -> None:
    """Set admin lock on the room via PATCH."""
    response = vis.api.http.patch(
        f"/v1/rooms/{vis.room}",
        json={"locked": True},
        headers=vis.api._headers(),
    )
    vis.api.raise_for_status(response)


# =============================================================================
# Basic Ownership Tests
# =============================================================================


def test_owner_claims_own_geometry(server_auth: str) -> None:
    """Owner can claim a geometry by setting owner to their own ID."""
    vis = ZnDraw(url=server_auth)
    user_id = _get_user_id(vis)

    vis.geometries["sphere"] = Sphere(owner=user_id)
    retrieved = vis.geometries["sphere"]
    assert retrieved.owner == user_id

    vis.disconnect()


def test_owner_releases_geometry(server_auth: str) -> None:
    """Owner can release a geometry by setting owner to None."""
    vis = ZnDraw(url=server_auth)
    user_id = _get_user_id(vis)

    vis.geometries["sphere"] = Sphere(owner=user_id)
    vis.geometries["sphere"] = Sphere(owner=None)

    retrieved = vis.geometries["sphere"]
    assert retrieved.owner is None

    vis.disconnect()


# =============================================================================
# Non-Owner Rejection Tests
# =============================================================================


def test_non_owner_cannot_edit_owned_geometry(server_auth: str) -> None:
    """Non-owner gets PermissionError when editing an owned geometry."""
    room_id = uuid.uuid4().hex
    vis_a = ZnDraw(url=server_auth, room=room_id)
    vis_b = ZnDraw(url=server_auth, room=room_id)
    user_a_id = _get_user_id(vis_a)

    vis_a.geometries["sphere"] = Sphere(owner=user_a_id)

    with pytest.raises(PermissionError):
        vis_b.geometries["sphere"] = Sphere(radius=[99.0])

    vis_a.disconnect()
    vis_b.disconnect()


def test_non_owner_cannot_delete_owned_geometry(server_auth: str) -> None:
    """Non-owner gets PermissionError when deleting an owned geometry."""
    room_id = uuid.uuid4().hex
    vis_a = ZnDraw(url=server_auth, room=room_id)
    vis_b = ZnDraw(url=server_auth, room=room_id)
    user_a_id = _get_user_id(vis_a)

    vis_a.geometries["sphere"] = Sphere(owner=user_a_id)

    with pytest.raises(PermissionError):
        del vis_b.geometries["sphere"]

    vis_a.disconnect()
    vis_b.disconnect()


# =============================================================================
# Selection Ownership Tests
# =============================================================================


def test_selection_on_owned_geometry_owner_succeeds(server_auth: str) -> None:
    """Owner can update selection on their own geometry."""
    vis = ZnDraw(url=server_auth)
    user_id = _get_user_id(vis)

    vis.geometries["sphere"] = Sphere(owner=user_id)
    vis.selections["sphere"] = (0, 1)

    assert vis.selections["sphere"] == (0, 1)

    vis.disconnect()


def test_selection_on_owned_geometry_non_owner_blocked(server_auth: str) -> None:
    """Non-owner gets PermissionError when updating selection on owned geometry."""
    room_id = uuid.uuid4().hex
    vis_a = ZnDraw(url=server_auth, room=room_id)
    vis_b = ZnDraw(url=server_auth, room=room_id)
    user_a_id = _get_user_id(vis_a)

    vis_a.geometries["sphere"] = Sphere(owner=user_a_id)

    with pytest.raises(PermissionError):
        vis_b.selections["sphere"] = (0, 1)

    vis_a.disconnect()
    vis_b.disconnect()


# =============================================================================
# Admin Lock + Unowned Geometry Tests
# =============================================================================


def test_admin_locked_unowned_non_superuser_blocked(server_auth: str) -> None:
    """Non-superuser cannot edit unowned geometries in admin-locked room."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    guest = ZnDraw(url=server_auth, room=admin.room)

    # Guest creates geometry first (room unlocked)
    guest.geometries["sphere"] = Sphere()

    # Admin locks room
    _lock_room(admin)

    # Guest tries to edit unowned geometry → blocked
    with pytest.raises(RoomLockedError):
        guest.geometries["sphere"] = Sphere(radius=[99.0])

    admin.disconnect()
    guest.disconnect()


def test_admin_locked_unowned_superuser_allowed(server_auth: str) -> None:
    """Superuser can edit unowned geometries in admin-locked room."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    guest = ZnDraw(url=server_auth, room=admin.room)

    # Guest creates geometry first (room unlocked)
    guest.geometries["sphere"] = Sphere()

    # Admin locks room
    _lock_room(admin)

    # Admin edits unowned geometry → allowed
    admin.geometries["sphere"] = Sphere(radius=[99.0])

    retrieved = admin.geometries["sphere"]
    assert isinstance(retrieved, Sphere)
    assert retrieved.radius == [99.0]

    admin.disconnect()
    guest.disconnect()


# =============================================================================
# Admin Lock + Owned Geometry Tests
# =============================================================================


def test_admin_locked_owner_can_edit_own(server_auth: str) -> None:
    """Owner can edit their own geometry even when room is admin-locked."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    guest = ZnDraw(url=server_auth, room=admin.room)
    guest_id = _get_user_id(guest)

    # Guest claims geometry (room unlocked)
    guest.geometries["sphere"] = Sphere(owner=guest_id)

    # Admin locks room
    _lock_room(admin)

    # Guest edits their own geometry → still allowed
    guest.geometries["sphere"] = Sphere(owner=guest_id, radius=[5.0])

    retrieved = guest.geometries["sphere"]
    assert isinstance(retrieved, Sphere)
    assert retrieved.radius == [5.0]

    admin.disconnect()
    guest.disconnect()


def test_admin_locked_admin_can_edit_others_geometry(server_auth: str) -> None:
    """Admin can edit a geometry owned by another user in admin-locked room."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    guest = ZnDraw(url=server_auth, room=admin.room)
    guest_id = _get_user_id(guest)

    # Guest claims geometry (room unlocked)
    guest.geometries["sphere"] = Sphere(owner=guest_id)

    # Admin locks room
    _lock_room(admin)

    # Admin edits guest's geometry → allowed (superuser bypass)
    admin.geometries["sphere"] = Sphere(owner=guest_id, radius=[10.0])

    retrieved = admin.geometries["sphere"]
    assert isinstance(retrieved, Sphere)
    assert retrieved.radius == [10.0]

    admin.disconnect()
    guest.disconnect()


def test_admin_can_claim_others_geometry(server_auth: str) -> None:
    """Admin can claim a geometry owned by another user."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    guest = ZnDraw(url=server_auth, room=admin.room)
    guest_id = _get_user_id(guest)
    admin_id = _get_user_id(admin)

    # Guest claims geometry
    guest.geometries["sphere"] = Sphere(owner=guest_id)

    # Admin claims it for themselves
    admin.geometries["sphere"] = Sphere(owner=admin_id)

    retrieved = admin.geometries["sphere"]
    assert retrieved.owner == admin_id

    admin.disconnect()
    guest.disconnect()


# =============================================================================
# Admin Lock + Claiming Tests
# =============================================================================


def test_admin_locked_claiming_non_superuser_blocked(server_auth: str) -> None:
    """Non-superuser cannot claim an unowned geometry in admin-locked room."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    guest = ZnDraw(url=server_auth, room=admin.room)
    guest_id = _get_user_id(guest)

    # Guest creates unowned geometry (room unlocked)
    guest.geometries["sphere"] = Sphere()

    # Admin locks room
    _lock_room(admin)

    # Guest tries to claim → blocked (admin lock blocks all non-superuser edits on unowned)
    with pytest.raises(RoomLockedError):
        guest.geometries["sphere"] = Sphere(owner=guest_id)

    admin.disconnect()
    guest.disconnect()
