"""E2E tests for register_job / register_extension convenience API.

Covers:
- register_job(cls, room='@global')
- register_extension(cls, public=True) deprecated shim
- register_job(cls, public=True) deprecated parameter
- ValueError when both room and public are specified
- Global extension visibility from any room
- Auth guards (guest cannot, admin can register global)
- register_extension(cls) room-scoped default path
- register_extension(cls, unexpected_kwarg=...) raises TypeError
"""

import warnings

import pytest

from zndraw import GLOBAL_ROOM, ZnDraw

# =============================================================================
# register_job with room='@global'
# =============================================================================


def test_register_job_global(server, Echo, get_job_list):
    """register_job(cls, room='@global') registers a global job."""
    worker = ZnDraw(url=server)
    try:
        worker.register_job(Echo, room=GLOBAL_ROOM)
        jobs = get_job_list(worker, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


# =============================================================================
# Deprecated register_extension shim
# =============================================================================


def test_register_extension_public(server, Echo, get_job_list):
    """Deprecated register_extension(cls, public=True) registers globally."""
    worker = ZnDraw(url=server)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            worker.register_extension(Echo, public=True)
            assert any(issubclass(x.category, DeprecationWarning) for x in w)
        jobs = get_job_list(worker, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


# =============================================================================
# Deprecated public= parameter on register_job
# =============================================================================


def test_register_job_public_deprecated(server, Echo, get_job_list):
    """register_job(cls, public=True) works but emits DeprecationWarning."""
    worker = ZnDraw(url=server)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            worker.register_job(Echo, public=True)
            dep_warnings = [x for x in w if issubclass(x.category, DeprecationWarning)]
            assert len(dep_warnings) >= 1
            assert "room='@global'" in str(dep_warnings[0].message)
        jobs = get_job_list(worker, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


# =============================================================================
# ValueError when both room and public are set
# =============================================================================


@pytest.mark.parametrize("public", [True, False])
def test_register_job_room_and_public_raises(server, Echo, public):
    """Passing both room= and public= (any value) raises ValueError."""
    vis = ZnDraw(url=server)
    try:
        with pytest.raises(ValueError, match="Cannot specify both"):
            vis.register_job(Echo, room=vis.room, public=public)
    finally:
        vis.disconnect()


# =============================================================================
# Deprecated register_extension — room-scoped default
# =============================================================================


def test_register_extension_room_scoped(server, Echo, get_job_list):
    """register_extension(cls) without public registers in the worker's room."""
    worker = ZnDraw(url=server)
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            worker.register_extension(Echo)
            assert any(issubclass(x.category, DeprecationWarning) for x in w)
        jobs = get_job_list(worker, room_id=worker.room)
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        worker.jobs.disconnect()
        worker.disconnect()


# =============================================================================
# register_extension — unexpected kwargs
# =============================================================================


def test_register_extension_unexpected_kwargs_raises(server, Echo):
    """register_extension(cls, foo='bar') raises TypeError."""
    worker = ZnDraw(url=server)
    try:
        with pytest.raises(TypeError, match="Unexpected keyword argument"):
            worker.register_extension(Echo, foo="bar")
    finally:
        worker.disconnect()


# =============================================================================
# Global extension visible in all rooms
# =============================================================================


def test_global_extension_visible_in_all_rooms(server, Echo, get_job_list):
    """A @global extension registered via register_job is visible from any room."""
    registrar = ZnDraw(url=server)
    room_a = ZnDraw(url=server)
    room_b = ZnDraw(url=server)
    try:
        registrar.register_job(Echo, room=GLOBAL_ROOM)

        # Verify visible from room_a
        jobs_a = get_job_list(room_a, room_id=room_a.room)
        global_a = [j for j in jobs_a if j.full_name.startswith("@global")]
        names_a = {j.name for j in global_a}
        assert "Echo" in names_a

        # Verify visible from room_b (different room)
        jobs_b = get_job_list(room_b, room_id=room_b.room)
        global_b = [j for j in jobs_b if j.full_name.startswith("@global")]
        names_b = {j.name for j in global_b}
        assert "Echo" in names_b
    finally:
        registrar.jobs.disconnect()
        registrar.disconnect()
        room_a.disconnect()
        room_b.disconnect()


# =============================================================================
# Auth Mode — register_job API
# =============================================================================


def test_guest_cannot_register_global_via_register_job(server_auth, Echo):
    """Guest user gets 403 when using register_job(room='@global')."""
    guest = ZnDraw(url=server_auth)
    try:
        with pytest.raises(PermissionError):
            guest.register_job(Echo, room=GLOBAL_ROOM)
    finally:
        guest.disconnect()


def test_admin_can_register_global_via_register_job(server_auth, Echo, get_job_list):
    """Admin user can register global jobs via register_job."""
    admin = ZnDraw(url=server_auth, user="admin@local.test", password="adminpassword")
    try:
        admin.register_job(Echo, room=GLOBAL_ROOM)
        jobs = get_job_list(admin, room_id="@global")
        names = {j.name for j in jobs}
        assert "Echo" in names
    finally:
        admin.jobs.disconnect()
        admin.disconnect()
