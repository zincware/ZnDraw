"""Tests for protected geometry flag functionality."""

import pytest

from zndraw import ZnDraw
from zndraw.geometries import Camera, Sphere


def test_protected_geometry_cannot_be_deleted(server):
    """Protected geometries cannot be deleted by default."""
    vis = ZnDraw(url=server, room="room-protected-delete", user="tester")

    vis.geometries["protected_cam"] = Camera(protected=True)

    with pytest.raises(PermissionError, match="Cannot delete protected geometry"):
        del vis.geometries["protected_cam"]


def test_protected_geometry_can_be_deleted_after_unprotect(server):
    """Protected geometries can be deleted after setting protected=False."""
    vis = ZnDraw(url=server, room="room-protected-unprotect", user="tester")

    vis.geometries["protected_cam"] = Camera(protected=True)

    # Update to unprotect
    cam = vis.geometries["protected_cam"]
    vis.geometries["protected_cam"] = Camera(
        position=cam.position,
        target=cam.target,
        protected=False,
    )

    # Now deletion should work
    del vis.geometries["protected_cam"]
    assert "protected_cam" not in vis.geometries


def test_unprotected_geometry_can_be_deleted(server):
    """Regular geometries can be deleted normally."""
    vis = ZnDraw(url=server, room="room-unprotected-delete", user="tester")

    vis.geometries["test_sphere"] = Sphere()
    del vis.geometries["test_sphere"]
    assert "test_sphere" not in vis.geometries


def test_protected_default_is_false(server):
    """Geometries are not protected by default."""
    vis = ZnDraw(url=server, room="room-protected-default", user="tester")

    vis.geometries["test_sphere"] = Sphere()
    sphere = vis.geometries["test_sphere"]
    assert sphere.protected is False


def test_protected_camera_can_be_set(server):
    """Camera can be created with protected=True."""
    vis = ZnDraw(url=server, room="room-protected-camera", user="tester")

    vis.geometries["my_cam"] = Camera(protected=True)
    cam = vis.geometries["my_cam"]
    assert cam.protected is True


def test_protected_flag_multi_client_sync(server):
    """Protected flag is synced correctly between clients."""
    vis1 = ZnDraw(url=server, room="room-protected-sync", user="tester1")
    vis2 = ZnDraw(url=server, room="room-protected-sync", user="tester2")

    vis1.geometries["protected_cam"] = Camera(protected=True)

    # Poll until geometry is synced to client 2 (max 5 seconds)
    max_wait = 5.0
    poll_interval = 0.1
    elapsed = 0.0
    while elapsed < max_wait:
        if "protected_cam" in vis2.geometries:
            break
        vis2.socket.sio.sleep(poll_interval)
        elapsed += poll_interval
    else:
        pytest.fail(
            f"Geometry 'protected_cam' not synced to client 2 within {max_wait}s"
        )

    # Client 2 should see the protected flag
    cam = vis2.geometries["protected_cam"]
    assert cam.protected is True

    # Client 2 should not be able to delete
    with pytest.raises(PermissionError, match="Cannot delete protected geometry"):
        del vis2.geometries["protected_cam"]
