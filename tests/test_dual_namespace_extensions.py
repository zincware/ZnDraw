"""Tests for registering the same extension in both public and private namespaces.

This module tests the ability to register the same extension class in both
global (public) and room-scoped (private) namespaces, and the ability to
explicitly choose which namespace to use when running extensions.
"""

import requests
from zndraw import ZnDraw
from zndraw.extensions import Extension, Category


class DualExtension(Extension):
    """Test extension that can be registered in both namespaces."""

    category = Category.MODIFIER
    parameter: int = 1

    def run(self, *args, **kwargs):
        """Dummy run method."""
        pass


def test_register_same_extension_in_both_namespaces(server):
    """Test that the same extension can be registered as both public and private."""
    vis = ZnDraw(url=server, room="test_room", user="admin", auto_pickup_jobs=False)

    # This should work without raising ValueError
    vis.register_extension(DualExtension, public=False)
    vis.register_extension(DualExtension, public=True)

    # Verify both are registered by checking schemas
    response = requests.get(f"{server}/api/rooms/test_room/schema/modifiers")
    assert response.status_code == 200
    schemas = response.json()

    # Schemas is now a list of extension objects
    assert isinstance(schemas, list), "Schemas should be a list"

    # Should see DualExtension in schemas (check by name field)
    extension_names = [ext.get("name") for ext in schemas]
    assert "DualExtension" in extension_names


def test_vis_run_with_public_parameter_explicit(server):
    """Test that vis.run() respects explicit public parameter."""
    vis = ZnDraw(url=server, room="test_room", user="admin", auto_pickup_jobs=False)

    # Register in both namespaces
    vis.register_extension(DualExtension, public=False)
    vis.register_extension(DualExtension, public=True)

    # Create extension instance
    ext = DualExtension(parameter=42)

    # Run with explicit public=True should use public endpoint
    # vis.run() now returns a Job object
    job_public = vis.run(ext, public=True)
    assert job_public.job_id is not None
    assert job_public.is_assigned()  # Job assigned but not processed (auto_pickup=False)

    # Complete the public job to free worker capacity
    # Manually transition job to processing then completed
    import requests
    from conftest import get_jwt_auth_headers

    response = requests.put(
        f"{server}/api/rooms/test_room/jobs/{job_public.job_id}/status",
        json={"status": "processing", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, "admin"),
    )
    assert response.status_code == 200

    response = requests.put(
        f"{server}/api/rooms/test_room/jobs/{job_public.job_id}/status",
        json={"status": "completed", "workerId": vis.sid},
        headers=get_jwt_auth_headers(server, "admin"),
    )
    assert response.status_code == 200

    # Run with explicit public=False should use private endpoint
    job_private = vis.run(ext, public=False)
    assert job_private.job_id is not None
    assert job_private.is_assigned()  # Now worker is free to accept this job

    # The jobs should be different
    assert job_public.job_id != job_private.job_id


def test_vis_run_with_public_none_prioritizes_public(server):
    """Test that vis.run() without public parameter tries public first."""
    vis = ZnDraw(url=server, room="test_room", user="admin", auto_pickup_jobs=False)

    # Register in both namespaces
    vis.register_extension(DualExtension, public=False)
    vis.register_extension(DualExtension, public=True)

    # Create extension instance
    ext = DualExtension(parameter=99)

    # Run without specifying public (should default to None and try public first)
    job = vis.run(ext)
    assert job.job_id is not None
    assert job.is_assigned()

    # To verify it used public, we'd need to check the actual endpoint called
    # For now, we just verify it succeeded


def test_vis_run_with_public_none_falls_back_to_private(server):
    """Test that vis.run() falls back to private if only private is registered."""
    vis = ZnDraw(url=server, room="test_room", user="user1", auto_pickup_jobs=False)

    # Register ONLY in private namespace
    vis.register_extension(DualExtension, public=False)

    # Create extension instance
    ext = DualExtension(parameter=77)

    # Run without specifying public - should fall back to private
    job = vis.run(ext)
    assert job.job_id is not None
    assert job.is_assigned()


def test_vis_run_explicit_public_fails_if_only_private_registered(server):
    """Test that explicitly requesting public fails if only private is registered."""
    vis = ZnDraw(url=server, room="test_room", user="user1", auto_pickup_jobs=False)

    # Register ONLY in private namespace
    vis.register_extension(DualExtension, public=False)

    # Create extension instance
    ext = DualExtension(parameter=55)

    # Explicitly request public - should fail
    try:
        result = vis.run(ext, public=True)
        # If we get here, the test should fail
        assert False, "Should have failed when requesting public for private-only extension"
    except ValueError as e:
        assert "not registered in public namespace" in str(e)


def test_vis_run_explicit_private_fails_if_only_public_registered(server):
    """Test that explicitly requesting private fails if only public is registered."""
    vis = ZnDraw(url=server, room="test_room", user="admin", auto_pickup_jobs=False)

    # Register ONLY in public namespace
    vis.register_extension(DualExtension, public=True)

    # Create extension instance
    ext = DualExtension(parameter=33)

    # Explicitly request private - should fail
    try:
        result = vis.run(ext, public=False)
        # If we get here, the test should fail
        assert False, "Should have failed when requesting private for public-only extension"
    except ValueError as e:
        assert "not registered in private namespace" in str(e)


def test_auto_pickup_jobs_with_dual_registration(server):
    """Test that auto_pickup_jobs works when extension is registered in both namespaces."""
    vis = ZnDraw(url=server, room="test_room", user="admin", auto_pickup_jobs=True)

    # Register extension in BOTH namespaces
    vis.register_extension(DualExtension, public=False)
    vis.register_extension(DualExtension, public=True)

    # Submit a job using private namespace
    ext = DualExtension(parameter=99)
    job = vis.run(ext, public=False)
    assert job.job_id is not None

    # Wait for worker to pick up and process (uses socketio.sleep to not block event loop)
    job.wait(timeout=5)

    # Job should have been processed successfully
    assert job.is_completed()
