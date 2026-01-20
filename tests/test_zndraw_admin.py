"""Tests for ZnDraw admin functionality.

These tests use subprocess-based servers to avoid eventlet monkey-patching issues
that occur when running the full test suite.
"""

import pytest

from zndraw import ZnDraw

# Must match conftest.py values
TEST_ADMIN_USERNAME = "test-admin"
TEST_ADMIN_PASSWORD = "test-admin-password"


def test_zndraw_local_mode_without_password(server):
    """Test ZnDraw in local mode without password - should be admin."""
    vis = ZnDraw(url=server, room="test_local", user="Alice")

    assert vis.is_admin is True
    assert vis.role == "admin"
    assert vis.user == "Alice"

    vis.disconnect()


def test_zndraw_deployment_mode_with_correct_password(server_admin_mode):
    """Test ZnDraw in deployment mode with correct password - should be admin."""
    vis = ZnDraw(
        url=server_admin_mode,
        room="test_admin",
        user=TEST_ADMIN_USERNAME,
        password=TEST_ADMIN_PASSWORD,
    )

    assert vis.is_admin is True
    assert vis.role == "admin"
    assert vis.user == TEST_ADMIN_USERNAME

    vis.disconnect()


def test_zndraw_deployment_mode_with_wrong_password(server_admin_mode):
    """Test ZnDraw in deployment mode with wrong password - should raise error."""
    with pytest.raises(RuntimeError, match=r"Login failed.*Authentication failed"):
        ZnDraw(
            url=server_admin_mode,
            room="test_wrong",
            user=TEST_ADMIN_USERNAME,
            password="wrongpassword",
        )


def test_zndraw_deployment_mode_without_password(server_admin_mode):
    """Test ZnDraw in deployment mode without password - should NOT be admin."""
    vis = ZnDraw(url=server_admin_mode, room="test_nopass", user="Bob")

    assert vis.is_admin is False
    assert vis.role == "guest"
    assert vis.user == "Bob"

    vis.disconnect()


def test_zndraw_password_parameter_is_optional(server):
    """Test that password parameter is optional and defaults to None."""
    vis = ZnDraw(url=server, room="test_optional", user="Charlie")

    # Should work without password parameter
    assert vis.password is None
    assert vis.is_admin is True  # Local mode - everyone is admin

    vis.disconnect()
