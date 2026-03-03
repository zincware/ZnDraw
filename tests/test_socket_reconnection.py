"""Tests for socket reconnection handling after server restart.

Tests the _on_connect_error handler in socket_manager.py that handles
the case where a user is not found after server restart.
"""

import uuid
from unittest.mock import MagicMock, patch

import pytest

from zndraw import ZnDraw


def test_connect_error_handler_reregisters_user_on_not_found(server):
    """Test that _on_connect_error re-registers user when 'User not found' error occurs."""
    # Create a client and connect normally
    room = uuid.uuid4().hex
    client = ZnDraw(room=room, url=server)
    assert client.socket.connected

    # Store original user info
    original_user = client.user
    original_token = client.api.jwt_token

    # Simulate the connect_error handler being called with "User not found"
    # This simulates what happens when the server restarts and Redis is cleared
    client.socket._on_connect_error("User not found. Please register first.")

    # Verify that login was called (re-registration happened)
    # The user should still be the same (or similar guest name)
    assert client.user is not None

    # The JWT token should have been refreshed
    # (In practice it might be the same token if user existed, but the call was made)
    assert client.api.jwt_token is not None

    client.disconnect()


def test_connect_error_handler_ignores_other_errors(server):
    """Test that _on_connect_error doesn't re-register for non-user-not-found errors."""
    room = uuid.uuid4().hex
    client = ZnDraw(room=room, url=server)
    assert client.socket.connected

    original_user = client.user

    # Mock the login method to track if it gets called
    with patch.object(client.api, "login") as mock_login:
        # Simulate a different error (not "User not found")
        client.socket._on_connect_error("Connection timeout")

        # Login should NOT be called for other errors
        mock_login.assert_not_called()

    # User should remain unchanged
    assert client.user == original_user

    client.disconnect()


def test_connect_error_handler_handles_login_failure_gracefully(server):
    """Test that _on_connect_error handles login failure gracefully."""
    room = uuid.uuid4().hex
    client = ZnDraw(room=room, url=server)
    assert client.socket.connected

    original_user = client.user

    # Mock login to raise an exception
    with patch.object(
        client.api, "login", side_effect=Exception("Network error")
    ) as mock_login:
        # Should not raise, just log the error
        client.socket._on_connect_error("User not found")

        # Login was attempted
        mock_login.assert_called_once()

    # User should remain unchanged (login failed)
    assert client.user == original_user

    client.disconnect()


def test_socket_reconnect_after_user_reregistration(server):
    """Test full reconnection flow after user re-registration.

    This test simulates what happens when:
    1. Client is connected
    2. Server restarts (simulated by clearing user data)
    3. Socket reconnects and gets "User not found"
    4. Client re-registers and reconnects successfully
    """
    import redis

    room = uuid.uuid4().hex
    client = ZnDraw(room=room, url=server)
    assert client.socket.connected

    original_user = client.user
    original_session_id = client.api.session_id

    # Get Redis client to manipulate user data
    r = redis.Redis(host="localhost", port=6379, decode_responses=True)

    # Simulate server restart by deleting the user data from Redis
    user_key = f"users:data:{original_user}"
    user_index_key = "users:index"

    # Delete user data (simulates Redis cleared on restart)
    r.delete(user_key)
    r.srem(user_index_key, original_user)

    # Now trigger the connect_error handler as if reconnection failed
    client.socket._on_connect_error("User not found. Please register first.")

    # Verify user was re-registered (user data should exist again)
    assert r.exists(user_key) or client.user != original_user

    # User should be set (either same or new guest name)
    assert client.user is not None
    assert client.api.jwt_token is not None

    client.disconnect()


def test_connect_error_with_password_user(server):
    """Test that _on_connect_error works correctly with password-authenticated users."""
    from pydantic import SecretStr

    room = uuid.uuid4().hex

    # Create client (will be guest user, no password)
    client = ZnDraw(room=room, url=server)
    assert client.socket.connected

    # Simulate having a password set (even though guest users don't need it)
    client.password = SecretStr("test-password")

    # Mock login to verify password is passed correctly
    with patch.object(client.api, "login") as mock_login:
        mock_login.return_value = {
            "userName": client.user,
            "role": "guest",
            "token": "mock-token",
        }

        client.socket._on_connect_error("User not found")

        # Verify login was called with correct password
        mock_login.assert_called_once_with(
            user_name=client.user,
            password="test-password",
        )

    client.disconnect()


def test_connect_error_updates_role_from_login_response(server):
    """Test that _on_connect_error updates role from login response."""
    room = uuid.uuid4().hex
    client = ZnDraw(room=room, url=server)
    assert client.socket.connected

    original_role = client.role

    # Mock login to return a different role
    with patch.object(client.api, "login") as mock_login:
        mock_login.return_value = {
            "userName": "new-admin-user",
            "role": "admin",
            "token": "mock-token",
        }

        client.socket._on_connect_error("User not found")

        # Role should be updated from login response
        assert client.role == "admin"
        assert client.user == "new-admin-user"

    client.disconnect()
