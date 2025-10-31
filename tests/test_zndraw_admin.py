"""Tests for ZnDraw admin functionality."""

import os
import threading
import time

import pytest

from zndraw import ZnDraw
from zndraw.server import create_app, socketio


@pytest.fixture
def clear_admin_env_vars():
    """Clear admin environment variables before and after test."""
    # Save original values
    original_username = os.environ.get("ZNDRAW_ADMIN_USERNAME")
    original_password = os.environ.get("ZNDRAW_ADMIN_PASSWORD")

    # Clear before test
    if "ZNDRAW_ADMIN_USERNAME" in os.environ:
        del os.environ["ZNDRAW_ADMIN_USERNAME"]
    if "ZNDRAW_ADMIN_PASSWORD" in os.environ:
        del os.environ["ZNDRAW_ADMIN_PASSWORD"]

    yield

    # Restore after test
    if original_username is not None:
        os.environ["ZNDRAW_ADMIN_USERNAME"] = original_username
    elif "ZNDRAW_ADMIN_USERNAME" in os.environ:
        del os.environ["ZNDRAW_ADMIN_USERNAME"]

    if original_password is not None:
        os.environ["ZNDRAW_ADMIN_PASSWORD"] = original_password
    elif "ZNDRAW_ADMIN_PASSWORD" in os.environ:
        del os.environ["ZNDRAW_ADMIN_PASSWORD"]


@pytest.fixture
def local_server(clear_admin_env_vars):
    """Start a local Flask-SocketIO server in local mode."""
    app = create_app(redis_url=None)
    app.config["TESTING"] = True
    app.config["SERVER_URL"] = "http://localhost:5555"

    # Run server in a background thread
    server_thread = threading.Thread(
        target=lambda: socketio.run(app, host="127.0.0.1", port=5555, debug=False),
        daemon=True,
    )
    server_thread.start()
    time.sleep(1)  # Give server time to start

    yield "http://localhost:5555"

    # Server will be stopped when test ends (daemon thread)


@pytest.fixture
def deployment_server(clear_admin_env_vars):
    """Start a local Flask-SocketIO server in deployment mode."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    app = create_app(redis_url=None)
    app.config["TESTING"] = True
    app.config["SERVER_URL"] = "http://localhost:5556"

    # Run server in a background thread
    server_thread = threading.Thread(
        target=lambda: socketio.run(app, host="127.0.0.1", port=5556, debug=False),
        daemon=True,
    )
    server_thread.start()
    time.sleep(1)  # Give server time to start

    yield "http://localhost:5556"

    # Server will be stopped when test ends (daemon thread)


def test_zndraw_local_mode_without_password(local_server):
    """Test ZnDraw in local mode without password - should be admin."""
    vis = ZnDraw(url=local_server, room="test_local", user="Alice")

    assert vis.is_admin is True
    assert vis.role == "admin"
    assert vis.user == "Alice"

    vis.disconnect()


def test_zndraw_deployment_mode_with_correct_password(deployment_server):
    """Test ZnDraw in deployment mode with correct password - should be admin."""
    vis = ZnDraw(
        url=deployment_server, room="test_admin", user="admin", password="secret123"
    )

    assert vis.is_admin is True
    assert vis.role == "admin"
    assert vis.user == "admin"

    vis.disconnect()


def test_zndraw_deployment_mode_with_wrong_password(deployment_server):
    """Test ZnDraw in deployment mode with wrong password - should raise error."""
    with pytest.raises(
        RuntimeError, match="Login failed.*Invalid username or password"
    ):
        ZnDraw(
            url=deployment_server,
            room="test_wrong",
            user="admin",
            password="wrongpassword",
        )


def test_zndraw_deployment_mode_without_password(deployment_server):
    """Test ZnDraw in deployment mode without password - should NOT be admin."""
    vis = ZnDraw(url=deployment_server, room="test_nopass", user="Bob")

    assert vis.is_admin is False
    assert vis.role == "guest"
    assert vis.user == "Bob"

    vis.disconnect()


def test_zndraw_password_parameter_is_optional(local_server):
    """Test that password parameter is optional and defaults to None."""
    vis = ZnDraw(url=local_server, room="test_optional", user="Charlie")

    # Should work without password parameter
    assert vis.password is None
    assert vis.is_admin is True  # Local mode - everyone is admin

    vis.disconnect()
