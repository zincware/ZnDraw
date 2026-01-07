"""Tests for ZnDraw admin functionality.

These tests use subprocess-based servers to avoid eventlet monkey-patching issues
that occur when running the full test suite.
"""

import os
import shutil
import signal
import subprocess

import pytest
import redis

from zndraw import ZnDraw
from zndraw.server_manager import remove_server_info


@pytest.fixture
def local_server(tmp_path, get_free_port, wait_for_server):
    """Start a zndraw server in local mode (no admin credentials) via subprocess."""
    from zndraw import config as config_module

    # Reset config singleton
    config_module._config = None

    port = get_free_port()
    storage_path = tmp_path / "zndraw-local"
    redis_url = "redis://localhost:6379"

    # Create clean environment without admin credentials (local mode)
    env = os.environ.copy()
    env.pop("ZNDRAW_ADMIN_USERNAME", None)
    env.pop("ZNDRAW_ADMIN_PASSWORD", None)

    # Start server via subprocess (unique port ensures new server starts)
    proc = subprocess.Popen(
        [
            "zndraw",
            "--port",
            str(port),
            "--no-celery",
            "--storage-path",
            str(storage_path),
            "--redis-url",
            redis_url,
            "--no-browser",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )

    # Wait for server to be ready
    if not wait_for_server(port):
        proc.kill()
        raise TimeoutError(f"Local server did not start on port {port}")

    try:
        yield f"http://127.0.0.1:{port}"
    finally:
        proc.send_signal(signal.SIGTERM)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()

        # Clean up
        shutil.rmtree(storage_path, ignore_errors=True)
        r = redis.Redis.from_url(redis_url, decode_responses=True)
        r.flushall()
        remove_server_info(port)
        config_module._config = None


@pytest.fixture
def deployment_server(tmp_path, get_free_port, wait_for_server):
    """Start a zndraw server in deployment mode (with admin credentials) via subprocess."""
    from zndraw import config as config_module

    # Reset config singleton
    config_module._config = None

    port = get_free_port()
    storage_path = tmp_path / "zndraw-deployment"
    redis_url = "redis://localhost:6379"

    # Create environment with admin credentials (deployment mode)
    env = os.environ.copy()
    env["ZNDRAW_ADMIN_USERNAME"] = "admin"
    env["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    # Start server via subprocess (unique port ensures new server starts)
    proc = subprocess.Popen(
        [
            "zndraw",
            "--port",
            str(port),
            "--no-celery",
            "--storage-path",
            str(storage_path),
            "--redis-url",
            redis_url,
            "--no-browser",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env=env,
    )

    # Wait for server to be ready
    if not wait_for_server(port):
        proc.kill()
        raise TimeoutError(f"Deployment server did not start on port {port}")

    try:
        yield f"http://127.0.0.1:{port}"
    finally:
        proc.send_signal(signal.SIGTERM)
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()

        # Clean up
        shutil.rmtree(storage_path, ignore_errors=True)
        r = redis.Redis.from_url(redis_url, decode_responses=True)
        r.flushall()
        remove_server_info(port)
        config_module._config = None


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
    with pytest.raises(RuntimeError, match=r"Login failed.*Authentication failed"):
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
