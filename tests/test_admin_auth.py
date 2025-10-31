"""Tests for admin authentication and @require_admin decorator."""

import os

import pytest

from zndraw.server import create_app


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
def app_local_mode(clear_admin_env_vars):
    """Create Flask app in local mode (no admin credentials)."""
    app = create_app(redis_url=None)
    app.config["TESTING"] = True
    return app


@pytest.fixture
def app_deployment_mode(clear_admin_env_vars):
    """Create Flask app in deployment mode (admin credentials set)."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    app = create_app(redis_url=None)
    app.config["TESTING"] = True
    return app


def test_login_local_mode_all_users_admin(app_local_mode):
    """Test that all users are admin in local mode."""
    client = app_local_mode.test_client()

    # Login as any user
    response = client.post("/api/login", json={"userName": "John Doe"})
    assert response.status_code == 200
    data = response.get_json()

    assert data["status"] == "ok"
    assert data["role"] == "admin"
    assert "token" in data
    assert "userName" in data
    assert data["userName"] == "John Doe"


def test_login_deployment_mode_admin_credentials(app_deployment_mode):
    """Test login with correct admin credentials in deployment mode."""
    client = app_deployment_mode.test_client()

    # Login with correct admin credentials
    response = client.post(
        "/api/login", json={"userName": "admin", "password": "secret123"}
    )
    assert response.status_code == 200
    data = response.get_json()

    assert data["status"] == "ok"
    assert data["role"] == "admin"
    assert "token" in data


def test_login_deployment_mode_wrong_credentials(app_deployment_mode):
    """Test login with wrong credentials in deployment mode."""
    client = app_deployment_mode.test_client()

    # Login with wrong password - should be rejected
    response = client.post(
        "/api/login", json={"userName": "admin", "password": "wrong"}
    )
    assert response.status_code == 401
    data = response.get_json()

    assert "error" in data
    assert data["error"] == "Invalid username or password"


def test_login_deployment_mode_non_admin_user(app_deployment_mode):
    """Test login as non-admin user in deployment mode."""
    client = app_deployment_mode.test_client()

    # Login without password (non-admin user)
    response = client.post("/api/login", json={"userName": "John Doe"})
    assert response.status_code == 200
    data = response.get_json()

    assert data["status"] == "ok"
    assert data["role"] == "guest"
    assert "token" in data


def test_shutdown_requires_auth(app_local_mode):
    """Test shutdown endpoint requires authentication."""
    client = app_local_mode.test_client()

    # Try to shutdown without token
    response = client.post("/api/shutdown")
    assert response.status_code == 401
    data = response.get_json()
    assert "error" in data
    assert data["type"] == "AuthError"


def test_shutdown_local_mode_any_user(app_local_mode):
    """Test shutdown works for any authenticated user in local mode."""
    client = app_local_mode.test_client()

    # Login
    login_response = client.post("/api/login", json={"userName": "John Doe"})
    token = login_response.get_json()["token"]

    # Try to shutdown with valid token
    # socketio.stop() raises SystemExit in test mode, which is expected
    with pytest.raises(SystemExit):
        response = client.post(
            "/api/shutdown", headers={"Authorization": f"Bearer {token}"}
        )


def test_shutdown_deployment_mode_admin_only(app_deployment_mode):
    """Test shutdown only works for admin in deployment mode."""
    client = app_deployment_mode.test_client()

    # Login as admin
    admin_response = client.post(
        "/api/login", json={"userName": "admin", "password": "secret123"}
    )
    admin_token = admin_response.get_json()["token"]

    # Admin can shutdown
    # socketio.stop() raises SystemExit in test mode, which is expected
    with pytest.raises(SystemExit):
        response = client.post(
            "/api/shutdown", headers={"Authorization": f"Bearer {admin_token}"}
        )


def test_shutdown_deployment_mode_non_admin_forbidden(app_deployment_mode):
    """Test non-admin user cannot shutdown in deployment mode."""
    client = app_deployment_mode.test_client()

    # Login as non-admin user
    user_response = client.post("/api/login", json={"userName": "John Doe"})
    user_token = user_response.get_json()["token"]

    # Non-admin cannot shutdown
    response = client.post(
        "/api/shutdown", headers={"Authorization": f"Bearer {user_token}"}
    )
    assert response.status_code == 403
    data = response.get_json()
    assert "error" in data
    assert data["type"] == "AdminAccessError"
    assert "Admin access required" in data["error"]


def test_login_missing_username_creates_anonymous_guest(app_local_mode):
    """Test login without username creates anonymous guest."""
    client = app_local_mode.test_client()

    response = client.post("/api/login", json={})
    assert response.status_code == 200
    data = response.get_json()
    assert data["status"] == "ok"
    # Anonymous guest gets auto-generated username
    assert "userName" in data
    assert data["userName"].startswith("user-")
    # In local mode, all users are admin
    assert data["role"] == "admin"
