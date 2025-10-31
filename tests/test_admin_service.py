"""Tests for AdminService."""

import os

import pytest
from znsocket import MemoryStorage

from zndraw.services import AdminService


@pytest.fixture
def redis_client():
    """Create a memory storage client for testing."""
    return MemoryStorage()


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


def test_admin_service_local_mode(redis_client, clear_admin_env_vars):
    """Test AdminService in local mode (no admin credentials)."""
    service = AdminService(redis_client)

    # Should be in local mode
    assert not service.is_deployment_mode()

    # All users should be admin in local mode
    assert service.is_admin("any-client-id")
    assert service.is_admin("another-client-id")


def test_admin_service_deployment_mode(redis_client, clear_admin_env_vars):
    """Test AdminService in deployment mode (admin credentials set)."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    # Should be in deployment mode
    assert service.is_deployment_mode()

    # No users are admin by default in deployment mode
    assert not service.is_admin("client-1")


def test_admin_service_validate_credentials(redis_client, clear_admin_env_vars):
    """Test credential validation in deployment mode."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    # Correct credentials
    assert service.validate_admin_credentials("admin", "secret123")

    # Wrong username
    assert not service.validate_admin_credentials("notadmin", "secret123")

    # Wrong password
    assert not service.validate_admin_credentials("admin", "wrong")

    # Both wrong
    assert not service.validate_admin_credentials("notadmin", "wrong")


def test_admin_service_grant_revoke(redis_client, clear_admin_env_vars):
    """Test granting and revoking admin privileges."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    client_id = "test-client-123"

    # Initially not admin
    assert not service.is_admin(client_id)

    # Grant admin
    service.grant_admin(client_id)
    assert service.is_admin(client_id)

    # Revoke admin
    service.revoke_admin(client_id)
    assert not service.is_admin(client_id)


def test_admin_service_get_all_admins(redis_client, clear_admin_env_vars):
    """Test getting all admin users."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    # Grant admin to multiple clients
    service.grant_admin("client-1")
    service.grant_admin("client-2")
    service.grant_admin("client-3")

    admins = service.get_all_admins()
    assert len(admins) == 3
    assert "client-1" in admins
    assert "client-2" in admins
    assert "client-3" in admins


def test_admin_service_env_validation_only_username(redis_client, clear_admin_env_vars):
    """Test that only setting username raises error."""
    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    # Password not set

    with pytest.raises(ValueError, match="Both ZNDRAW_ADMIN_USERNAME and ZNDRAW_ADMIN_PASSWORD"):
        AdminService(redis_client)


def test_admin_service_env_validation_only_password(redis_client, clear_admin_env_vars):
    """Test that only setting password raises error."""
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"
    # Username not set

    with pytest.raises(ValueError, match="Both ZNDRAW_ADMIN_USERNAME and ZNDRAW_ADMIN_PASSWORD"):
        AdminService(redis_client)


def test_admin_service_validate_credentials_in_local_mode(redis_client, clear_admin_env_vars):
    """Test that validate_admin_credentials returns False in local mode."""
    service = AdminService(redis_client)

    # In local mode, validation should always return False
    # (because admin mode is not enforced)
    assert not service.validate_admin_credentials("any", "credentials")
