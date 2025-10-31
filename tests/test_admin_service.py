"""Tests for AdminService."""

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
    import os

    original_username = os.environ.get("ZNDRAW_ADMIN_USERNAME")
    original_password = os.environ.get("ZNDRAW_ADMIN_PASSWORD")

    if "ZNDRAW_ADMIN_USERNAME" in os.environ:
        del os.environ["ZNDRAW_ADMIN_USERNAME"]
    if "ZNDRAW_ADMIN_PASSWORD" in os.environ:
        del os.environ["ZNDRAW_ADMIN_PASSWORD"]

    yield

    if original_username is not None:
        os.environ["ZNDRAW_ADMIN_USERNAME"] = original_username
    elif "ZNDRAW_ADMIN_USERNAME" in os.environ:
        del os.environ["ZNDRAW_ADMIN_USERNAME"]

    if original_password is not None:
        os.environ["ZNDRAW_ADMIN_PASSWORD"] = original_password
    elif "ZNDRAW_ADMIN_PASSWORD" in os.environ:
        del os.environ["ZNDRAW_ADMIN_PASSWORD"]


def test_admin_service_local_mode(redis_client, clear_admin_env_vars):
    """Test AdminService in local mode (no admin credentials configured)."""
    service = AdminService(redis_client)

    assert not service.is_deployment_mode()

    # In local mode, all users are considered admin
    assert service.is_admin("any-username")
    assert service.is_admin("another-username")
    assert service.is_admin("guest-user-123")


def test_admin_service_deployment_mode(redis_client, clear_admin_env_vars):
    """Test AdminService in deployment mode (admin credentials set)."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    assert service.is_deployment_mode()

    # No users are admin by default in deployment mode
    assert not service.is_admin("regular-user")
    assert not service.is_admin("another-user")


def test_admin_service_validate_credentials_correct(redis_client, clear_admin_env_vars):
    """Test credential validation with correct credentials."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    assert service.validate_admin_credentials("admin", "secret123")


def test_admin_service_validate_credentials_wrong_username(
    redis_client, clear_admin_env_vars
):
    """Test credential validation with wrong username."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    assert not service.validate_admin_credentials("notadmin", "secret123")
    assert not service.validate_admin_credentials("", "secret123")


def test_admin_service_validate_credentials_wrong_password(
    redis_client, clear_admin_env_vars
):
    """Test credential validation with wrong password."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    assert not service.validate_admin_credentials("admin", "wrong")
    assert not service.validate_admin_credentials("admin", "")
    assert not service.validate_admin_credentials("admin", "secret124")


def test_admin_service_validate_credentials_both_wrong(
    redis_client, clear_admin_env_vars
):
    """Test credential validation with both username and password wrong."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    assert not service.validate_admin_credentials("notadmin", "wrong")
    assert not service.validate_admin_credentials("", "")


def test_admin_service_grant_revoke(redis_client, clear_admin_env_vars):
    """Test granting and revoking admin privileges to userName."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    user_name = "test-user"

    # Initially not admin
    assert not service.is_admin(user_name)

    # Grant admin
    service.grant_admin(user_name)
    assert service.is_admin(user_name)

    # Verify Redis key structure is correct
    assert redis_client.get(f"admin:user:{user_name}") in ("1", b"1")

    # Revoke admin
    service.revoke_admin(user_name)
    assert not service.is_admin(user_name)

    # Verify Redis key was deleted
    assert not redis_client.exists(f"admin:user:{user_name}")


def test_admin_service_grant_multiple_admins(redis_client, clear_admin_env_vars):
    """Test granting admin to multiple users."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    # Grant admin to multiple users
    service.grant_admin("alice")
    service.grant_admin("bob")
    service.grant_admin("charlie")

    assert service.is_admin("alice")
    assert service.is_admin("bob")
    assert service.is_admin("charlie")

    # Non-admin should still not be admin
    assert not service.is_admin("dave")


def test_admin_service_grant_idempotent(redis_client, clear_admin_env_vars):
    """Test that granting admin multiple times is idempotent."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    user_name = "test-user"

    service.grant_admin(user_name)
    assert service.is_admin(user_name)

    # Grant again - should not cause issues
    service.grant_admin(user_name)
    assert service.is_admin(user_name)


def test_admin_service_revoke_non_admin(redis_client, clear_admin_env_vars):
    """Test revoking admin from user who isn't admin (should not error)."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    user_name = "test-user"

    assert not service.is_admin(user_name)

    # Revoking non-admin should not raise error
    service.revoke_admin(user_name)
    assert not service.is_admin(user_name)


def test_admin_service_get_all_admins(redis_client, clear_admin_env_vars):
    """Test getting all admin users."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    # Grant admin to multiple users
    service.grant_admin("alice")
    service.grant_admin("bob")
    service.grant_admin("charlie")

    admins = service.get_all_admins()
    assert len(admins) == 3
    assert "alice" in admins
    assert "bob" in admins
    assert "charlie" in admins


def test_admin_service_get_all_admins_empty(redis_client, clear_admin_env_vars):
    """Test getting all admins when none exist."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    admins = service.get_all_admins()
    assert len(admins) == 0


def test_admin_service_get_all_admins_after_revoke(redis_client, clear_admin_env_vars):
    """Test getting all admins after revoking some."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    service.grant_admin("alice")
    service.grant_admin("bob")
    service.grant_admin("charlie")

    # Revoke bob's admin
    service.revoke_admin("bob")

    admins = service.get_all_admins()
    assert len(admins) == 2
    assert "alice" in admins
    assert "charlie" in admins
    assert "bob" not in admins


def test_admin_service_env_validation_only_username(redis_client, clear_admin_env_vars):
    """Test that only setting username raises error."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    # Password not set

    with pytest.raises(
        ValueError, match="Both ZNDRAW_ADMIN_USERNAME and ZNDRAW_ADMIN_PASSWORD"
    ):
        AdminService(redis_client)


def test_admin_service_env_validation_only_password(redis_client, clear_admin_env_vars):
    """Test that only setting password raises error."""
    import os

    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"
    # Username not set

    with pytest.raises(
        ValueError, match="Both ZNDRAW_ADMIN_USERNAME and ZNDRAW_ADMIN_PASSWORD"
    ):
        AdminService(redis_client)


def test_admin_service_validate_credentials_in_local_mode(
    redis_client, clear_admin_env_vars
):
    """Test that validate_admin_credentials returns False in local mode."""
    service = AdminService(redis_client)

    # In local mode, validation should always return False
    # (because admin mode is not enforced - everyone is admin)
    assert not service.validate_admin_credentials("any", "credentials")
    assert not service.validate_admin_credentials("admin", "password")


def test_admin_service_local_mode_grant_revoke_ignored(
    redis_client, clear_admin_env_vars
):
    """Test that grant/revoke work even in local mode (for consistency)."""
    service = AdminService(redis_client)

    user_name = "test-user"

    # In local mode, is_admin always returns True
    assert service.is_admin(user_name)

    # But grant should still work (for when switching to deployment mode)
    service.grant_admin(user_name)
    assert service.is_admin(user_name)

    # And revoke should work
    service.revoke_admin(user_name)

    # But is_admin still returns True because in local mode
    assert service.is_admin(user_name)


def test_admin_service_deployment_mode_redis_key_format(
    redis_client, clear_admin_env_vars
):
    """Test that admin status uses correct Redis key format."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    user_name = "test-user"

    service.grant_admin(user_name)

    # Should use admin:user:<userName> format (NOT admin:client:<clientId>)
    assert redis_client.exists(f"admin:user:{user_name}")
    assert not redis_client.exists(f"admin:client:{user_name}")


def test_admin_service_case_sensitive_username(redis_client, clear_admin_env_vars):
    """Test that usernames are case-sensitive for admin checks."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    service.grant_admin("Alice")

    assert service.is_admin("Alice")
    assert not service.is_admin("alice")  # Different case
    assert not service.is_admin("ALICE")  # Different case


def test_admin_service_special_characters_in_username(
    redis_client, clear_admin_env_vars
):
    """Test admin service with special characters in username."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret123"

    service = AdminService(redis_client)

    # Usernames with special characters
    special_usernames = [
        "user-with-dashes",
        "user_with_underscores",
        "user.with.dots",
        "user@email",
        "user123",
    ]

    for user_name in special_usernames:
        service.grant_admin(user_name)
        assert service.is_admin(user_name)

        service.revoke_admin(user_name)
        assert not service.is_admin(user_name)
