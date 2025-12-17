"""Tests for AdminService."""

import pytest
from znsocket import MemoryStorage

from zndraw.services import AdminService


@pytest.fixture
def redis_client():
    """Create a memory storage client for testing."""
    return MemoryStorage()


def test_admin_service_local_mode(redis_client):
    """Test AdminService in local mode (no admin credentials configured)."""
    service = AdminService(redis_client)

    assert not service.is_deployment_mode()

    # In local mode, all users are considered admin
    assert service.is_admin("any-username")
    assert service.is_admin("another-username")
    assert service.is_admin("guest-user-123")


def test_admin_service_deployment_mode(redis_client):
    """Test AdminService in deployment mode (admin credentials set)."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    assert service.is_deployment_mode()

    # No users are admin by default in deployment mode
    assert not service.is_admin("regular-user")
    assert not service.is_admin("another-user")


def test_admin_service_validate_credentials_correct(redis_client):
    """Test credential validation with correct credentials."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    assert service.validate_admin_credentials("admin", "secret123")


def test_admin_service_validate_credentials_wrong_username(redis_client):
    """Test credential validation with wrong username."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    assert not service.validate_admin_credentials("notadmin", "secret123")
    assert not service.validate_admin_credentials("", "secret123")


def test_admin_service_validate_credentials_wrong_password(redis_client):
    """Test credential validation with wrong password."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    assert not service.validate_admin_credentials("admin", "wrong")
    assert not service.validate_admin_credentials("admin", "")
    assert not service.validate_admin_credentials("admin", "secret124")


def test_admin_service_validate_credentials_both_wrong(redis_client):
    """Test credential validation with both username and password wrong."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    assert not service.validate_admin_credentials("notadmin", "wrong")
    assert not service.validate_admin_credentials("", "")


def test_admin_service_grant_revoke(redis_client):
    """Test granting and revoking admin privileges to userName."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    user_name = "test-user"

    # Initially not admin
    assert not service.is_admin(user_name)

    # Grant admin
    service.grant_admin(user_name)
    assert service.is_admin(user_name)

    # Verify Redis key structure is correct
    from zndraw.app.redis_keys import UserKeys

    keys = UserKeys(user_name)
    assert redis_client.get(keys.admin_key()) in ("1", b"1")

    # Revoke admin
    service.revoke_admin(user_name)
    assert not service.is_admin(user_name)

    # Verify Redis key was deleted
    assert not redis_client.exists(keys.admin_key())


def test_admin_service_grant_multiple_admins(redis_client):
    """Test granting admin to multiple users."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    # Grant admin to multiple users
    service.grant_admin("alice")
    service.grant_admin("bob")
    service.grant_admin("charlie")

    assert service.is_admin("alice")
    assert service.is_admin("bob")
    assert service.is_admin("charlie")

    # Non-admin should still not be admin
    assert not service.is_admin("dave")


def test_admin_service_grant_idempotent(redis_client):
    """Test that granting admin multiple times is idempotent."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    user_name = "test-user"

    service.grant_admin(user_name)
    assert service.is_admin(user_name)

    # Grant again - should not cause issues
    service.grant_admin(user_name)
    assert service.is_admin(user_name)


def test_admin_service_revoke_non_admin(redis_client):
    """Test revoking admin from user who isn't admin (should not error)."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    user_name = "test-user"

    assert not service.is_admin(user_name)

    # Revoking non-admin should not raise error
    service.revoke_admin(user_name)
    assert not service.is_admin(user_name)


def test_admin_service_get_all_admins(redis_client):
    """Test getting all admin users."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    # Grant admin to multiple users
    service.grant_admin("alice")
    service.grant_admin("bob")
    service.grant_admin("charlie")

    admins = service.get_all_admins()
    assert len(admins) == 3
    assert "alice" in admins
    assert "bob" in admins
    assert "charlie" in admins


def test_admin_service_get_all_admins_empty(redis_client):
    """Test getting all admins when none exist."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    admins = service.get_all_admins()
    assert len(admins) == 0


def test_admin_service_get_all_admins_after_revoke(redis_client):
    """Test getting all admins after revoking some."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

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


def test_admin_service_validate_credentials_in_local_mode(redis_client):
    """Test that validate_admin_credentials returns False in local mode."""
    service = AdminService(redis_client)

    # In local mode, validation should always return False
    # (because admin mode is not enforced - everyone is admin)
    assert not service.validate_admin_credentials("any", "credentials")
    assert not service.validate_admin_credentials("admin", "password")


def test_admin_service_local_mode_grant_revoke_ignored(redis_client):
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


def test_admin_service_deployment_mode_redis_key_format(redis_client):
    """Test that admin status uses correct Redis key format."""
    from zndraw.app.redis_keys import UserKeys

    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    user_name = "test-user"
    keys = UserKeys(user_name)

    service.grant_admin(user_name)

    # Should use users:admin:{username} format
    assert redis_client.exists(keys.admin_key())
    assert not redis_client.exists(f"admin:client:{user_name}")


def test_admin_service_case_sensitive_username(redis_client):
    """Test that usernames are case-sensitive for admin checks."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

    service.grant_admin("Alice")

    assert service.is_admin("Alice")
    assert not service.is_admin("alice")  # Different case
    assert not service.is_admin("ALICE")  # Different case


def test_admin_service_special_characters_in_username(redis_client):
    """Test admin service with special characters in username."""
    service = AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret123",
    )

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
