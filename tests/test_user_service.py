"""Tests for UserService."""

import pytest
from znsocket import MemoryStorage

from zndraw.services import AdminService, UserService
from zndraw.services.user_service import UserRole


@pytest.fixture
def redis_client():
    """Create a memory storage client for testing."""
    return MemoryStorage()


@pytest.fixture
def user_service(redis_client):
    """Create UserService instance for testing."""
    return UserService(redis_client)


@pytest.fixture
def admin_service(redis_client):
    """Create AdminService instance for testing."""
    return AdminService(redis_client)


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


def test_username_exists_false(user_service):
    """Test that non-existent username returns False."""
    assert not user_service.username_exists("nonexistent-user")


def test_username_exists_true(user_service):
    """Test that existing username returns True."""
    user_name = "existing-user"
    user_service.create_user(user_name)
    assert user_service.username_exists(user_name)


def test_create_user_success(user_service, redis_client):
    """Test creating a new user."""
    user_name = "new-user"
    result = user_service.create_user(user_name)

    assert result is True
    assert redis_client.exists(f"user:{user_name}")
    assert redis_client.hget(f"user:{user_name}", "userName") == user_name
    assert redis_client.hexists(f"user:{user_name}", "createdAt")
    assert redis_client.hexists(f"user:{user_name}", "lastLogin")


def test_create_user_duplicate_raises_error(user_service):
    """Test that creating duplicate user raises ValueError."""
    user_name = "duplicate-user"
    user_service.create_user(user_name)

    with pytest.raises(ValueError, match="already exists"):
        user_service.create_user(user_name)


def test_get_user_role_guest(user_service):
    """Test getting role for guest user (no password set)."""
    user_name = "guest-user"
    user_service.create_user(user_name)

    role = user_service.get_user_role(user_name)
    assert role == UserRole.GUEST


def test_get_user_role_user(user_service):
    """Test getting role for registered user (has password)."""
    user_name = "registered-user"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password123")

    role = user_service.get_user_role(user_name)
    assert role == UserRole.USER


def test_get_user_role_admin(user_service, admin_service, clear_admin_env_vars):
    """Test getting role for admin user."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret"

    admin_service_deployed = AdminService(admin_service.r)
    user_name = "admin-user"

    user_service.create_user(user_name)
    admin_service_deployed.grant_admin(user_name)

    role = user_service.get_user_role(user_name)
    assert role == UserRole.ADMIN


def test_get_user_role_nonexistent_user(user_service):
    """Test getting role for non-existent user returns GUEST."""
    role = user_service.get_user_role("nonexistent")
    assert role == UserRole.GUEST


def test_is_registered_false(user_service):
    """Test is_registered for guest user."""
    user_name = "guest-user"
    user_service.create_user(user_name)

    assert not user_service.is_registered(user_name)


def test_is_registered_true(user_service):
    """Test is_registered for registered user."""
    user_name = "registered-user"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password")

    assert user_service.is_registered(user_name)


def test_is_registered_nonexistent_user(user_service):
    """Test is_registered for non-existent user returns False."""
    assert not user_service.is_registered("nonexistent")


def test_register_user_same_username(user_service, redis_client):
    """Test registering user with same username (guest -> user)."""
    user_name = "guest-user"
    password = "mypassword123"

    user_service.create_user(user_name)
    assert not user_service.is_registered(user_name)

    result = user_service.register_user(user_name, user_name, password)
    assert result is True

    assert user_service.is_registered(user_name)
    assert redis_client.hexists(f"user:{user_name}", "passwordHash")
    assert redis_client.hexists(f"user:{user_name}", "passwordSalt")

    role = user_service.get_user_role(user_name)
    assert role == UserRole.USER


def test_register_user_new_username(user_service, redis_client):
    """Test registering with new username (guest -> user with name change)."""
    guest_name = "user-abc123"
    new_username = "CoolUsername"
    password = "mypassword123"

    user_service.create_user(guest_name)
    result = user_service.register_user(guest_name, new_username, password)
    assert result is True

    # New username should exist and be registered
    assert user_service.username_exists(new_username)
    assert user_service.is_registered(new_username)
    assert redis_client.hexists(f"user:{new_username}", "passwordHash")
    assert redis_client.hexists(f"user:{new_username}", "passwordSalt")

    # Old guest username should be deleted
    assert not user_service.username_exists(guest_name)

    role = user_service.get_user_role(new_username)
    assert role == UserRole.USER


def test_register_user_already_registered(user_service):
    """Test that registering already registered user raises error."""
    user_name = "registered-user"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password1")

    with pytest.raises(ValueError, match="already registered"):
        user_service.register_user(user_name, user_name, "password2")


def test_register_user_username_taken(user_service):
    """Test registering guest with existing username fails."""
    existing_user = "existing-user"
    guest_user = "guest-user"

    user_service.create_user(existing_user)
    user_service.register_user(existing_user, existing_user, "password")

    user_service.create_user(guest_user)

    with pytest.raises(ValueError, match="already taken"):
        user_service.register_user(guest_user, existing_user, "password")


def test_register_user_empty_username(user_service):
    """Test that empty username is rejected during registration."""
    guest_name = "guest-user"
    user_service.create_user(guest_name)

    with pytest.raises(ValueError, match="Username cannot be empty"):
        user_service.register_user(guest_name, "", "password")

    with pytest.raises(ValueError, match="Username cannot be empty"):
        user_service.register_user(guest_name, "   ", "password")


def test_register_user_preserves_timestamps(user_service, redis_client):
    """Test that registration preserves original createdAt timestamp."""
    import time

    guest_name = "guest-user"
    new_username = "registered-user"

    user_service.create_user(guest_name)
    original_created_at = redis_client.hget(f"user:{guest_name}", "createdAt")

    # Wait a bit to ensure timestamp would be different
    time.sleep(0.1)

    user_service.register_user(guest_name, new_username, "password")

    new_created_at = redis_client.hget(f"user:{new_username}", "createdAt")

    # Decode bytes if necessary
    if isinstance(original_created_at, bytes):
        original_created_at = original_created_at.decode("utf-8")
    if isinstance(new_created_at, bytes):
        new_created_at = new_created_at.decode("utf-8")

    assert new_created_at == original_created_at


def test_verify_password_correct(user_service):
    """Test password verification with correct password."""
    user_name = "user"
    password = "correctpassword"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, password)

    assert user_service.verify_password(user_name, password)


def test_verify_password_incorrect(user_service):
    """Test password verification with wrong password."""
    user_name = "user"
    password = "correctpassword"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, password)

    assert not user_service.verify_password(user_name, "wrongpassword")


def test_verify_password_not_registered(user_service):
    """Test password verification for unregistered user."""
    user_name = "guest-user"
    user_service.create_user(user_name)

    assert not user_service.verify_password(user_name, "anypassword")


def test_verify_password_nonexistent_user(user_service):
    """Test password verification for non-existent user."""
    assert not user_service.verify_password("nonexistent", "anypassword")


def test_change_password_success(user_service):
    """Test changing password successfully."""
    user_name = "user"
    old_password = "oldpass123"
    new_password = "newpass456"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, old_password)

    result = user_service.change_password(user_name, old_password, new_password)
    assert result is True

    assert not user_service.verify_password(user_name, old_password)
    assert user_service.verify_password(user_name, new_password)


def test_change_password_wrong_old_password(user_service):
    """Test that wrong old password fails."""
    user_name = "user"
    old_password = "oldpass123"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, old_password)

    with pytest.raises(ValueError, match="Current password is incorrect"):
        user_service.change_password(user_name, "wrongold", "newpass")


def test_change_password_not_registered(user_service):
    """Test that unregistered user cannot change password."""
    user_name = "guest-user"
    user_service.create_user(user_name)

    with pytest.raises(ValueError, match="not registered"):
        user_service.change_password(user_name, "old", "new")


def test_change_password_nonexistent_user(user_service):
    """Test that non-existent user cannot change password."""
    with pytest.raises(ValueError, match="not registered"):
        user_service.change_password("nonexistent", "old", "new")


def test_reset_password_admin(user_service):
    """Test admin resetting user password without old password."""
    user_name = "user"
    original_password = "original123"
    new_password = "reset456"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, original_password)

    result = user_service.reset_password(user_name, new_password)
    assert result is True

    assert not user_service.verify_password(user_name, original_password)
    assert user_service.verify_password(user_name, new_password)


def test_reset_password_not_registered(user_service):
    """Test that cannot reset password for unregistered user."""
    user_name = "guest-user"
    user_service.create_user(user_name)

    with pytest.raises(ValueError, match="not registered"):
        user_service.reset_password(user_name, "newpass")


def test_reset_password_nonexistent_user(user_service):
    """Test that cannot reset password for non-existent user."""
    with pytest.raises(ValueError, match="not registered"):
        user_service.reset_password("nonexistent", "newpass")


def test_update_last_login(user_service, redis_client):
    """Test updating last login timestamp."""
    import time

    user_name = "user"
    user_service.create_user(user_name)

    original_last_login = redis_client.hget(f"user:{user_name}", "lastLogin")

    time.sleep(0.1)

    user_service.update_last_login(user_name)

    new_last_login = redis_client.hget(f"user:{user_name}", "lastLogin")

    if isinstance(original_last_login, bytes):
        original_last_login = original_last_login.decode("utf-8")
    if isinstance(new_last_login, bytes):
        new_last_login = new_last_login.decode("utf-8")

    assert new_last_login != original_last_login


def test_delete_user(user_service, clear_admin_env_vars, redis_client):
    """Test deleting a user."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret"

    admin_service = AdminService(redis_client)
    user_name = "DeleteMe"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password")
    admin_service.grant_admin(user_name)

    assert redis_client.exists(f"user:{user_name}")
    assert admin_service.is_admin(user_name)

    result = user_service.delete_user(user_name)
    assert result is True

    assert not redis_client.exists(f"user:{user_name}")
    assert not admin_service.is_admin(user_name)


def test_delete_user_with_admin_status(
    user_service, clear_admin_env_vars, redis_client
):
    """Test deleting a user also removes admin status."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret"

    admin_service = AdminService(redis_client)
    user_name = "AdminUser"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password")
    admin_service.grant_admin(user_name)

    assert admin_service.is_admin(user_name)

    user_service.delete_user(user_name)

    assert not admin_service.is_admin(user_name)


def test_delete_nonexistent_user(user_service):
    """Test deleting non-existent user (should not raise error)."""
    result = user_service.delete_user("nonexistent")
    assert result is True


def test_list_all_users(user_service, clear_admin_env_vars, redis_client):
    """Test listing all users with different roles."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret"

    admin_service = AdminService(redis_client)

    guest_name = "GuestUser"
    user_name = "RegisteredUser"
    admin_name = "AdminUser"

    user_service.create_user(guest_name)

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password")

    user_service.create_user(admin_name)
    user_service.register_user(admin_name, admin_name, "password")
    admin_service.grant_admin(admin_name)

    users = user_service.list_all_users()

    assert len(users) == 3

    guest = next(u for u in users if u["userName"] == guest_name)
    user = next(u for u in users if u["userName"] == user_name)
    admin = next(u for u in users if u["userName"] == admin_name)

    assert guest["role"] == "guest"
    assert "createdAt" in guest
    assert "lastLogin" in guest

    assert user["role"] == "user"
    assert "createdAt" in user
    assert "lastLogin" in user

    assert admin["role"] == "admin"
    assert "createdAt" in admin
    assert "lastLogin" in admin


def test_list_all_users_empty(user_service):
    """Test listing users when none exist."""
    users = user_service.list_all_users()
    assert users == []


def test_password_hashing_different_salts(user_service, redis_client):
    """Test that same password with different salts produces different hashes."""
    user_name_1 = "user-1"
    user_name_2 = "user-2"
    same_password = "samepassword123"

    user_service.create_user(user_name_1)
    user_service.register_user(user_name_1, user_name_1, same_password)

    user_service.create_user(user_name_2)
    user_service.register_user(user_name_2, user_name_2, same_password)

    assert user_service.verify_password(user_name_1, same_password)
    assert user_service.verify_password(user_name_2, same_password)

    hash_1 = redis_client.hget(f"user:{user_name_1}", "passwordHash")
    hash_2 = redis_client.hget(f"user:{user_name_2}", "passwordHash")

    if isinstance(hash_1, bytes):
        hash_1 = hash_1.decode("utf-8")
    if isinstance(hash_2, bytes):
        hash_2 = hash_2.decode("utf-8")

    assert hash_1 != hash_2


def test_no_password_requirements(user_service):
    """Test that any password is accepted (no requirements enforced)."""
    base_name = "test-user"

    # Very short password
    user_name_1 = f"{base_name}-1"
    user_service.create_user(user_name_1)
    user_service.register_user(user_name_1, user_name_1, "1")
    assert user_service.verify_password(user_name_1, "1")

    # Password with only spaces
    user_name_2 = f"{base_name}-2"
    user_service.create_user(user_name_2)
    user_service.register_user(user_name_2, user_name_2, "   ")
    assert user_service.verify_password(user_name_2, "   ")

    # Very long password
    user_name_3 = f"{base_name}-3"
    long_pass = "x" * 1000
    user_service.create_user(user_name_3)
    user_service.register_user(user_name_3, user_name_3, long_pass)
    assert user_service.verify_password(user_name_3, long_pass)


def test_register_user_transfers_admin_status(
    user_service, clear_admin_env_vars, redis_client
):
    """Test that registering with new username transfers admin status."""
    import os

    os.environ["ZNDRAW_ADMIN_USERNAME"] = "admin"
    os.environ["ZNDRAW_ADMIN_PASSWORD"] = "secret"

    admin_service = AdminService(redis_client)

    guest_name = "guest-admin"
    new_username = "admin-user"

    user_service.create_user(guest_name)
    admin_service.grant_admin(guest_name)
    assert admin_service.is_admin(guest_name)

    user_service.register_user(guest_name, new_username, "password")

    assert not admin_service.is_admin(guest_name)
    assert admin_service.is_admin(new_username)
