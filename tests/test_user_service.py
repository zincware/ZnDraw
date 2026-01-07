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
    """Create AdminService instance in local mode for testing."""
    return AdminService(redis_client)


@pytest.fixture
def admin_service_deployment(redis_client):
    """Create AdminService instance in deployment mode for testing."""
    return AdminService(
        redis_client,
        admin_username="admin",
        admin_password="secret",
    )


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
    from zndraw.app.redis_keys import UserKeys

    user_name = "new-user"
    result = user_service.create_user(user_name)
    keys = UserKeys(user_name)

    assert result is True
    assert redis_client.exists(keys.hash_key())
    assert redis_client.hget(keys.hash_key(), "userName") == user_name
    assert redis_client.hexists(keys.hash_key(), "createdAt")
    assert redis_client.hexists(keys.hash_key(), "lastLogin")


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


def test_get_user_role_admin(user_service, admin_service_deployment):
    """Test getting role for admin user."""
    user_name = "admin-user"

    user_service.create_user(user_name)
    admin_service_deployment.grant_admin(user_name)

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
    from zndraw.app.redis_keys import UserKeys

    user_name = "guest-user"
    password = "mypassword123"
    keys = UserKeys(user_name)

    user_service.create_user(user_name)
    assert not user_service.is_registered(user_name)

    result = user_service.register_user(user_name, user_name, password)
    assert result is True

    assert user_service.is_registered(user_name)
    assert redis_client.hexists(keys.hash_key(), "passwordHash")
    # Note: Argon2 stores salt in the hash itself, no separate passwordSalt field

    role = user_service.get_user_role(user_name)
    assert role == UserRole.USER


def test_register_user_new_username(user_service, redis_client):
    """Test registering with new username (guest -> user with name change)."""
    from zndraw.app.redis_keys import UserKeys

    guest_name = "user-abc123"
    new_username = "CoolUsername"
    password = "mypassword123"
    new_keys = UserKeys(new_username)

    user_service.create_user(guest_name)
    result = user_service.register_user(guest_name, new_username, password)
    assert result is True

    # New username should exist and be registered
    assert user_service.username_exists(new_username)
    assert user_service.is_registered(new_username)
    assert redis_client.hexists(new_keys.hash_key(), "passwordHash")
    # Note: Argon2 stores salt in the hash itself, no separate passwordSalt field

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

    from zndraw.app.redis_keys import UserKeys

    guest_name = "guest-user"
    new_username = "registered-user"
    guest_keys = UserKeys(guest_name)
    new_keys = UserKeys(new_username)

    user_service.create_user(guest_name)
    original_created_at = redis_client.hget(guest_keys.hash_key(), "createdAt")

    # Wait a bit to ensure timestamp would be different
    time.sleep(0.1)

    user_service.register_user(guest_name, new_username, "password")

    new_created_at = redis_client.hget(new_keys.hash_key(), "createdAt")

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
        user_service.change_password(user_name, "wrongold1", "newpass123")


def test_change_password_not_registered(user_service):
    """Test that unregistered user cannot change password."""
    user_name = "guest-user"
    user_service.create_user(user_name)

    with pytest.raises(ValueError, match="not registered"):
        user_service.change_password(user_name, "oldpass123", "newpass123")


def test_change_password_nonexistent_user(user_service):
    """Test that non-existent user cannot change password."""
    with pytest.raises(ValueError, match="not registered"):
        user_service.change_password("nonexistent", "oldpass123", "newpass123")


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
        user_service.reset_password(user_name, "newpass123")


def test_reset_password_nonexistent_user(user_service):
    """Test that cannot reset password for non-existent user."""
    with pytest.raises(ValueError, match="not registered"):
        user_service.reset_password("nonexistent", "newpass123")


def test_update_last_login(user_service, redis_client):
    """Test updating last login timestamp."""
    import time

    from zndraw.app.redis_keys import UserKeys

    user_name = "user"
    keys = UserKeys(user_name)
    user_service.create_user(user_name)

    original_last_login = redis_client.hget(keys.hash_key(), "lastLogin")

    time.sleep(0.1)

    user_service.update_last_login(user_name)

    new_last_login = redis_client.hget(keys.hash_key(), "lastLogin")

    if isinstance(original_last_login, bytes):
        original_last_login = original_last_login.decode("utf-8")
    if isinstance(new_last_login, bytes):
        new_last_login = new_last_login.decode("utf-8")

    assert new_last_login != original_last_login


def test_delete_user(user_service, admin_service_deployment):
    """Test deleting a user."""
    user_name = "DeleteMe"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password")
    admin_service_deployment.grant_admin(user_name)

    assert user_service.username_exists(user_name)
    assert admin_service_deployment.is_admin(user_name)

    result = user_service.delete_user(user_name)
    assert result is True

    assert not user_service.username_exists(user_name)
    assert not admin_service_deployment.is_admin(user_name)


def test_delete_user_with_admin_status(user_service, admin_service_deployment):
    """Test deleting a user also removes admin status."""
    user_name = "AdminUser"

    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password")
    admin_service_deployment.grant_admin(user_name)

    assert admin_service_deployment.is_admin(user_name)

    user_service.delete_user(user_name)

    assert not admin_service_deployment.is_admin(user_name)


def test_delete_nonexistent_user(user_service):
    """Test deleting non-existent user (should not raise error)."""
    result = user_service.delete_user("nonexistent")
    assert result is True


def test_list_all_users(user_service, admin_service_deployment):
    """Test listing all users with different roles."""
    admin_service = admin_service_deployment

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


def test_list_all_users_with_visited_rooms(user_service, redis_client):
    """Test listing users when visited_rooms keys exist.

    This tests the bug where scanning for 'user:*' also matches
    'user:{username}:visited_rooms' keys, which are Sets not Hashes,
    causing WRONGTYPE errors when hgetall() is called on them.

    Note: Username must NOT contain 'user' substring to trigger the bug,
    because str.replace('user:', '') replaces all occurrences.
    """
    from zndraw.app.redis_keys import UserKeys

    # Use 'admin' - a username without 'user' substring to properly trigger bug
    user_name = "admin"
    user_service.create_user(user_name)

    # Add visited rooms (creates a Set key that matches 'user:*' pattern)
    keys = UserKeys(user_name)
    redis_client.sadd(keys.visited_rooms(), "room1")
    redis_client.sadd(keys.visited_rooms(), "room2")

    # This should NOT fail even though user:admin:visited_rooms exists
    users = user_service.list_all_users()

    assert len(users) == 1
    assert users[0]["userName"] == user_name


def test_password_hashing_produces_unique_hashes(user_service, redis_client):
    """Test that same password produces different hashes (Argon2 random salt)."""
    from zndraw.app.redis_keys import UserKeys

    user_name_1 = "user-1"
    user_name_2 = "user-2"
    same_password = "samepassword123"
    keys_1 = UserKeys(user_name_1)
    keys_2 = UserKeys(user_name_2)

    user_service.create_user(user_name_1)
    user_service.register_user(user_name_1, user_name_1, same_password)

    user_service.create_user(user_name_2)
    user_service.register_user(user_name_2, user_name_2, same_password)

    # Both users should be able to verify their password
    assert user_service.verify_password(user_name_1, same_password)
    assert user_service.verify_password(user_name_2, same_password)

    # Argon2 includes random salt in the hash, so same password produces different hashes
    hash_1 = redis_client.hget(keys_1.hash_key(), "passwordHash")
    hash_2 = redis_client.hget(keys_2.hash_key(), "passwordHash")

    if isinstance(hash_1, bytes):
        hash_1 = hash_1.decode("utf-8")
    if isinstance(hash_2, bytes):
        hash_2 = hash_2.decode("utf-8")

    # Hashes should be different due to random salt
    assert hash_1 != hash_2
    # Argon2 hashes start with $argon2
    assert hash_1.startswith("$argon2")
    assert hash_2.startswith("$argon2")


def test_password_requirements_enforced(user_service):
    """Test that password requirements are enforced."""
    from zndraw.services.user_service import PasswordValidationError

    base_name = "test-user"

    # Very short password (< 8 chars) should be rejected
    user_name_1 = f"{base_name}-1"
    user_service.create_user(user_name_1)
    with pytest.raises(PasswordValidationError, match="at least 8 characters"):
        user_service.register_user(user_name_1, user_name_1, "1")

    # Password with only spaces (< 8 chars) should be rejected
    user_name_2 = f"{base_name}-2"
    user_service.create_user(user_name_2)
    with pytest.raises(PasswordValidationError, match="at least 8 characters"):
        user_service.register_user(user_name_2, user_name_2, "   ")

    # Valid password (8+ chars) should be accepted
    user_name_3 = f"{base_name}-3"
    user_service.create_user(user_name_3)
    user_service.register_user(user_name_3, user_name_3, "validpass123")
    assert user_service.verify_password(user_name_3, "validpass123")

    # Very long password should be accepted
    user_name_4 = f"{base_name}-4"
    long_pass = "x" * 1000
    user_service.create_user(user_name_4)
    user_service.register_user(user_name_4, user_name_4, long_pass)
    assert user_service.verify_password(user_name_4, long_pass)


def test_register_user_transfers_admin_status(user_service, admin_service_deployment):
    """Test that registering with new username transfers admin status."""
    admin_service = admin_service_deployment

    guest_name = "guest-admin"
    new_username = "admin-user"

    user_service.create_user(guest_name)
    admin_service.grant_admin(guest_name)
    assert admin_service.is_admin(guest_name)

    user_service.register_user(guest_name, new_username, "password")

    assert not admin_service.is_admin(guest_name)
    assert admin_service.is_admin(new_username)


def test_ensure_user_exists_creates_new_user(user_service, redis_client):
    """Test ensure_user_exists creates user when missing."""
    from zndraw.app.redis_keys import UserKeys

    user_name = "new-user"
    keys = UserKeys(user_name)

    # User should not exist initially
    assert not user_service.username_exists(user_name)

    # ensure_user_exists should create the user
    created = user_service.ensure_user_exists(user_name)
    assert created is True

    # User should now exist
    assert user_service.username_exists(user_name)
    assert redis_client.exists(keys.hash_key())
    assert redis_client.hget(keys.hash_key(), "userName") == user_name
    assert redis_client.hexists(keys.hash_key(), "createdAt")
    assert redis_client.hexists(keys.hash_key(), "lastLogin")


def test_ensure_user_exists_idempotent(user_service):
    """Test ensure_user_exists is idempotent - doesn't create duplicate."""
    user_name = "test-user"

    # First call creates user
    created1 = user_service.ensure_user_exists(user_name)
    assert created1 is True

    # Second call should return False (already exists)
    created2 = user_service.ensure_user_exists(user_name)
    assert created2 is False

    # User still exists
    assert user_service.username_exists(user_name)


def test_ensure_user_exists_preserves_existing_data(user_service, redis_client):
    """Test ensure_user_exists doesn't modify existing user data."""
    from zndraw.app.redis_keys import UserKeys

    user_name = "existing-user"
    keys = UserKeys(user_name)

    # Create and register a user
    user_service.create_user(user_name)
    user_service.register_user(user_name, user_name, "password123")

    # Get original data
    original_hash = redis_client.hget(keys.hash_key(), "passwordHash")
    original_created_at = redis_client.hget(keys.hash_key(), "createdAt")

    # ensure_user_exists should not modify existing user
    created = user_service.ensure_user_exists(user_name)
    assert created is False

    # Data should be unchanged
    assert redis_client.hget(keys.hash_key(), "passwordHash") == original_hash
    assert redis_client.hget(keys.hash_key(), "createdAt") == original_created_at
    assert user_service.verify_password(user_name, "password123")
