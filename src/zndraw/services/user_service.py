"""User management service.

Manages user accounts, roles, and authentication.
Uses userName as the primary identifier (no client_id).
"""

import logging
import re
from enum import Enum

from argon2 import PasswordHasher
from argon2.exceptions import VerifyMismatchError
from redis import Redis

from zndraw.app.redis_keys import GlobalIndexKeys, UserKeys
from zndraw.utils.time import utc_now_iso

log = logging.getLogger(__name__)

# Argon2 password hasher with secure defaults
_ph = PasswordHasher()

# Password requirements
MIN_PASSWORD_LENGTH = 8


class PasswordValidationError(ValueError):
    """Raised when password doesn't meet requirements."""

    pass


def validate_password(password: str) -> None:
    """Validate password meets minimum requirements.

    Requirements:
    - At least 8 characters

    Parameters
    ----------
    password : str
        Password to validate

    Raises
    ------
    PasswordValidationError
        If password doesn't meet requirements
    """
    if not password or len(password) < MIN_PASSWORD_LENGTH:
        raise PasswordValidationError(
            f"Password must be at least {MIN_PASSWORD_LENGTH} characters"
        )


# Username requirements
MAX_USERNAME_LENGTH = 64
# Pattern: 1 alphanumeric char + 0-63 more chars = 1-64 total (matches MAX_USERNAME_LENGTH)
USERNAME_PATTERN = re.compile(r"^[a-zA-Z0-9][a-zA-Z0-9_-]{0,63}$")


class UsernameValidationError(ValueError):
    """Raised when username doesn't meet requirements."""

    pass


def validate_username(user_name: str) -> None:
    """Validate username format.

    Requirements:
    - 1-64 characters
    - Alphanumeric, underscore, hyphen only
    - Must start with alphanumeric character

    Parameters
    ----------
    user_name : str
        Username to validate

    Raises
    ------
    UsernameValidationError
        If username doesn't meet requirements
    """
    if not user_name:
        raise UsernameValidationError("Username is required")

    if len(user_name) > MAX_USERNAME_LENGTH:
        raise UsernameValidationError(
            f"Username must be at most {MAX_USERNAME_LENGTH} characters"
        )

    if not USERNAME_PATTERN.match(user_name):
        raise UsernameValidationError(
            "Username must contain only letters, numbers, underscores, and hyphens, "
            "and must start with a letter or number"
        )


class UserRole(str, Enum):
    """User role levels."""

    GUEST = "guest"  # Anonymous user, no password
    USER = "user"  # Registered with password
    ADMIN = "admin"  # Elevated privileges


class UserService:
    """Handles user account management and authentication.

    Simplified architecture:
    - userName is the primary identifier (unique, immutable)
    - Guests get auto-generated usernames (user-xyz)
    - Registration selects permanent username + password
    - No client_id concept

    Parameters
    ----------
    redis_client : Redis
        Redis client instance for storing user data

    Notes
    -----
    Password storage uses Argon2id (winner of Password Hashing Competition).
    Passwords must be at least 8 characters (see validate_password()).
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client
        log.debug("UserService initialized")

    def username_exists(self, user_name: str) -> bool:
        """Check if a username already exists.

        Parameters
        ----------
        user_name : str
            Username to check

        Returns
        -------
        bool
            True if username exists
        """
        keys = UserKeys(user_name)
        return self.r.exists(keys.hash_key())

    def get_user_role(self, user_name: str) -> UserRole:
        """Get the role of a user.

        Parameters
        ----------
        user_name : str
            Username

        Returns
        -------
        UserRole
            Role of the user (guest, user, or admin)
        """
        keys = UserKeys(user_name)

        # Check admin first
        admin_value = self.r.get(keys.admin_key())
        if admin_value == "1" or admin_value == b"1":
            return UserRole.ADMIN

        # Check if user has password (registered)
        has_password = self.r.hexists(keys.hash_key(), "passwordHash")

        if has_password:
            return UserRole.USER

        return UserRole.GUEST

    def is_registered(self, user_name: str) -> bool:
        """Check if user has registered (set a password).

        Parameters
        ----------
        user_name : str
            Username

        Returns
        -------
        bool
            True if user has set a password
        """
        keys = UserKeys(user_name)
        return self.r.hexists(keys.hash_key(), "passwordHash")

    def create_user(self, user_name: str) -> bool:
        """Create a new user account (guest - no password yet).

        Parameters
        ----------
        user_name : str
            Username to create

        Returns
        -------
        bool
            True if created successfully

        Raises
        ------
        ValueError
            If username already exists
        """
        if self.username_exists(user_name):
            raise ValueError(f"Username '{user_name}' already exists")

        # Create user entry
        keys = UserKeys(user_name)
        current_time = utc_now_iso()

        pipe = self.r.pipeline()
        pipe.hset(
            keys.hash_key(),
            mapping={
                "userName": user_name,
                "createdAt": current_time,
                "lastLogin": current_time,
            },
        )
        # Add to global users index for O(1) listing
        pipe.sadd(GlobalIndexKeys.users_index(), user_name)
        pipe.execute()

        log.debug(f"Created user {user_name}")
        return True

    def ensure_user_exists(self, user_name: str) -> bool:
        """Ensure user exists in Redis, creating if necessary.

        Idempotent operation - safe to call multiple times.
        Used by /api/user/register to create guest users.

        Parameters
        ----------
        user_name : str
            Username to ensure exists

        Returns
        -------
        bool
            True if user was created, False if already existed
        """
        if self.username_exists(user_name):
            return False

        keys = UserKeys(user_name)
        current_time = utc_now_iso()

        pipe = self.r.pipeline()
        pipe.hset(
            keys.hash_key(),
            mapping={
                "userName": user_name,
                "createdAt": current_time,
                "lastLogin": current_time,
            },
        )
        pipe.sadd(GlobalIndexKeys.users_index(), user_name)
        pipe.execute()

        log.debug(f"Created guest user {user_name}")
        return True

    def register_user(
        self, old_user_name: str, new_user_name: str, password: str
    ) -> bool:
        """Register a guest user with chosen username + password.

        This promotes guest → user by:
        1. Validating password requirements
        2. Checking new username availability
        3. Creating new user entry with password
        4. Deleting old guest entry

        Parameters
        ----------
        old_user_name : str
            Current guest username
        new_user_name : str
            Desired username
        password : str
            Password to set (must be at least 8 characters)

        Returns
        -------
        bool
            True if registration successful

        Raises
        ------
        ValueError
            If new username already exists or is invalid
        PasswordValidationError
            If password doesn't meet requirements
        """
        # Validate password first
        validate_password(password)

        # Validate new username
        if not new_user_name or not new_user_name.strip():
            raise ValueError("Username cannot be empty")

        new_user_name = new_user_name.strip()

        # Check if new username is available
        if self.username_exists(new_user_name) and new_user_name != old_user_name:
            raise ValueError(f"Username '{new_user_name}' is already taken")

        # If same username, just add password
        if new_user_name == old_user_name:
            if self.is_registered(old_user_name):
                raise ValueError("User is already registered")

            password_hash = _ph.hash(password)

            keys = UserKeys(old_user_name)
            self.r.hset(keys.hash_key(), "passwordHash", password_hash)

            log.debug(f"User {old_user_name} registered (guest → user)")
            return True

        # Different username - create new entry and delete old
        old_keys = UserKeys(old_user_name)
        new_keys = UserKeys(new_user_name)

        password_hash = _ph.hash(password)

        old_data = self.r.hgetall(old_keys.hash_key())
        current_time = utc_now_iso()

        created_at = old_data.get(b"createdAt") or old_data.get("createdAt")
        if isinstance(created_at, bytes):
            created_at = created_at.decode("utf-8")

        pipe = self.r.pipeline()
        pipe.hset(
            new_keys.hash_key(),
            mapping={
                "userName": new_user_name,
                "passwordHash": password_hash,
                "createdAt": created_at or current_time,
                "lastLogin": current_time,
            },
        )

        pipe.sadd(GlobalIndexKeys.users_index(), new_user_name)
        pipe.srem(GlobalIndexKeys.users_index(), old_user_name)

        is_admin = self.r.get(old_keys.admin_key())
        if is_admin:
            pipe.set(new_keys.admin_key(), "1")
            pipe.delete(old_keys.admin_key())
            pipe.sadd(GlobalIndexKeys.admins_index(), new_user_name)
            pipe.srem(GlobalIndexKeys.admins_index(), old_user_name)

        pipe.delete(old_keys.hash_key())
        pipe.execute()

        log.debug(
            f"User {old_user_name} registered as {new_user_name} (guest → user, username changed)"
        )
        return True

    def verify_password(self, user_name: str, password: str) -> bool:
        """Verify a user's password using Argon2.

        Parameters
        ----------
        user_name : str
            Username
        password : str
            Password to verify

        Returns
        -------
        bool
            True if password matches
        """
        if not self.is_registered(user_name):
            return False

        keys = UserKeys(user_name)
        stored_hash = self.r.hget(keys.hash_key(), "passwordHash")

        if not stored_hash:
            return False

        if isinstance(stored_hash, bytes):
            stored_hash = stored_hash.decode("utf-8")

        try:
            _ph.verify(stored_hash, password)
            return True
        except VerifyMismatchError:
            return False

    def change_password(
        self, user_name: str, old_password: str, new_password: str
    ) -> bool:
        """Change a user's password.

        Parameters
        ----------
        user_name : str
            Username
        old_password : str
            Current password for verification
        new_password : str
            New password to set (must be at least 8 characters)

        Returns
        -------
        bool
            True if password changed successfully

        Raises
        ------
        ValueError
            If user is not registered or old password is wrong
        PasswordValidationError
            If new password doesn't meet requirements
        """
        # Validate new password first
        validate_password(new_password)

        if not self.is_registered(user_name):
            raise ValueError("User is not registered")

        if not self.verify_password(user_name, old_password):
            raise ValueError("Current password is incorrect")

        password_hash = _ph.hash(new_password)

        keys = UserKeys(user_name)
        self.r.hset(keys.hash_key(), "passwordHash", password_hash)

        log.debug(f"User {user_name} changed password")
        return True

    def reset_password(self, user_name: str, new_password: str) -> bool:
        """Admin function to reset a user's password without old password.

        Parameters
        ----------
        user_name : str
            Username
        new_password : str
            New password to set (must be at least 8 characters)

        Returns
        -------
        bool
            True if password reset successfully

        Raises
        ------
        ValueError
            If user is not registered
        PasswordValidationError
            If new password doesn't meet requirements
        """
        # Validate new password first
        validate_password(new_password)

        if not self.is_registered(user_name):
            raise ValueError("User is not registered")

        password_hash = _ph.hash(new_password)

        keys = UserKeys(user_name)
        self.r.hset(keys.hash_key(), "passwordHash", password_hash)

        log.debug(f"Admin reset password for user {user_name}")
        return True

    def update_last_login(self, user_name: str) -> None:
        """Update the last login timestamp for a user.

        Parameters
        ----------
        user_name : str
            Username
        """
        keys = UserKeys(user_name)
        self.r.hset(keys.hash_key(), "lastLogin", utc_now_iso())

    def delete_user(self, user_name: str) -> bool:
        """Delete a user (hard deletion).

        Removes all user data from Redis.

        Parameters
        ----------
        user_name : str
            Username

        Returns
        -------
        bool
            True if user deleted successfully
        """
        keys = UserKeys(user_name)

        pipe = self.r.pipeline()
        pipe.delete(keys.hash_key())
        pipe.delete(keys.admin_key())
        pipe.delete(keys.visited_rooms())
        # Remove from indices
        pipe.srem(GlobalIndexKeys.users_index(), user_name)
        pipe.srem(GlobalIndexKeys.admins_index(), user_name)
        pipe.execute()

        log.debug(f"User {user_name} deleted")
        return True

    def list_all_users(self) -> list[dict]:
        """List all users with their roles.

        Returns
        -------
        list[dict]
            List of user dictionaries with userName, role
        """
        users = []

        # Get all usernames from the global index
        all_usernames = self.r.smembers(GlobalIndexKeys.users_index())

        for user_name in all_usernames:
            if isinstance(user_name, bytes):
                user_name = user_name.decode("utf-8")

            # Get user data
            keys = UserKeys(user_name)
            user_data = self.r.hgetall(keys.hash_key())
            if not user_data:
                # User in index but data missing - clean up stale index entry
                self.r.srem(GlobalIndexKeys.users_index(), user_name)
                continue

            # Handle bytes from Redis
            user_dict = {}
            for k, v in user_data.items():
                if isinstance(k, bytes):
                    k = k.decode("utf-8")
                if isinstance(v, bytes):
                    v = v.decode("utf-8")
                user_dict[k] = v

            role = self.get_user_role(user_name)

            users.append(
                {
                    "userName": user_name,
                    "role": role.value,
                    "createdAt": user_dict.get("createdAt"),
                    "lastLogin": user_dict.get("lastLogin"),
                }
            )

        return users
