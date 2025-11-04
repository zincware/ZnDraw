"""User management service.

Manages user accounts, roles, and authentication.
Uses userName as the primary identifier (no client_id).
"""

import hashlib
import logging
import secrets
from enum import Enum

from redis import Redis

from zndraw.app.redis_keys import UserKeys

log = logging.getLogger(__name__)


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
    Password storage uses SHA-256 hashing with salt.
    No password requirements enforced (as per spec).
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client
        log.info("UserService initialized")

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
        import datetime

        keys = UserKeys(user_name)
        current_time = datetime.datetime.utcnow().isoformat()
        self.r.hset(
            keys.hash_key(),
            mapping={
                "userName": user_name,
                "createdAt": current_time,
                "lastLogin": current_time,
            },
        )

        log.info(f"Created user {user_name}")
        return True

    def register_user(
        self, old_user_name: str, new_user_name: str, password: str
    ) -> bool:
        """Register a guest user with chosen username + password.

        This promotes guest → user by:
        1. Checking new username availability
        2. Creating new user entry with password
        3. Deleting old guest entry

        Parameters
        ----------
        old_user_name : str
            Current guest username
        new_user_name : str
            Desired username
        password : str
            Password to set

        Returns
        -------
        bool
            True if registration successful

        Raises
        ------
        ValueError
            If new username already exists or is invalid
        """
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

            # Generate salt and hash password
            salt = secrets.token_hex(16)
            password_hash = self._hash_password(password, salt)

            # Add password to existing user
            keys = UserKeys(old_user_name)
            self.r.hset(keys.hash_key(), "passwordHash", password_hash)
            self.r.hset(keys.hash_key(), "passwordSalt", salt)

            log.info(f"User {old_user_name} registered (guest → user)")
            return True

        # Different username - create new entry and delete old
        import datetime

        old_keys = UserKeys(old_user_name)
        new_keys = UserKeys(new_user_name)

        # Generate salt and hash password
        salt = secrets.token_hex(16)
        password_hash = self._hash_password(password, salt)

        # Get old user data
        old_data = self.r.hgetall(old_keys.hash_key())

        # Create new user entry
        current_time = datetime.datetime.utcnow().isoformat()

        # Decode bytes if needed
        created_at = old_data.get(b"createdAt") or old_data.get("createdAt")
        if isinstance(created_at, bytes):
            created_at = created_at.decode("utf-8")

        self.r.hset(
            new_keys.hash_key(),
            mapping={
                "userName": new_user_name,
                "passwordHash": password_hash,
                "passwordSalt": salt,
                "createdAt": created_at or current_time,
                "lastLogin": current_time,
            },
        )

        # Transfer admin status if exists
        if self.r.get(old_keys.admin_key()):
            self.r.set(new_keys.admin_key(), "1")
            self.r.delete(old_keys.admin_key())

        # Delete old user
        self.r.delete(old_keys.hash_key())

        log.info(
            f"User {old_user_name} registered as {new_user_name} (guest → user, username changed)"
        )
        return True

    def verify_password(self, user_name: str, password: str) -> bool:
        """Verify a user's password.

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
        salt = self.r.hget(keys.hash_key(), "passwordSalt")

        if not stored_hash or not salt:
            return False

        # Handle bytes from Redis
        if isinstance(stored_hash, bytes):
            stored_hash = stored_hash.decode("utf-8")
        if isinstance(salt, bytes):
            salt = salt.decode("utf-8")

        computed_hash = self._hash_password(password, salt)
        return computed_hash == stored_hash

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
            New password to set

        Returns
        -------
        bool
            True if password changed successfully

        Raises
        ------
        ValueError
            If user is not registered or old password is wrong
        """
        if not self.is_registered(user_name):
            raise ValueError("User is not registered")

        if not self.verify_password(user_name, old_password):
            raise ValueError("Current password is incorrect")

        # Generate new salt and hash
        salt = secrets.token_hex(16)
        password_hash = self._hash_password(new_password, salt)

        # Update password
        keys = UserKeys(user_name)
        self.r.hset(keys.hash_key(), "passwordHash", password_hash)
        self.r.hset(keys.hash_key(), "passwordSalt", salt)

        log.info(f"User {user_name} changed password")
        return True

    def reset_password(self, user_name: str, new_password: str) -> bool:
        """Admin function to reset a user's password without old password.

        Parameters
        ----------
        user_name : str
            Username
        new_password : str
            New password to set

        Returns
        -------
        bool
            True if password reset successfully

        Raises
        ------
        ValueError
            If user is not registered
        """
        if not self.is_registered(user_name):
            raise ValueError("User is not registered")

        # Generate new salt and hash
        salt = secrets.token_hex(16)
        password_hash = self._hash_password(new_password, salt)

        # Update password
        keys = UserKeys(user_name)
        self.r.hset(keys.hash_key(), "passwordHash", password_hash)
        self.r.hset(keys.hash_key(), "passwordSalt", salt)

        log.info(f"Admin reset password for user {user_name}")
        return True

    def update_last_login(self, user_name: str) -> None:
        """Update the last login timestamp for a user.

        Parameters
        ----------
        user_name : str
            Username
        """
        import datetime

        keys = UserKeys(user_name)
        current_time = datetime.datetime.utcnow().isoformat()
        self.r.hset(keys.hash_key(), "lastLogin", current_time)

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

        # Delete user data
        self.r.delete(keys.hash_key())

        # Delete admin status if exists
        self.r.delete(keys.admin_key())

        log.info(f"User {user_name} deleted")
        return True

    def list_all_users(self) -> list[dict]:
        """List all users with their roles.

        Returns
        -------
        list[dict]
            List of user dictionaries with userName, role
        """
        users = []

        # Scan for all user keys
        for key in self.r.scan_iter(match="user:*"):
            # Extract userName from key
            if isinstance(key, bytes):
                key = key.decode("utf-8")
            user_name = key.replace("user:", "")

            # Get user data
            keys = UserKeys(user_name)
            user_data = self.r.hgetall(keys.hash_key())
            if not user_data:
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

    def _hash_password(self, password: str, salt: str) -> str:
        """Hash password with salt using SHA-256.

        Parameters
        ----------
        password : str
            Password to hash
        salt : str
            Salt for hashing

        Returns
        -------
        str
            Hex-encoded hash
        """
        combined = f"{password}{salt}".encode("utf-8")
        return hashlib.sha256(combined).hexdigest()
