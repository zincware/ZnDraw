"""Admin access control service.

Manages admin user privileges and deployment mode detection.
"""

import logging
import secrets

from redis import Redis

from zndraw.app.redis_keys import GlobalIndexKeys, UserKeys

log = logging.getLogger(__name__)


class AdminService:
    """Handles admin access control and deployment mode.

    Deployment mode is automatically enabled when both admin_username
    and admin_password are configured in the ZnDrawConfig.

    In local mode (no admin credentials), all users are admins.
    In deployment mode (admin credentials set), only the admin account has privileges.

    Parameters
    ----------
    redis_client : Redis
        Redis client instance for storing admin status
    admin_username : str | None
        Admin username from config (None for local mode)
    admin_password : str | None
        Admin password from config (None for local mode)
    """

    def __init__(
        self,
        redis_client: Redis,
        admin_username: str | None = None,
        admin_password: str | None = None,
    ):
        self.r = redis_client
        self._admin_username = admin_username
        self._admin_password = admin_password
        self._deployment_mode = bool(admin_username and admin_password)

        if self._deployment_mode:
            log.debug("Admin service initialized in DEPLOYMENT mode")
            log.debug(f"Admin username: {admin_username} (password configured)")
        else:
            log.debug("Admin service initialized in LOCAL mode (all users are admin)")

    def is_deployment_mode(self) -> bool:
        """Check if running in deployment mode.

        Returns
        -------
        bool
            True if deployment mode is active
        """
        return self._deployment_mode

    def validate_admin_credentials(self, username: str, password: str) -> bool:
        """Validate admin credentials against configured values.

        Parameters
        ----------
        username : str
            Username to validate
        password : str
            Password to validate

        Returns
        -------
        bool
            True if credentials match admin configuration
        """
        if not self._deployment_mode:
            return False

        # Use constant-time comparison to prevent timing attacks
        username_match = secrets.compare_digest(
            username.encode("utf-8"),
            (self._admin_username or "").encode("utf-8"),
        )
        password_match = secrets.compare_digest(
            password.encode("utf-8"),
            (self._admin_password or "").encode("utf-8"),
        )
        return username_match and password_match

    def is_admin_username(self, username: str) -> bool:
        """Check if username matches the configured admin username.

        Used to detect wrong password attempts for admin user.

        Parameters
        ----------
        username : str
            Username to check

        Returns
        -------
        bool
            True if username matches configured admin username
        """
        return self._deployment_mode and username == self._admin_username

    def grant_admin(self, user_name: str) -> None:
        """Grant admin privileges to a user.

        Parameters
        ----------
        user_name : str
            Username to grant admin status
        """
        keys = UserKeys(user_name)
        pipe = self.r.pipeline()
        pipe.set(keys.admin_key(), "1")
        pipe.sadd(GlobalIndexKeys.admins_index(), user_name)
        pipe.execute()
        log.debug(f"Granted admin privileges to user {user_name}")

    def revoke_admin(self, user_name: str) -> None:
        """Revoke admin privileges from a user.

        Parameters
        ----------
        user_name : str
            Username to revoke admin status
        """
        keys = UserKeys(user_name)
        pipe = self.r.pipeline()
        pipe.delete(keys.admin_key())
        pipe.srem(GlobalIndexKeys.admins_index(), user_name)
        pipe.execute()
        log.debug(f"Revoked admin privileges from user {user_name}")

    def is_admin(self, user_name: str) -> bool:
        """Check if a user has admin privileges.

        In local mode, always returns True.
        In deployment mode, checks Redis for admin status.

        Parameters
        ----------
        user_name : str
            Username to check

        Returns
        -------
        bool
            True if user has admin privileges
        """
        # In local mode, everyone is admin
        if not self._deployment_mode:
            return True

        # In deployment mode, check Redis
        keys = UserKeys(user_name)
        result = self.r.get(keys.admin_key())
        return result == "1" or result == b"1"

    def get_all_admins(self) -> set[str]:
        """Get all usernames with admin privileges.

        Returns
        -------
        set[str]
            Set of admin usernames
        """
        # Get all admin usernames from the global index
        return self.r.smembers(GlobalIndexKeys.admins_index())
