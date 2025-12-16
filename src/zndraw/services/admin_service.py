"""Admin access control service.

Manages admin user privileges and deployment mode detection.
"""

import logging

from redis import Redis

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
            log.info("Admin service initialized in DEPLOYMENT mode")
            log.info(f"Admin username: {admin_username} (password configured)")
        else:
            log.info("Admin service initialized in LOCAL mode (all users are admin)")

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

        return username == self._admin_username and password == self._admin_password

    def grant_admin(self, user_name: str) -> None:
        """Grant admin privileges to a user.

        Parameters
        ----------
        user_name : str
            Username to grant admin status
        """
        # Use hash-based storage for MemoryStorage compatibility
        key = f"admin:user:{user_name}"
        self.r.set(key, "1")
        log.info(f"Granted admin privileges to user {user_name}")

    def revoke_admin(self, user_name: str) -> None:
        """Revoke admin privileges from a user.

        Parameters
        ----------
        user_name : str
            Username to revoke admin status
        """
        key = f"admin:user:{user_name}"
        self.r.delete(key)
        log.info(f"Revoked admin privileges from user {user_name}")

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
        key = f"admin:user:{user_name}"
        result = self.r.get(key)
        return result == "1" or result == b"1"

    def get_all_admins(self) -> set[str]:
        """Get all usernames with admin privileges.

        Returns
        -------
        set[str]
            Set of admin usernames
        """
        admins = set()
        # Scan for all admin:user:* keys
        for key in self.r.scan_iter(match="admin:user:*"):
            user_name = key.replace("admin:user:", "")
            if self.r.get(key) == "1":
                admins.add(user_name)
        return admins
