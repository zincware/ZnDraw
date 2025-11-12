"""Centralized configuration management for ZnDraw.

Reads from environment variables with sensible defaults.
All components (CLI, server, Celery) use this module for configuration.

Environment variables follow the pattern ZNDRAW_* for application settings
and FLASK_* for Flask-specific settings.

Example:
    >>> from zndraw.config import get_config
    >>> config = get_config()
    >>> print(config.redis_url)
    redis://localhost:6379
"""

import logging
import os
from dataclasses import dataclass, field

log = logging.getLogger(__name__)


def _parse_bool(value: str) -> bool:
    """Parse boolean from environment variable string."""
    return value.lower() in ("true", "1", "yes", "on")


def _getenv_int(key: str, default: int) -> int:
    """Get integer from environment with fallback to default."""
    value = os.getenv(key)
    if value is None:
        return default
    try:
        return int(value)
    except ValueError:
        log.warning(f"Invalid integer value for {key}={value}, using default {default}")
        return default


@dataclass
class ZnDrawConfig:
    """ZnDraw configuration loaded from environment variables.

    All fields have defaults that work for local development.
    Production deployments should override via environment variables.

    Attributes
    ----------
    redis_url : str | None
        Redis connection URL. None means in-memory mode (single process only).
    storage_path : str
        Base directory for LMDB storage files (one .lmdb per room).
    lmdb_map_size : int
        Maximum size per LMDB database in bytes (virtual allocation).
        Default: 1 GB. Increase for rooms with large trajectories.
    server_host : str
        Server bind host address.
    server_port : int
        Server bind port number.
    server_url : str | None
        Full server URL for callbacks. Auto-generated from host/port if None.
    log_level : str
        Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    flask_secret_key : str
        Flask session secret key. MUST change in production!
    admin_username : str | None
        Admin username for protected endpoints. None disables admin mode.
    admin_password : str | None
        Admin password for protected endpoints. Must be set with username.
    upload_temp : str
        Temporary directory for file uploads.
    max_upload_mb : int
        Maximum upload size in megabytes.
    simgen_enabled : bool
        Enable SiMGen molecular generation features.
    file_browser_enabled : bool
        Enable file browser feature.
    file_browser_root : str
        Root directory for file browser.
    celery_enabled : bool
        Enable Celery background task processing.
    """

    # Core server configuration
    redis_url: str | None = field(
        default_factory=lambda: os.getenv("ZNDRAW_REDIS_URL")
    )
    storage_path: str = field(
        default_factory=lambda: os.getenv("ZNDRAW_STORAGE_PATH", "./zndraw-data")
    )
    server_host: str = field(
        default_factory=lambda: os.getenv("ZNDRAW_SERVER_HOST", "localhost")
    )
    server_port: int = field(
        default_factory=lambda: _getenv_int("ZNDRAW_SERVER_PORT", 5000)
    )
    server_url: str | None = field(
        default_factory=lambda: os.getenv("ZNDRAW_SERVER_URL")
    )

    # Logging
    log_level: str = field(
        default_factory=lambda: os.getenv("ZNDRAW_LOG_LEVEL", "WARNING")
    )

    # Security
    flask_secret_key: str = field(
        default_factory=lambda: os.getenv(
            "FLASK_SECRET_KEY", "dev-secret-key-change-in-production"
        )
    )
    admin_username: str | None = field(
        default_factory=lambda: os.getenv("ZNDRAW_ADMIN_USERNAME")
    )
    admin_password: str | None = field(
        default_factory=lambda: os.getenv("ZNDRAW_ADMIN_PASSWORD")
    )

    # Upload & Storage
    upload_temp: str = field(
        default_factory=lambda: os.getenv("ZNDRAW_UPLOAD_TEMP", "/tmp/zndraw_uploads")
    )
    max_upload_mb: int = field(
        default_factory=lambda: _getenv_int("ZNDRAW_MAX_UPLOAD_MB", 500)
    )
    lmdb_map_size: int = field(
        default_factory=lambda: _getenv_int("ZNDRAW_LMDB_MAP_SIZE", 1_073_741_824)
    )

    # Optional features
    simgen_enabled: bool = field(
        default_factory=lambda: _parse_bool(
            os.getenv("ZNDRAW_SIMGEN_ENABLED", "false")
        )
    )
    file_browser_enabled: bool = field(
        default_factory=lambda: _parse_bool(
            os.getenv("ZNDRAW_FILE_BROWSER_ENABLED", "false")
        )
    )
    file_browser_root: str = field(
        default_factory=lambda: os.getenv("ZNDRAW_FILE_BROWSER_ROOT", os.getcwd())
    )
    celery_enabled: bool = field(
        default_factory=lambda: _parse_bool(os.getenv("ZNDRAW_CELERY_ENABLED", "true"))
    )

    def __post_init__(self):
        """Validate configuration after initialization."""
        self._validate()
        self._log_config()

    def _validate(self):
        """Validate configuration values.

        Raises
        ------
        ValueError
            If configuration is invalid.
        """
        # Validate admin credentials (both must be set or both unset)
        if (self.admin_username is None) != (self.admin_password is None):
            raise ValueError(
                "ZNDRAW_ADMIN_USERNAME and ZNDRAW_ADMIN_PASSWORD must both be set or both be unset"
            )

        # Validate port range
        if not 1 <= self.server_port <= 65535:
            raise ValueError(
                f"Invalid port number: {self.server_port}. Must be between 1 and 65535"
            )

        # Validate upload size
        if self.max_upload_mb < 1:
            raise ValueError(
                f"Invalid max upload size: {self.max_upload_mb}MB. Must be at least 1MB"
            )

        # Validate log level
        valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        if self.log_level.upper() not in valid_levels:
            log.warning(
                f"Invalid log level '{self.log_level}', using WARNING. "
                f"Valid levels: {', '.join(valid_levels)}"
            )
            self.log_level = "WARNING"

        # Auto-generate server_url if not explicitly set via environment variable
        # Check if server_url came from environment or was auto-generated
        env_server_url = os.getenv("ZNDRAW_SERVER_URL")
        if env_server_url is None:
            # Not set via environment, so regenerate from host:port
            # Don't use "localhost" in URL if host is 0.0.0.0 (Docker case)
            url_host = self.server_host
            if url_host == "0.0.0.0":
                url_host = "localhost"
            self.server_url = f"http://{url_host}:{self.server_port}"

    def _log_config(self):
        """Log configuration for debugging (excludes sensitive data)."""
        log.info("=" * 80)
        log.info("ZnDraw Configuration:")
        log.info(f"  Redis URL: {self.redis_url or 'None (in-memory mode)'}")
        log.info(f"  Storage Path: {self.storage_path}")
        log.info(f"  LMDB Map Size: {self.lmdb_map_size / 1024**3:.2f} GB per room")
        log.info(f"  Server: {self.server_url}")
        log.info(f"  Log Level: {self.log_level}")
        log.info(
            f"  Admin Mode: {'Enabled' if self.admin_username else 'Disabled'}"
        )
        log.info(f"  Max Upload: {self.max_upload_mb}MB")
        log.info(f"  SiMGen: {'Enabled' if self.simgen_enabled else 'Disabled'}")
        log.info(
            f"  File Browser: {'Enabled' if self.file_browser_enabled else 'Disabled'}"
        )
        log.info(f"  Celery: {'Enabled' if self.celery_enabled else 'Disabled'}")
        log.info("=" * 80)


# Global config instance (singleton pattern)
_config: ZnDrawConfig | None = None


def get_config() -> ZnDrawConfig:
    """Get or create the global configuration instance.

    Returns
    -------
    ZnDrawConfig
        Global configuration instance loaded from environment variables.

    Example
    -------
    >>> from zndraw.config import get_config
    >>> config = get_config()
    >>> print(config.storage_path)
    ./zndraw-data
    """
    global _config
    if _config is None:
        _config = ZnDrawConfig()
    return _config


def reload_config() -> ZnDrawConfig:
    """Reload configuration from environment.

    Useful for testing or when environment variables change at runtime.

    Returns
    -------
    ZnDrawConfig
        Newly created configuration instance.

    Example
    -------
    >>> import os
    >>> os.environ["ZNDRAW_STORAGE_PATH"] = "/new/path"
    >>> config = reload_config()
    >>> print(config.storage_path)
    /new/path
    """
    global _config
    _config = ZnDrawConfig()
    return _config
