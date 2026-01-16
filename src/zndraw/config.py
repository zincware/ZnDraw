"""Centralized configuration management for ZnDraw.

Uses Pydantic Settings for type-safe configuration with environment variable support.
All components (CLI, server, Celery) use this module for configuration.

Environment variables:
    ZNDRAW_STORAGE__TYPE: "memory", "lmdb", or "mongodb"
    ZNDRAW_STORAGE__PATH: Path for LMDB storage (when type=lmdb)
    ZNDRAW_STORAGE__MAP_SIZE: LMDB map size in bytes (when type=lmdb)
    ZNDRAW_STORAGE__URL: MongoDB connection URL (when type=mongodb)
    ZNDRAW_STORAGE__DATABASE: MongoDB database name (when type=mongodb)
    ZNDRAW_MEDIA_PATH: Path for local media storage (screenshots, etc.)
    ZNDRAW_REDIS_URL: Redis connection URL
    ZNDRAW_SERVER_HOST: Server host
    ZNDRAW_SERVER_PORT: Server port
    ... and more (see ZnDrawConfig fields)

Example:
    >>> from zndraw.config import get_config
    >>> config = get_config()
    >>> print(config.storage)
    InMemoryStorageConfig(type='memory')
"""

import logging
import os
from typing import Annotated, Literal, Union
from urllib.parse import urlparse, urlunparse

from pydantic import BaseModel, Field, SecretStr, field_validator, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

log = logging.getLogger(__name__)


class InMemoryStorageConfig(BaseModel):
    """In-memory storage configuration (no persistence)."""

    type: Literal["memory"] = "memory"


class LMDBStorageConfig(BaseModel):
    """LMDB storage configuration (local file-based storage)."""

    type: Literal["lmdb"] = "lmdb"
    path: str = Field(
        default="./zndraw-data", description="Base directory for LMDB files"
    )
    map_size: int = Field(
        default=1_073_741_824,
        description="Maximum size per LMDB database in bytes (default: 1 GB)",
    )


class MongoDBStorageConfig(BaseModel):
    """MongoDB storage configuration (distributed storage)."""

    type: Literal["mongodb"] = "mongodb"
    url: str = Field(description="MongoDB connection URI")
    database: str = Field(default="zndraw", description="MongoDB database name")

    @field_validator("url")
    @classmethod
    def validate_url(cls, v: str) -> str:
        """Validate MongoDB URL format."""
        if not (v.startswith("mongodb://") or v.startswith("mongodb+srv://")):
            raise ValueError(
                "MongoDB URL must start with 'mongodb://' or 'mongodb+srv://'"
            )
        return v

    def get_masked_url(self) -> str:
        """Return URL with credentials masked for logging."""
        parsed = urlparse(self.url)
        if parsed.password:
            # Replace password with ***
            netloc = f"{parsed.username}:***@{parsed.hostname}"
            if parsed.port:
                netloc += f":{parsed.port}"
            masked = parsed._replace(netloc=netloc)
            return urlunparse(masked)
        return self.url


StorageConfig = Annotated[
    Union[InMemoryStorageConfig, LMDBStorageConfig, MongoDBStorageConfig],
    Field(discriminator="type"),
]


class ZnDrawConfig(BaseSettings):
    """ZnDraw configuration with Pydantic Settings.

    Loads from environment variables with ZNDRAW_ prefix.
    Nested models use __ delimiter (e.g., ZNDRAW_STORAGE__TYPE).
    """

    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_",
        env_nested_delimiter="__",
        case_sensitive=False,
    )

    # Storage configuration (discriminated union)
    storage: StorageConfig = Field(default_factory=InMemoryStorageConfig)

    # Media storage (always local filesystem, for screenshots etc.)
    media_path: str = Field(
        default="./zndraw-media",
        description="Local path for media files (screenshots, etc.)",
    )

    # Redis configuration
    redis_url: str | None = Field(
        default=None,
        description="Redis connection URL. None means in-memory mode.",
    )

    # Server configuration
    server_host: str = Field(default="0.0.0.0", description="Server bind host address")
    server_port: int = Field(
        default=5000, ge=1, le=65535, description="Server bind port"
    )
    server_url: str | None = Field(
        default=None,
        description="Full server URL for callbacks. Auto-generated if None.",
    )

    # Logging
    log_level: str = Field(default="WARNING", description="Logging level")

    # Security
    flask_secret_key: str = Field(
        default="dev-secret-key-change-in-production",
        description="Flask session secret key",
        alias="FLASK_SECRET_KEY",
    )
    admin_username: str | None = Field(default=None, description="Admin username")
    admin_password: SecretStr | None = Field(default=None, description="Admin password")

    # Upload settings
    upload_temp: str = Field(
        default="/tmp/zndraw_uploads",
        description="Temporary directory for file uploads. "
        "In production, configure a secure dedicated directory.",
    )
    max_upload_mb: int = Field(
        default=500, ge=1, description="Maximum upload size in MB"
    )

    # Feature flags
    simgen_enabled: bool = Field(default=False, description="Enable SiMGen features")
    file_browser_enabled: bool = Field(default=False, description="Enable file browser")
    file_browser_root: str = Field(
        default_factory=os.getcwd, description="Root directory for file browser"
    )
    lock_template_room: bool = Field(
        default=False, description="Auto-lock rooms created from CLI file loading"
    )

    @field_validator("log_level")
    @classmethod
    def validate_log_level(cls, v: str) -> str:
        """Validate and normalize log level."""
        valid_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        upper_v = v.upper()
        if upper_v not in valid_levels:
            log.warning(
                f"Invalid log level '{v}', using WARNING. Valid: {', '.join(valid_levels)}"
            )
            return "WARNING"
        return upper_v

    @model_validator(mode="after")
    def validate_admin_credentials(self) -> "ZnDrawConfig":
        """Validate that admin credentials are both set or both unset."""
        if (self.admin_username is None) != (self.admin_password is None):
            raise ValueError(
                "admin_username and admin_password must both be set or both be unset"
            )
        return self

    @model_validator(mode="after")
    def auto_generate_server_url(self) -> "ZnDrawConfig":
        """Auto-generate server_url if not explicitly set."""
        if self.server_url is None:
            url_host = (
                self.server_host if self.server_host != "0.0.0.0" else "localhost"
            )
            self.server_url = f"http://{url_host}:{self.server_port}"
        return self

    def log_config(self) -> None:
        """Log configuration for debugging (excludes sensitive data)."""
        log.debug("=" * 80)
        log.debug("ZnDraw Configuration:")
        log.debug(f"  Redis URL: {self.redis_url or 'None (in-memory mode)'}")

        match self.storage:
            case InMemoryStorageConfig():
                log.debug("  Storage Backend: In-Memory (no persistence)")
            case MongoDBStorageConfig():
                log.debug("  Storage Backend: MongoDB")
                log.debug(f"    URL: {self.storage.get_masked_url()}")
                log.debug(f"    Database: {self.storage.database}")
            case LMDBStorageConfig():
                log.debug("  Storage Backend: LMDB")
                log.debug(f"    Path: {self.storage.path}")
                log.debug(
                    f"    Map Size: {self.storage.map_size / 1024**3:.2f} GB per room"
                )

        log.debug(f"  Media Path: {self.media_path}")
        log.debug(f"  Server: {self.server_url}")
        log.debug(f"  Log Level: {self.log_level}")
        log.debug(f"  Admin Mode: {'Enabled' if self.admin_username else 'Disabled'}")
        log.debug(f"  Max Upload: {self.max_upload_mb}MB")
        log.debug(f"  SiMGen: {'Enabled' if self.simgen_enabled else 'Disabled'}")
        log.debug(
            f"  File Browser: {'Enabled' if self.file_browser_enabled else 'Disabled'}"
        )
        log.debug(
            f"  Lock Template Room: {'Enabled' if self.lock_template_room else 'Disabled'}"
        )
        log.debug("=" * 80)


# Global config instance (singleton pattern)
_config: ZnDrawConfig | None = None


def get_config() -> ZnDrawConfig:
    """Get or create the global configuration instance.

    Returns
    -------
    ZnDrawConfig
        Global configuration instance loaded from environment variables.
    """
    global _config
    if _config is None:
        _config = ZnDrawConfig()
    return _config


def set_config(config: ZnDrawConfig) -> None:
    """Set the global configuration instance.

    Parameters
    ----------
    config : ZnDrawConfig
        Configuration instance to use globally.
    """
    global _config
    _config = config


def reload_config() -> ZnDrawConfig:
    """Reload configuration from environment.

    Useful for testing or when environment variables change at runtime.

    Returns
    -------
    ZnDrawConfig
        Newly created configuration instance.
    """
    global _config
    _config = ZnDrawConfig()
    return _config
