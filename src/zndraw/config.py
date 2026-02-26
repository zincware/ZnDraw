"""Application configuration with pluggable storage backends.

This module provides a Settings class that can be configured via environment
variables with the ZNDRAW_ prefix. Storage backends use discriminated unions
to support multiple storage types (memory, LMDB, MongoDB).

Examples:
    Default (memory storage):
        No additional env vars needed.

    LMDB storage:
        ZNDRAW_STORAGE__TYPE=lmdb
        ZNDRAW_STORAGE__PATH=/path/to/lmdb
        ZNDRAW_STORAGE__MAP_SIZE=1073741824

    MongoDB storage:
        ZNDRAW_STORAGE__TYPE=mongodb
        ZNDRAW_STORAGE__URL=mongodb://localhost:27017
        ZNDRAW_STORAGE__DATABASE=zndraw
"""

from pathlib import Path
from typing import Annotated, Literal

from fastapi import Depends, Request
from pydantic import BaseModel, Field, SecretStr
from pydantic_settings import BaseSettings, SettingsConfigDict

# =============================================================================
# Storage Configuration Models (Discriminated Union)
# =============================================================================


class MemoryStorage(BaseModel):
    """In-memory storage configuration (default).

    Suitable for development and testing. Data is not persisted.
    """

    type: Literal["memory"] = "memory"


class LMDBStorage(BaseModel):
    """LMDB storage configuration.

    Lightning Memory-Mapped Database for fast, persistent key-value storage.

    Attributes:
        path: Directory path for LMDB database files.
        map_size: Maximum database size in bytes. Default is 1GB.
    """

    type: Literal["lmdb"] = "lmdb"
    path: Path
    map_size: int = 1024 * 1024 * 1024  # 1GB default


class MongoDBStorage(BaseModel):
    """MongoDB storage configuration.

    Suitable for distributed deployments and large datasets.

    Attributes:
        url: MongoDB connection URL.
        database: Database name to use. Default is 'zndraw'.
    """

    type: Literal["mongodb"] = "mongodb"
    url: str
    database: str = "zndraw"


# Type alias for the discriminated union
StorageConfig = Annotated[
    MemoryStorage | LMDBStorage | MongoDBStorage,
    Field(discriminator="type"),
]


# =============================================================================
# Main Settings Class
# =============================================================================


class Settings(BaseSettings):
    """Application settings loaded from environment variables.

    All settings use the ZNDRAW_ prefix. Nested settings (like storage)
    use double underscores as delimiters.

    Attributes:
        guest_password: Password used for anonymous guest accounts.
        redis_url: Optional Redis URL for pub/sub and presence.
        presence_ttl: Presence key TTL in seconds.
        storage: Frame storage configuration (memory, lmdb, or mongodb).
        media_path: Directory for storing media files.
        host: Server bind host.
        port: Server bind port.
    """

    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_",
        env_nested_delimiter="__",
    )

    # Auth
    guest_password: SecretStr = SecretStr("zndraw")
    worker_password: SecretStr = SecretStr("zndraw-worker")

    # Database
    database_url: str = "sqlite+aiosqlite://"  # In-memory SQLite by default
    init_db_on_startup: bool = True  # False for multi-worker production

    # Core settings
    redis_url: str | None = None
    presence_ttl: int = 60
    edit_lock_ttl: int = 10  # seconds — Redis TTL for edit locks

    # Storage configuration
    storage: StorageConfig = Field(default_factory=MemoryStorage)

    # Server configuration
    media_path: Path = Path("zndraw-media")
    host: str = "0.0.0.0"
    port: int = 5000


def get_zndraw_settings(request: Request) -> Settings:
    """Retrieve ZnDraw settings from app.state."""
    return request.app.state.settings


SettingsDep = Annotated[Settings, Depends(get_zndraw_settings)]
