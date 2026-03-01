"""Application configuration.

Settings class configured via environment variables with the ZNDRAW_ prefix.
Storage backend is selected by URI string (e.g. ``memory://``,
``/path/to/data.lmdb``, ``mongodb://host:port/db``).
"""

from pathlib import Path
from typing import Annotated

from fastapi import Depends, Request
from pydantic import SecretStr
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    """Application settings loaded from environment variables.

    All settings use the ZNDRAW_ prefix.
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
    edit_lock_ttl: int = 10  # seconds — Redis TTL for edit locks

    # Storage backend URI (memory://, *.lmdb, mongodb://host/db)
    storage: str = "memory://"

    # Server configuration
    media_path: Path = Path("zndraw-media")
    host: str = "0.0.0.0"
    port: int = 8000

    # Feature flags
    simgen_enabled: bool = False

    # Worker
    worker_enabled: bool = True  # False in Docker (dedicated workers)
    server_url: str | None = None  # For TaskIQ workers to reach FastAPI


def get_zndraw_settings(request: Request) -> Settings:
    """Retrieve ZnDraw settings from app.state."""
    return request.app.state.settings


SettingsDep = Annotated[Settings, Depends(get_zndraw_settings)]
