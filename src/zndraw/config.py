"""Application configuration.

Settings class configured via environment variables with the ZNDRAW_SERVER_ prefix.
Storage backend is selected by URI string (e.g. ``memory://``,
``/path/to/data.lmdb``, ``mongodb://host:port/db``).
"""

from pathlib import Path
from typing import Annotated

from fastapi import Depends, Request
from pydantic import SecretStr
from pydantic_settings import (
    BaseSettings,
    PydanticBaseSettingsSource,
    PyprojectTomlConfigSettingsSource,
    SettingsConfigDict,
)


class Settings(BaseSettings):
    """Application settings loaded from environment variables.

    All settings use the ``ZNDRAW_SERVER_`` prefix.

    The default ``@internal`` filesystem provider is gated on
    ``filebrowser_enabled`` (no DB row, no task registration, frontend icon
    hidden when disabled) and rooted at ``filebrowser_path``.
    """

    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_SERVER_",
        pyproject_toml_table_header=("tool", "zndraw", "server"),
    )

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,  # noqa: ARG003
        file_secret_settings: PydanticBaseSettingsSource,  # noqa: ARG003
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        """Return settings sources in priority order.

        Parameters
        ----------
        settings_cls : type[BaseSettings]
            The settings class.
        init_settings : PydanticBaseSettingsSource
            Init kwargs source.
        env_settings : PydanticBaseSettingsSource
            Environment variable source.
        dotenv_settings : PydanticBaseSettingsSource
            Dotenv file source.
        file_secret_settings : PydanticBaseSettingsSource
            Secret file source.

        Returns
        -------
        tuple[PydanticBaseSettingsSource, ...]
            Ordered sources: init, env, pyproject.toml.
        """
        return (
            init_settings,
            env_settings,
            PyprojectTomlConfigSettingsSource(settings_cls),
        )

    # Auth
    guest_password: SecretStr = SecretStr("zndraw")
    internal_worker_email: str = "worker@internal.user"

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
    internal_url: str | None = None  # For TaskIQ workers to reach FastAPI

    # Filesystem provider
    filebrowser_enabled: bool = True
    filebrowser_path: str = "."

    # Taskiq broker / result backend isolation (per-server namespacing)
    task_queue_name: str = "zndraw:tasks"
    result_backend_key_prefix: str = "zndraw"

    # Provider executor
    provider_executor_timeout: float = 30.0


def get_zndraw_settings(request: Request) -> Settings:
    """Retrieve ZnDraw settings from app.state."""
    return request.app.state.settings


SettingsDep = Annotated[Settings, Depends(get_zndraw_settings)]
