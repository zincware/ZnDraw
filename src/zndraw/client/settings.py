"""Client-side settings resolved via pydantic-settings source chain.

Sources (highest to lowest priority):
    init args > env vars (ZNDRAW_*) > pyproject.toml [tool.zndraw] > StateFileSource
"""

from __future__ import annotations

from pydantic import (
    SecretStr,  # noqa: TC002 — required at runtime for pydantic validation
)
from pydantic_settings import (
    BaseSettings,
    PydanticBaseSettingsSource,
    PyprojectTomlConfigSettingsSource,
    SettingsConfigDict,
)

from zndraw.settings_sources import StateFileSource


class ClientSettings(BaseSettings):
    """Connection settings for the ZnDraw client.

    Attributes
    ----------
    url : str | None
        Server URL (e.g. http://localhost:8000).
    room : str | None
        Room ID to connect to.
    user : str | None
        User email for authentication.
    password : SecretStr | None
        Password for login.
    token : str | None
        Authentication token (JWT or local_token).
    """

    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_",
        pyproject_toml_table_header=("tool", "zndraw"),
    )

    url: str | None = None
    room: str | None = None
    user: str | None = None
    password: SecretStr | None = None
    token: str | None = None

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,  # noqa: ARG003
        file_secret_settings: PydanticBaseSettingsSource,  # noqa: ARG003
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        """Return sources: init > env > pyproject.toml > state file."""
        return (
            init_settings,
            env_settings,
            PyprojectTomlConfigSettingsSource(settings_cls),
            StateFileSource(settings_cls),
        )
