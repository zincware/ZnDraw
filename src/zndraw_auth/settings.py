"""Configuration settings for zndraw-auth."""

from typing import Annotated

from fastapi import Depends, Request
from pydantic import SecretStr, computed_field
from pydantic_settings import BaseSettings, SettingsConfigDict


class AuthSettings(BaseSettings):
    """Authentication settings loaded from environment variables.

    All settings can be overridden with ZNDRAW_AUTH_ prefix.
    Example: ZNDRAW_AUTH_SECRET_KEY=your-secret-key

    Admin mode vs Dev mode:
    - If DEFAULT_ADMIN_EMAIL and DEFAULT_ADMIN_PASSWORD are set, the system
      runs in "production mode": only the configured admin is a superuser,
      and new users are created as regular users.
    - If they are NOT set, the system runs in "dev mode": all newly
      registered users are automatically granted superuser privileges.
    """

    model_config = SettingsConfigDict(
        env_prefix="ZNDRAW_AUTH_",
        env_file=".env",
        extra="ignore",
    )

    # JWT settings
    secret_key: SecretStr = SecretStr("CHANGE-ME-IN-PRODUCTION-SECRET!")
    token_lifetime_seconds: int = 3600  # 1 hour

    # Password reset / verification tokens
    reset_password_token_secret: SecretStr = SecretStr("CHANGE-ME-RESET")
    verification_token_secret: SecretStr = SecretStr("CHANGE-ME-VERIFY")

    # Default admin user (production mode)
    default_admin_email: str | None = None
    """Email for the default admin user. If set, enables production mode."""

    default_admin_password: SecretStr | None = None
    """Password for the default admin user."""

    @computed_field  # type: ignore[prop-decorator]
    @property
    def is_dev_mode(self) -> bool:
        """True if no admin credentials configured (all users become superusers)."""
        return self.default_admin_email is None


def get_auth_settings(request: Request) -> AuthSettings:
    """Retrieve auth settings from app.state."""
    return request.app.state.auth_settings


AuthSettingsDep = Annotated[AuthSettings, Depends(get_auth_settings)]
