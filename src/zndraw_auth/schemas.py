"""Pydantic schemas for user operations."""

import uuid

from fastapi_users import schemas
from pydantic import BaseModel


class UserRead(schemas.BaseUser[uuid.UUID]):
    """Schema for reading user data (responses).

    Overrides ``email`` to plain ``str`` so that addresses with reserved
    TLDs (``.local``, ``.test``) stored by ``ensure_default_admin`` can
    be serialized without a pydantic ``EmailStr`` validation error.
    """

    email: str  # type: ignore[assignment]


class UserCreate(schemas.BaseUserCreate):
    """Schema for creating a new user."""


class UserUpdate(schemas.BaseUserUpdate):
    """Schema for updating an existing user."""


# --- OAuth2 Schemas ---


class TokenResponse(BaseModel):
    """OAuth2 bearer token response.

    This is the standard OAuth2 token response format returned by
    the /auth/jwt/login endpoint. Useful for type-safe testing and
    client implementations.

    Note: FastAPI and fastapi-users don't export a standard schema for this,
    so we provide it here for convenience.
    """

    access_token: str
    token_type: str


class CLILoginCreateResponse(BaseModel):
    """Response from creating a CLI login challenge."""

    code: str
    secret: str
    approve_url: str


class CLILoginStatusResponse(BaseModel):
    """Response from polling a CLI login challenge."""

    status: str
    token: str | None = None


class ImpersonationTokenResponse(BaseModel):
    """Response from admin token minting."""

    access_token: str
    token_type: str = "bearer"
