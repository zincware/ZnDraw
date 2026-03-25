"""Test fixtures for zndraw-auth."""

from collections.abc import AsyncGenerator
from typing import Annotated

import pytest
from fastapi import Depends, FastAPI
from httpx import ASGITransport, AsyncClient
from pydantic import BaseModel
from sqlalchemy import text
from sqlalchemy.ext.asyncio import async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool

from zndraw_auth import (
    SessionDep,
    User,
    UserCreate,
    UserRead,
    UserUpdate,
    auth_backend,
    current_active_user,
    current_optional_user,
    current_superuser,
    current_user_scoped_session,
    ensure_default_admin,
    fastapi_users,
)
from zndraw_auth.admin import admin_token_router
from zndraw_auth.cli_login import cli_login_router
from zndraw_auth.db import Base
from zndraw_auth.settings import AuthSettings

# --- Shared Test Models ---


class LoginForm(BaseModel):
    """OAuth2 password login form data for testing.

    Uses grant_type=password for RFC 6749 OAuth2 password flow.
    """

    username: str  # email
    password: str
    grant_type: str = "password"
    scope: str = ""
    client_id: str | None = None
    client_secret: str | None = None


# --- Settings Fixtures ---


@pytest.fixture
def login_form_class() -> type:
    """Fixture that returns the LoginForm class for dependency injection."""
    return LoginForm


@pytest.fixture
def test_settings() -> AuthSettings:
    """Settings with in-memory database (production mode with admin configured)."""
    return AuthSettings(
        secret_key="test-secret-key",
        reset_password_token_secret="test-reset-secret",
        verification_token_secret="test-verify-secret",
        # Production mode: admin configured, new users are NOT superusers
        default_admin_email="admin@test.com",
        default_admin_password="admin-password",
    )


@pytest.fixture
def test_settings_dev_mode() -> AuthSettings:
    """Settings in dev mode (no admin configured, all users become superusers)."""
    return AuthSettings(
        secret_key="test-secret-key",
        reset_password_token_secret="test-reset-secret",
        verification_token_secret="test-verify-secret",
        # Dev mode: no admin configured
    )


# --- App Helper ---


async def _create_test_app(
    settings: AuthSettings, *, create_admin: bool = True
) -> FastAPI:
    """Create a FastAPI test app with auth routes and test database.

    Parameters
    ----------
    settings : AuthSettings
        Auth settings to store in app.state.
    create_admin : bool
        Whether to create the default admin user.
    """
    test_engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    test_session_maker = async_sessionmaker(test_engine, expire_on_commit=False)

    app = FastAPI()

    # Store state for DI
    app.state.engine = test_engine
    app.state.session_maker = test_session_maker
    app.state.auth_settings = settings

    # Create all tables
    async with test_engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)

    # Optionally create default admin
    if create_admin:
        async with test_session_maker() as session:
            await ensure_default_admin(session, settings)

    # Include auth routers
    app.include_router(
        fastapi_users.get_auth_router(auth_backend),
        prefix="/auth/jwt",
        tags=["auth"],
    )
    app.include_router(
        fastapi_users.get_register_router(UserRead, UserCreate),
        prefix="/auth",
        tags=["auth"],
    )
    app.include_router(
        fastapi_users.get_users_router(UserRead, UserUpdate),
        prefix="/users",
        tags=["users"],
    )
    app.include_router(cli_login_router, prefix="/auth/cli-login", tags=["auth"])
    app.include_router(admin_token_router, prefix="/admin", tags=["admin"])

    # Test routes for dependency injection
    @app.get("/test/protected")
    async def protected_route(
        user: Annotated[User, Depends(current_active_user)],
    ) -> dict[str, str]:
        """Route requiring authenticated active user."""
        return {"user_id": str(user.id), "email": user.email}

    @app.get("/test/superuser")
    async def superuser_route(
        user: Annotated[User, Depends(current_superuser)],
    ) -> dict[str, str]:
        """Route requiring superuser."""
        return {"user_id": str(user.id), "is_superuser": str(user.is_superuser)}

    @app.get("/test/optional")
    async def optional_route(
        user: Annotated[User | None, Depends(current_optional_user)],
    ) -> dict[str, str | None]:
        """Route with optional authentication."""
        if user:
            return {"user_id": str(user.id), "authenticated": "true"}
        return {"user_id": None, "authenticated": "false"}

    @app.get("/test/session")
    async def session_route(
        session: SessionDep,
    ) -> dict[str, str]:
        """Route using async session dependency."""
        result = await session.execute(text("SELECT 1"))
        value = result.scalar()
        return {"db_check": str(value)}

    @app.get("/test/scoped-session")
    async def scoped_session_route(
        user: Annotated[User, Depends(current_user_scoped_session)],
    ) -> dict[str, str]:
        """Route using scoped-session auth (session closed before return)."""
        return {"user_id": str(user.id), "email": user.email}

    return app


# --- App Fixtures ---


@pytest.fixture
async def app(test_settings: AuthSettings) -> AsyncGenerator[FastAPI, None]:
    """Create test FastAPI app with dependency overrides."""
    app = await _create_test_app(test_settings)

    yield app

    # Cleanup
    await app.state.engine.dispose()


@pytest.fixture
async def client(app: FastAPI) -> AsyncGenerator[AsyncClient, None]:
    """Async test client.

    HTTPX AsyncClient with ASGITransport automatically triggers the app's lifespan.
    """
    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://test",
    ) as client:
        yield client


@pytest.fixture
async def app_dev_mode(
    test_settings_dev_mode: AuthSettings,
) -> AsyncGenerator[FastAPI, None]:
    """Create test FastAPI app in dev mode (all users become superusers)."""
    app = await _create_test_app(test_settings_dev_mode, create_admin=False)

    yield app

    # Cleanup
    await app.state.engine.dispose()


@pytest.fixture
async def client_dev_mode(
    app_dev_mode: FastAPI,
) -> AsyncGenerator[AsyncClient, None]:
    """Async test client in dev mode (all users become superusers)."""
    async with AsyncClient(
        transport=ASGITransport(app=app_dev_mode),
        base_url="http://test",
    ) as client:
        yield client
