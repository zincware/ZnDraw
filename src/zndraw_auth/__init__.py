"""ZnDraw Auth - Shared authentication for ZnDraw ecosystem.

Example usage:
    from zndraw_auth import (
        current_active_user,
        current_superuser,
        fastapi_users,
        auth_backend,
        get_session,
        User,
        UserRead,
        UserCreate,
    )

    # In your FastAPI app:
    app.include_router(
        fastapi_users.get_auth_router(auth_backend),
        prefix="/auth/jwt",
        tags=["auth"],
    )

    @app.get("/protected")
    async def protected(user: User = Depends(current_active_user)):
        return {"user_id": str(user.id)}
"""

from zndraw._version import __version__ as __version__
from zndraw_auth.admin import admin_token_router
from zndraw_auth.cli_login import cli_login_router
from zndraw_auth.db import (
    Base,
    CLILoginChallenge,
    SessionDep,
    User,
    create_engine_for_url,
    ensure_default_admin,
    get_engine,
    get_session,
    get_session_maker,
    get_user_db,
)
from zndraw_auth.schemas import (
    CLILoginCreateResponse,
    CLILoginStatusResponse,
    ImpersonationTokenResponse,
    TokenResponse,
    UserCreate,
    UserRead,
    UserUpdate,
)
from zndraw_auth.settings import AuthSettings, AuthSettingsDep, get_auth_settings
from zndraw_auth.users import (
    UserManager,
    auth_backend,
    current_active_user,
    current_optional_user,
    current_superuser,
    current_user_scoped_session,
    fastapi_users,
    get_user_manager,
)

__all__ = [
    "AuthSettings",
    "AuthSettingsDep",
    "Base",
    "CLILoginChallenge",
    "CLILoginCreateResponse",
    "CLILoginStatusResponse",
    "ImpersonationTokenResponse",
    "SessionDep",
    "TokenResponse",
    "User",
    "UserCreate",
    "UserManager",
    "UserRead",
    "UserUpdate",
    "admin_token_router",
    "auth_backend",
    "cli_login_router",
    "create_engine_for_url",
    "current_active_user",
    "current_optional_user",
    "current_superuser",
    "current_user_scoped_session",
    "ensure_default_admin",
    "fastapi_users",
    "get_auth_settings",
    "get_engine",
    "get_session",
    "get_session_maker",
    "get_user_db",
    "get_user_manager",
]
