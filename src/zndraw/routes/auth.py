"""Authentication REST API endpoints using zndraw-auth."""

from typing import Annotated
from uuid import uuid4

from fastapi import APIRouter, Depends, HTTPException
from fastapi.security import OAuth2PasswordRequestForm
from fastapi_users.authentication import JWTStrategy

from zndraw.config import Settings, get_zndraw_settings
from zndraw_auth import (
    AuthSettingsDep,
    UserCreate,
    UserManager,
    UserRead,
    UserUpdate,
    auth_backend,
    cli_login_router,
    fastapi_users,
    get_user_manager,
)

router = APIRouter(prefix="/v1/auth", tags=["auth"])


@router.post("/guest")
async def create_guest_session(
    auth_settings: AuthSettingsDep,
    user_manager: Annotated[UserManager, Depends(get_user_manager)],
    settings: Annotated[Settings, Depends(get_zndraw_settings)],
) -> dict:
    """Create anonymous guest user and return JWT token."""
    email = f"{uuid4().hex[:8]}@guest.user"
    password = settings.guest_password.get_secret_value()

    user = await user_manager.create(UserCreate(email=email, password=password))

    strategy = JWTStrategy(
        secret=auth_settings.secret_key.get_secret_value(),
        lifetime_seconds=auth_settings.token_lifetime_seconds,
    )
    token = await strategy.write_token(user)

    return {"access_token": token, "token_type": "bearer", "email": email}


@router.post("/jwt/login")
async def login(
    settings: Annotated[Settings, Depends(get_zndraw_settings)],
    auth_settings: AuthSettingsDep,
    user_manager: Annotated[UserManager, Depends(get_user_manager)],
    credentials: Annotated[OAuth2PasswordRequestForm, Depends()],
) -> dict:
    """Login endpoint with internal worker email guard.

    Shadows the fastapi-users login route (first-match wins in FastAPI)
    to block the internal service account from public login.
    """
    if credentials.username == settings.internal_worker_email:
        raise HTTPException(
            status_code=403,
            detail="Internal service accounts cannot log in via the public API",
        )

    user = await user_manager.authenticate(credentials)
    if user is None or not user.is_active:
        raise HTTPException(status_code=400, detail="LOGIN_BAD_CREDENTIALS")

    strategy = JWTStrategy(
        secret=auth_settings.secret_key.get_secret_value(),
        lifetime_seconds=auth_settings.token_lifetime_seconds,
    )
    token = await strategy.write_token(user)
    return {"access_token": token, "token_type": "bearer"}


# Include fastapi-users routers
router.include_router(
    fastapi_users.get_auth_router(auth_backend),
    prefix="/jwt",
)
router.include_router(
    fastapi_users.get_register_router(UserRead, UserCreate),
    prefix="",
)
router.include_router(
    fastapi_users.get_users_router(UserRead, UserUpdate),
    prefix="/users",
)
router.include_router(cli_login_router, prefix="/cli-login")
