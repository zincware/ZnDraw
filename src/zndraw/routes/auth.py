"""Authentication REST API endpoints using zndraw-auth."""

from uuid import uuid4

from fastapi import APIRouter, Depends
from fastapi_users.authentication import JWTStrategy
from zndraw_auth import (
    AuthSettingsDep,
    UserCreate,
    UserManager,
    UserRead,
    UserUpdate,
    auth_backend,
    fastapi_users,
    get_user_manager,
)

from zndraw.config import Settings, get_zndraw_settings

router = APIRouter(prefix="/v1/auth", tags=["auth"])


@router.post("/guest")
async def create_guest_session(
    auth_settings: AuthSettingsDep,
    user_manager: UserManager = Depends(get_user_manager),
    settings: Settings = Depends(get_zndraw_settings),
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
