"""Authentication REST API endpoints using zndraw-auth."""

from typing import Annotated, Literal
from uuid import uuid4

from fastapi import APIRouter, Depends
from fastapi_users.authentication import JWTStrategy
from pydantic import BaseModel

from zndraw.config import Settings, get_zndraw_settings
from zndraw_auth import (
    AuthSettingsDep,
    SessionDep,
    UserCreate,
    UserManager,
    UserRead,
    UserUpdate,
    auth_backend,
    cli_login_router,
    fastapi_users,
    generate_unique_display_name,
    get_user_manager,
)

router = APIRouter(prefix="/v1/auth", tags=["auth"])


class GuestSessionResponse(BaseModel):
    """Response shape for ``POST /v1/auth/guest`` — a fresh anonymous session.

    Carries the freshly issued JWT plus the user's auto-generated identity
    (``email`` for backwards-readable logs; ``display_name`` for the URL/UI).
    """

    access_token: str
    token_type: Literal["bearer"] = "bearer"
    email: str
    display_name: str


@router.post("/guest")
async def create_guest_session(
    auth_settings: AuthSettingsDep,
    user_manager: Annotated[UserManager, Depends(get_user_manager)],
    settings: Annotated[Settings, Depends(get_zndraw_settings)],
    session: SessionDep,
) -> GuestSessionResponse:
    """Create anonymous guest user (is_guest=True) and return JWT token."""
    email = f"{uuid4().hex[:8]}@guest.user"
    password = settings.guest_password.get_secret_value()
    display_name = await generate_unique_display_name(session)

    user = await user_manager.create(
        UserCreate(
            email=email,
            password=password,
            is_guest=True,
            display_name=display_name,
        )
    )

    strategy = JWTStrategy(
        secret=auth_settings.secret_key.get_secret_value(),
        lifetime_seconds=auth_settings.token_lifetime_seconds,
    )
    token = await strategy.write_token(user)

    return GuestSessionResponse(
        access_token=token, email=email, display_name=user.display_name
    )


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
