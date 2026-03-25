"""Admin token minting router."""

import logging
import uuid
from typing import Annotated

from fastapi import APIRouter, Depends, HTTPException

from zndraw_auth.cli_login import _mint_jwt
from zndraw_auth.db import SessionDep, User
from zndraw_auth.schemas import ImpersonationTokenResponse
from zndraw_auth.settings import AuthSettings, get_auth_settings
from zndraw_auth.users import current_superuser

log = logging.getLogger(__name__)

admin_token_router = APIRouter()


@admin_token_router.post(
    "/users/{user_id}/token",
    response_model=ImpersonationTokenResponse,
)
async def mint_token_for_user(
    user_id: uuid.UUID,
    admin: Annotated[User, Depends(current_superuser)],
    settings: Annotated[AuthSettings, Depends(get_auth_settings)],
    session: SessionDep,
) -> ImpersonationTokenResponse:
    """Mint a JWT for the target user (superuser only)."""
    target = await session.get(User, user_id)

    if target is None or not target.is_active:
        raise HTTPException(status_code=404, detail="User not found")

    token = _mint_jwt(
        user_id=target.id,
        secret=settings.secret_key.get_secret_value(),
        lifetime_seconds=settings.token_lifetime_seconds,
        extra_claims={"impersonated_by": str(admin.id)},
    )

    log.info("Admin %s minted token for user %s", admin.id, target.id)

    return ImpersonationTokenResponse(access_token=token)
