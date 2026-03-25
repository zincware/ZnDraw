"""CLI login (device-code flow) router."""

import logging
import secrets
import string
import uuid
from datetime import UTC, datetime, timedelta
from typing import Annotated

import jwt as pyjwt
from fastapi import APIRouter, Depends, HTTPException, Query, Response
from sqlalchemy import select

from zndraw_auth.db import CLILoginChallenge, SessionDep, User
from zndraw_auth.schemas import CLILoginCreateResponse, CLILoginStatusResponse
from zndraw_auth.settings import AuthSettings, get_auth_settings
from zndraw_auth.users import current_active_user

log = logging.getLogger(__name__)

cli_login_router = APIRouter()

CHALLENGE_LIFETIME_SECONDS = 300  # 5 minutes
CODE_LENGTH = 8
CODE_ALPHABET = string.ascii_uppercase + string.digits


def _generate_code() -> str:
    return "".join(secrets.choice(CODE_ALPHABET) for _ in range(CODE_LENGTH))


@cli_login_router.post("", response_model=CLILoginCreateResponse)
async def create_cli_login_challenge(
    session: SessionDep,
) -> CLILoginCreateResponse:
    """Create a CLI login challenge (no auth required)."""
    now = datetime.now(UTC)
    code = _generate_code()
    secret = secrets.token_urlsafe(32)

    challenge = CLILoginChallenge(
        code=code,
        secret=secret,
        status="pending",
        created_at=now,
        expires_at=now + timedelta(seconds=CHALLENGE_LIFETIME_SECONDS),
    )
    session.add(challenge)
    await session.commit()

    return CLILoginCreateResponse(
        code=code,
        secret=secret,
        approve_url=f"/auth/cli-login/{code}",
    )


@cli_login_router.get("/{code}", response_model=CLILoginStatusResponse)
async def poll_cli_login_challenge(
    code: str,
    secret: str = Query(...),
    *,
    session: SessionDep,
) -> CLILoginStatusResponse:
    """Poll a CLI login challenge status."""
    now = datetime.now(UTC).replace(tzinfo=None)

    result = await session.execute(
        select(CLILoginChallenge).where(CLILoginChallenge.code == code)
    )
    challenge = result.scalar_one_or_none()

    if challenge is None or challenge.secret != secret:
        raise HTTPException(status_code=404, detail="Challenge not found")

    if now > challenge.expires_at:
        raise HTTPException(status_code=410, detail="Challenge expired")

    if challenge.status == "redeemed":
        raise HTTPException(status_code=404, detail="Challenge not found")

    if challenge.status == "approved" and challenge.token is not None:
        token = challenge.token
        challenge.token = None
        challenge.secret = None
        challenge.status = "redeemed"
        await session.commit()
        return CLILoginStatusResponse(status="approved", token=token)

    return CLILoginStatusResponse(status="pending")


def _mint_jwt(
    user_id: uuid.UUID,
    secret: str,
    lifetime_seconds: int,
    *,
    extra_claims: dict | None = None,
) -> str:
    """Mint a JWT compatible with fastapi-users JWTStrategy."""
    now = datetime.now(UTC)
    payload: dict = {
        "sub": str(user_id),
        "aud": "fastapi-users:auth",
        "iat": now,
        "exp": now + timedelta(seconds=lifetime_seconds),
    }
    if extra_claims:
        payload.update(extra_claims)
    return pyjwt.encode(payload, secret, algorithm="HS256")


@cli_login_router.patch("/{code}", status_code=200)
async def approve_cli_login_challenge(
    code: str,
    user: Annotated[User, Depends(current_active_user)],
    settings: Annotated[AuthSettings, Depends(get_auth_settings)],
    session: SessionDep,
) -> dict[str, str]:
    """Approve a CLI login challenge (browser user, auth required)."""
    result = await session.execute(
        select(CLILoginChallenge).where(CLILoginChallenge.code == code)
    )
    challenge = result.scalar_one_or_none()

    if challenge is None or challenge.status != "pending":
        raise HTTPException(status_code=404, detail="Challenge not found")

    now = datetime.now(UTC).replace(tzinfo=None)
    if now > challenge.expires_at:
        raise HTTPException(status_code=410, detail="Challenge expired")

    token = _mint_jwt(
        user_id=user.id,
        secret=settings.secret_key.get_secret_value(),
        lifetime_seconds=settings.token_lifetime_seconds,
    )

    challenge.token = token
    challenge.user_id = user.id
    challenge.status = "approved"
    await session.commit()

    log.info("CLI login approved: user %s, code %s", user.id, code)
    return {"status": "approved"}


@cli_login_router.delete("/{code}", status_code=204)
async def reject_cli_login_challenge(
    code: str,
    user: Annotated[User, Depends(current_active_user)],
    session: SessionDep,
) -> Response:
    """Reject a CLI login challenge (browser user, auth required)."""
    result = await session.execute(
        select(CLILoginChallenge).where(CLILoginChallenge.code == code)
    )
    challenge = result.scalar_one_or_none()

    if challenge is None:
        raise HTTPException(status_code=404, detail="Challenge not found")

    await session.delete(challenge)
    await session.commit()

    log.info("CLI login rejected: by user %s, code %s", user.id, code)
    return Response(status_code=204)
