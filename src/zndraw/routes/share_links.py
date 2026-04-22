"""Room share-link management."""

import secrets
from datetime import UTC, datetime
from uuid import UUID

from fastapi import APIRouter, status
from sqlmodel import select

from zndraw.dependencies import (
    AccessManageDep,
    CurrentUserDep,
    SessionDep,
)
from zndraw.exceptions import (
    Forbidden,
    RoomNotFound,
    ShareLinkNotFound,
    problem_responses,
)
from zndraw.models import RoomShareLink
from zndraw.schemas import (
    CollectionResponse,
    ShareLinkCreate,
    ShareLinkResponse,
)

router = APIRouter(prefix="/v1/rooms/{room_id}/share-links", tags=["share-links"])


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(RoomNotFound, Forbidden),
)
async def create_share_link(
    session: SessionDep,
    access: AccessManageDep,  # noqa: ARG001
    current_user: CurrentUserDep,
    room_id: str,
    payload: ShareLinkCreate,
) -> ShareLinkResponse:
    """Create a new share link for a room. Requires manage capability.

    Parameters
    ----------
    session
        Async database session.
    access
        Access context verifying manage capability (side effect only).
    current_user
        Authenticated user creating the share link.
    room_id
        Path parameter identifying the room.
    payload
        Share link creation payload with access level and optional expiry.

    Returns
    -------
    ShareLinkResponse
        The created share link.
    """
    link = RoomShareLink(
        room_id=room_id,
        token=secrets.token_urlsafe(32),
        access=payload.access,
        expires_at=payload.expires_at,
        created_by_id=current_user.id,
    )
    session.add(link)
    await session.commit()
    await session.refresh(link)
    return ShareLinkResponse.model_validate(link)


@router.get(
    "",
    responses=problem_responses(RoomNotFound, Forbidden),
)
async def list_share_links(
    session: SessionDep,
    access: AccessManageDep,  # noqa: ARG001
    room_id: str,
) -> CollectionResponse[ShareLinkResponse]:
    """List all non-revoked share links for a room. Requires manage capability.

    Parameters
    ----------
    session
        Async database session.
    access
        Access context verifying manage capability (side effect only).
    room_id
        Path parameter identifying the room.

    Returns
    -------
    CollectionResponse[ShareLinkResponse]
        Collection of active share links.
    """
    result = await session.exec(
        select(RoomShareLink).where(
            RoomShareLink.room_id == room_id,
            RoomShareLink.revoked_at.is_(None),
        )
    )
    return CollectionResponse(
        items=[ShareLinkResponse.model_validate(link) for link in result.all()]
    )


@router.delete(
    "/{link_id}",
    status_code=status.HTTP_204_NO_CONTENT,
    responses=problem_responses(RoomNotFound, Forbidden, ShareLinkNotFound),
)
async def revoke_share_link(
    session: SessionDep,
    access: AccessManageDep,  # noqa: ARG001
    room_id: str,
    link_id: UUID,
) -> None:
    """Revoke a share link (sets revoked_at, keeps the row for audit).

    Parameters
    ----------
    session
        Async database session.
    access
        Access context verifying manage capability (side effect only).
    room_id
        Path parameter identifying the room.
    link_id
        UUID of the share link to revoke.
    """
    link = await session.get(RoomShareLink, link_id)
    if link is None or link.room_id != room_id or link.revoked_at is not None:
        raise ShareLinkNotFound.exception("Share link not found or already revoked")
    link.revoked_at = datetime.now(UTC)
    await session.commit()
