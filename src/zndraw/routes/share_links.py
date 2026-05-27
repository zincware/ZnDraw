"""Room share-link management."""

import secrets
from datetime import UTC, datetime
from uuid import UUID

from fastapi import APIRouter, status
from sqlmodel import col, select

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

router = APIRouter(
    prefix="/v1/rooms/{owner}/{room_name}/share-links", tags=["share-links"]
)


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(RoomNotFound, Forbidden),
)
async def create_share_link(
    session: SessionDep,
    access: AccessManageDep,
    current_user: CurrentUserDep,
    payload: ShareLinkCreate,
) -> ShareLinkResponse:
    """Create a new share link for a room. Requires manage capability."""
    room_id = access.room.id
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
    access: AccessManageDep,
) -> CollectionResponse[ShareLinkResponse]:
    """List all non-revoked share links for a room. Requires manage capability."""
    room_id = access.room.id
    result = await session.exec(
        select(RoomShareLink).where(
            RoomShareLink.room_id == room_id,
            col(RoomShareLink.revoked_at).is_(None),
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
    access: AccessManageDep,
    link_id: UUID,
) -> None:
    """Revoke a share link (sets revoked_at, keeps the row for audit)."""
    room_id = access.room.id
    link = await session.get(RoomShareLink, link_id)
    if link is None or link.room_id != room_id or link.revoked_at is not None:
        raise ShareLinkNotFound.exception("Share link not found or already revoked")
    link.revoked_at = datetime.now(UTC)
    await session.commit()
