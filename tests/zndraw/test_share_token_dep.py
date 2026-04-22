"""Integration test for the share-token resolver."""

from datetime import UTC, datetime, timedelta
from uuid import uuid4

import pytest
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.access import ShareAccess, ShareContext
from zndraw.dependencies import resolve_share_token
from zndraw.models import RoomShareLink


@pytest.mark.asyncio
async def test_resolve_valid_token(session: AsyncSession) -> None:
    room_id = "room-x"
    uid = uuid4()
    link = RoomShareLink(
        room_id=room_id, token="tkn-1", access=ShareAccess.EDIT, created_by_id=uid
    )
    session.add(link)
    await session.commit()

    ctx = await resolve_share_token(session, "tkn-1", room_id)
    assert ctx == ShareContext(room_id=room_id, access=ShareAccess.EDIT)


@pytest.mark.asyncio
async def test_resolve_returns_none_for_missing(session: AsyncSession) -> None:
    assert await resolve_share_token(session, None, "room") is None
    assert await resolve_share_token(session, "nope", "room") is None


@pytest.mark.asyncio
async def test_resolve_returns_none_for_wrong_room(session: AsyncSession) -> None:
    uid = uuid4()
    link = RoomShareLink(
        room_id="a", token="tkn-2", access=ShareAccess.VIEW, created_by_id=uid
    )
    session.add(link)
    await session.commit()

    assert await resolve_share_token(session, "tkn-2", "different") is None


@pytest.mark.asyncio
async def test_resolve_returns_none_for_revoked(session: AsyncSession) -> None:
    uid = uuid4()
    link = RoomShareLink(
        room_id="a",
        token="tkn-3",
        access=ShareAccess.VIEW,
        created_by_id=uid,
        revoked_at=datetime.now(UTC),
    )
    session.add(link)
    await session.commit()

    assert await resolve_share_token(session, "tkn-3", "a") is None


@pytest.mark.asyncio
async def test_resolve_returns_none_for_expired(session: AsyncSession) -> None:
    uid = uuid4()
    link = RoomShareLink(
        room_id="a",
        token="tkn-4",
        access=ShareAccess.VIEW,
        created_by_id=uid,
        expires_at=datetime.now(UTC) - timedelta(hours=1),
    )
    session.add(link)
    await session.commit()

    assert await resolve_share_token(session, "tkn-4", "a") is None
