"""Smoke tests for the consolidated room_lookup helpers."""

import pytest

from zndraw_joblib.room_lookup import fetch_room, room_address_for


@pytest.mark.asyncio
async def test_sigils_return_none_and_self(async_session_factory) -> None:
    async with async_session_factory() as session:
        assert await fetch_room(session, "@global") is None
        assert await fetch_room(session, "@internal") is None
        assert await room_address_for(session, "@global") == "@global"
        assert await room_address_for(session, "@internal") == "@internal"


@pytest.mark.asyncio
async def test_unknown_uuid_returns_none_and_self(async_session_factory) -> None:
    unknown = "99999999-9999-9999-9999-999999999999"
    async with async_session_factory() as session:
        assert await fetch_room(session, unknown) is None
        assert await room_address_for(session, unknown) == unknown


@pytest.mark.asyncio
async def test_garbage_string_returns_none_and_self(async_session_factory) -> None:
    async with async_session_factory() as session:
        assert await fetch_room(session, "not-a-uuid") is None
        assert await room_address_for(session, "not-a-uuid") == "not-a-uuid"
