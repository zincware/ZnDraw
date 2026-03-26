"""Tests for Chat REST API endpoints."""

from datetime import UTC, datetime
from typing import Any

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.models import Message, Room

# =============================================================================
# Helpers unique to this test file
# =============================================================================


async def _add_message(
    session: AsyncSession,
    room_id: str,
    user_id: Any,
    content: str,
    created_at: datetime | None = None,
) -> Message:
    msg = Message(room_id=room_id, user_id=user_id, content=content)
    if created_at is not None:
        msg.created_at = created_at
    session.add(msg)
    await session.commit()
    await session.refresh(msg)
    return msg


# =============================================================================
# POST (create message)
# =============================================================================


@pytest.mark.asyncio
async def test_create_message(
    client: AsyncClient, session: AsyncSession, mock_sio: MockSioServer
) -> None:
    """POST creates message, returns correct fields, and broadcasts."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/chat/messages",
        json={"content": "Hello!"},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    data = response.json()
    assert data["content"] == "Hello!"
    assert data["room_id"] == room.id
    assert data["email"] == user.email
    assert data["updated_at"] is None
    assert isinstance(data["id"], int)
    assert "created_at" in data

    # Verify socket broadcast
    assert len(mock_sio.emitted) == 1
    assert mock_sio.emitted[0]["event"] == "message_new"


@pytest.mark.asyncio
async def test_create_message_empty_content(
    client: AsyncClient, session: AsyncSession
) -> None:
    """POST with empty content returns 422."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/chat/messages",
        json={"content": ""},
        headers=auth_header(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_create_message_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """POST without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/chat/messages",
        json={"content": "Hello!"},
    )
    assert response.status_code == 401


# =============================================================================
# GET (list messages)
# =============================================================================


@pytest.mark.asyncio
async def test_list_messages_empty(
    client: AsyncClient, session: AsyncSession
) -> None:
    """GET returns empty list for room with no messages."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/chat/messages", headers=auth_header(token)
    )
    assert response.status_code == 200
    data = response.json()
    assert data["items"] == []
    assert data["metadata"]["has_more"] is False
    assert data["metadata"]["total_count"] == 0


@pytest.mark.asyncio
async def test_list_messages_returns_newest_first(
    client: AsyncClient, session: AsyncSession
) -> None:
    """GET returns messages ordered by created_at descending."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await _add_message(
        session,
        room.id,
        user.id,
        "First",
        datetime(2024, 1, 1, tzinfo=UTC),
    )
    await _add_message(
        session,
        room.id,
        user.id,
        "Second",
        datetime(2024, 1, 2, tzinfo=UTC),
    )
    await _add_message(
        session,
        room.id,
        user.id,
        "Third",
        datetime(2024, 1, 3, tzinfo=UTC),
    )

    response = await client.get(
        f"/v1/rooms/{room.id}/chat/messages", headers=auth_header(token)
    )
    assert response.status_code == 200
    data = response.json()
    contents = [m["content"] for m in data["items"]]
    assert contents == ["Third", "Second", "First"]
    assert data["metadata"]["total_count"] == 3
    assert data["metadata"]["has_more"] is False


@pytest.mark.asyncio
async def test_list_messages_pagination(
    client: AsyncClient, session: AsyncSession
) -> None:
    """before cursor returns older messages, has_more is correct."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    timestamps = [datetime(2024, 1, i, tzinfo=UTC) for i in range(1, 6)]
    for i, ts in enumerate(timestamps):
        await _add_message(session, room.id, user.id, f"msg-{i}", ts)

    # Fetch first page (limit=2)
    response = await client.get(
        f"/v1/rooms/{room.id}/chat/messages?limit=2", headers=auth_header(token)
    )
    assert response.status_code == 200
    page1 = response.json()
    assert len(page1["items"]) == 2
    assert page1["metadata"]["has_more"] is True
    assert page1["items"][0]["content"] == "msg-4"  # newest

    # Fetch second page using oldest_timestamp cursor
    cursor = page1["metadata"]["oldest_timestamp"]
    response = await client.get(
        f"/v1/rooms/{room.id}/chat/messages?limit=2&before={cursor}",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    page2 = response.json()
    assert len(page2["items"]) == 2
    assert page2["items"][0]["content"] == "msg-2"


@pytest.mark.asyncio
async def test_list_messages_includes_email(
    client: AsyncClient, session: AsyncSession
) -> None:
    """GET populates the email field from the User table."""
    user, token = await create_test_user_in_db(session, email="alice@test.com")
    room = await create_test_room(session, user)

    await client.post(
        f"/v1/rooms/{room.id}/chat/messages",
        json={"content": "Hi"},
        headers=auth_header(token),
    )

    response = await client.get(
        f"/v1/rooms/{room.id}/chat/messages", headers=auth_header(token)
    )
    assert response.status_code == 200
    assert response.json()["items"][0]["email"] == "alice@test.com"


# =============================================================================
# PATCH (edit message)
# =============================================================================


@pytest.mark.asyncio
async def test_edit_message(
    client: AsyncClient, session: AsyncSession, mock_sio: MockSioServer
) -> None:
    """PATCH updates content and sets updated_at."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    msg = await _add_message(session, room.id, user.id, "Original")

    response = await client.patch(
        f"/v1/rooms/{room.id}/chat/messages/{msg.id}",
        json={"content": "Edited"},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["content"] == "Edited"
    assert data["updated_at"] is not None

    # Verify broadcast
    edit_events = [e for e in mock_sio.emitted if e["event"] == "message_edited"]
    assert len(edit_events) == 1


@pytest.mark.asyncio
async def test_edit_message_ownership(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Cannot edit another user's message (403)."""
    user1, _ = await create_test_user_in_db(session, email="user1@test.com")
    _, token2 = await create_test_user_in_db(session, email="user2@test.com")
    room = await create_test_room(session, user1)

    msg = await _add_message(session, room.id, user1.id, "User1's msg")

    response = await client.patch(
        f"/v1/rooms/{room.id}/chat/messages/{msg.id}",
        json={"content": "Hacked"},
        headers=auth_header(token2),
    )
    assert response.status_code == 403
    assert "not-message-owner" in response.json()["type"]


@pytest.mark.asyncio
async def test_edit_message_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """404 for nonexistent message."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.patch(
        f"/v1/rooms/{room.id}/chat/messages/99999",
        json={"content": "Edit"},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "message-not-found" in response.json()["type"]


# =============================================================================
# Room not found
# =============================================================================


@pytest.mark.asyncio
async def test_list_messages_room_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """GET for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/nonexistent/chat/messages", headers=auth_header(token)
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
