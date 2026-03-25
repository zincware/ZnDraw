"""Tests for Chat REST API endpoints."""

from collections.abc import AsyncIterator
from datetime import UTC, datetime
from typing import Any
from unittest.mock import AsyncMock

import pytest
import pytest_asyncio
from conftest import MockSioServer, create_test_token, create_test_user_model
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel
from zndraw_auth import User
from zndraw_auth.settings import AuthSettings

from zndraw.config import Settings
from zndraw.models import MemberRole, Message, Room, RoomMembership

# =============================================================================
# Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="chat_session")
async def chat_session_fixture() -> AsyncIterator[AsyncSession]:
    """Create a fresh in-memory async database session for each test."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)
        factory = async_sessionmaker(
            bind=engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            yield session
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="mock_sio")
async def mock_sio_fixture() -> MockSioServer:
    return MockSioServer()


@pytest_asyncio.fixture(name="chat_client")
async def chat_client_fixture(
    chat_session: AsyncSession, mock_sio: MockSioServer
) -> AsyncIterator[AsyncClient]:
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield chat_session

    def get_sio_override() -> MockSioServer:
        return mock_sio

    # Mock Redis for WritableRoomDep (returns None = no edit lock)
    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_tsio] = get_sio_override
    app.dependency_overrides[get_redis] = lambda: mock_redis
    app.state.settings = Settings()
    app.state.auth_settings = AuthSettings()

    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        yield client

    app.dependency_overrides.clear()


async def _create_user(
    session: AsyncSession, email: str = "testuser@local.test"
) -> tuple[User, str]:
    user = create_test_user_model(email=email)
    session.add(user)
    await session.commit()
    await session.refresh(user)
    token = create_test_token(user)
    return user, token


async def _create_room(session: AsyncSession, user: User) -> Room:
    room = Room(description="Test Room", created_by_id=user.id, is_public=True)  # type: ignore[arg-type]
    session.add(room)
    await session.commit()
    await session.refresh(room)
    membership = RoomMembership(
        room_id=room.id,
        user_id=user.id,
        role=MemberRole.OWNER,  # type: ignore[arg-type]
    )
    session.add(membership)
    await session.commit()
    return room


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


def _auth(token: str) -> dict[str, str]:
    return {"Authorization": f"Bearer {token}"}


# =============================================================================
# POST (create message)
# =============================================================================


@pytest.mark.asyncio
async def test_create_message(
    chat_client: AsyncClient, chat_session: AsyncSession, mock_sio: MockSioServer
) -> None:
    """POST creates message, returns correct fields, and broadcasts."""
    user, token = await _create_user(chat_session)
    room = await _create_room(chat_session, user)

    response = await chat_client.post(
        f"/v1/rooms/{room.id}/chat/messages",
        json={"content": "Hello!"},
        headers=_auth(token),
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
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """POST with empty content returns 422."""
    user, token = await _create_user(chat_session)
    room = await _create_room(chat_session, user)

    response = await chat_client.post(
        f"/v1/rooms/{room.id}/chat/messages",
        json={"content": ""},
        headers=_auth(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_create_message_requires_auth(
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """POST without auth returns 401."""
    user, _ = await _create_user(chat_session)
    room = await _create_room(chat_session, user)

    response = await chat_client.post(
        f"/v1/rooms/{room.id}/chat/messages",
        json={"content": "Hello!"},
    )
    assert response.status_code == 401


# =============================================================================
# GET (list messages)
# =============================================================================


@pytest.mark.asyncio
async def test_list_messages_empty(
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """GET returns empty list for room with no messages."""
    user, token = await _create_user(chat_session)
    room = await _create_room(chat_session, user)

    response = await chat_client.get(
        f"/v1/rooms/{room.id}/chat/messages", headers=_auth(token)
    )
    assert response.status_code == 200
    data = response.json()
    assert data["items"] == []
    assert data["metadata"]["has_more"] is False
    assert data["metadata"]["total_count"] == 0


@pytest.mark.asyncio
async def test_list_messages_returns_newest_first(
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """GET returns messages ordered by created_at descending."""
    user, token = await _create_user(chat_session)
    room = await _create_room(chat_session, user)

    await _add_message(
        chat_session,
        room.id,
        user.id,
        "First",
        datetime(2024, 1, 1, tzinfo=UTC),
    )
    await _add_message(
        chat_session,
        room.id,
        user.id,
        "Second",
        datetime(2024, 1, 2, tzinfo=UTC),
    )
    await _add_message(
        chat_session,
        room.id,
        user.id,
        "Third",
        datetime(2024, 1, 3, tzinfo=UTC),
    )

    response = await chat_client.get(
        f"/v1/rooms/{room.id}/chat/messages", headers=_auth(token)
    )
    assert response.status_code == 200
    data = response.json()
    contents = [m["content"] for m in data["items"]]
    assert contents == ["Third", "Second", "First"]
    assert data["metadata"]["total_count"] == 3
    assert data["metadata"]["has_more"] is False


@pytest.mark.asyncio
async def test_list_messages_pagination(
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """before cursor returns older messages, has_more is correct."""
    user, token = await _create_user(chat_session)
    room = await _create_room(chat_session, user)

    timestamps = [datetime(2024, 1, i, tzinfo=UTC) for i in range(1, 6)]
    for i, ts in enumerate(timestamps):
        await _add_message(chat_session, room.id, user.id, f"msg-{i}", ts)

    # Fetch first page (limit=2)
    response = await chat_client.get(
        f"/v1/rooms/{room.id}/chat/messages?limit=2", headers=_auth(token)
    )
    assert response.status_code == 200
    page1 = response.json()
    assert len(page1["items"]) == 2
    assert page1["metadata"]["has_more"] is True
    assert page1["items"][0]["content"] == "msg-4"  # newest

    # Fetch second page using oldest_timestamp cursor
    cursor = page1["metadata"]["oldest_timestamp"]
    response = await chat_client.get(
        f"/v1/rooms/{room.id}/chat/messages?limit=2&before={cursor}",
        headers=_auth(token),
    )
    assert response.status_code == 200
    page2 = response.json()
    assert len(page2["items"]) == 2
    assert page2["items"][0]["content"] == "msg-2"


@pytest.mark.asyncio
async def test_list_messages_includes_email(
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """GET populates the email field from the User table."""
    user, token = await _create_user(chat_session, email="alice@test.com")
    room = await _create_room(chat_session, user)

    await chat_client.post(
        f"/v1/rooms/{room.id}/chat/messages",
        json={"content": "Hi"},
        headers=_auth(token),
    )

    response = await chat_client.get(
        f"/v1/rooms/{room.id}/chat/messages", headers=_auth(token)
    )
    assert response.status_code == 200
    assert response.json()["items"][0]["email"] == "alice@test.com"


# =============================================================================
# PATCH (edit message)
# =============================================================================


@pytest.mark.asyncio
async def test_edit_message(
    chat_client: AsyncClient, chat_session: AsyncSession, mock_sio: MockSioServer
) -> None:
    """PATCH updates content and sets updated_at."""
    user, token = await _create_user(chat_session)
    room = await _create_room(chat_session, user)

    msg = await _add_message(chat_session, room.id, user.id, "Original")

    response = await chat_client.patch(
        f"/v1/rooms/{room.id}/chat/messages/{msg.id}",
        json={"content": "Edited"},
        headers=_auth(token),
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
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """Cannot edit another user's message (403)."""
    user1, _ = await _create_user(chat_session, email="user1@test.com")
    _, token2 = await _create_user(chat_session, email="user2@test.com")
    room = await _create_room(chat_session, user1)

    msg = await _add_message(chat_session, room.id, user1.id, "User1's msg")

    response = await chat_client.patch(
        f"/v1/rooms/{room.id}/chat/messages/{msg.id}",
        json={"content": "Hacked"},
        headers=_auth(token2),
    )
    assert response.status_code == 403
    assert "not-message-owner" in response.json()["type"]


@pytest.mark.asyncio
async def test_edit_message_not_found(
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """404 for nonexistent message."""
    user, token = await _create_user(chat_session)
    room = await _create_room(chat_session, user)

    response = await chat_client.patch(
        f"/v1/rooms/{room.id}/chat/messages/99999",
        json={"content": "Edit"},
        headers=_auth(token),
    )
    assert response.status_code == 404
    assert "message-not-found" in response.json()["type"]


# =============================================================================
# Room not found
# =============================================================================


@pytest.mark.asyncio
async def test_list_messages_room_not_found(
    chat_client: AsyncClient, chat_session: AsyncSession
) -> None:
    """GET for non-existent room returns 404."""
    _, token = await _create_user(chat_session)

    response = await chat_client.get(
        "/v1/rooms/nonexistent/chat/messages", headers=_auth(token)
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
