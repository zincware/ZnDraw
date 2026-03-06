"""Tests for Socket.IO room functionality."""

import asyncio
import uuid

import pytest
import socketio
from httpx import AsyncClient
from zndraw_socketio import wrap

from zndraw.exceptions import NotInRoom, RoomNotFound
from zndraw.schemas import PresenceResponse, RoomCreateResponse
from zndraw.socket_events import (
    Heartbeat,
    HeartbeatResponse,
    RoomJoin,
    RoomJoinResponse,
    RoomLeave,
    RoomLeaveResponse,
    SessionJoined,
    SessionLeft,
    TypingResponse,
    TypingStart,
    TypingStop,
    UserGet,
    UserGetResponse,
)


async def _get_user_token(
    http_client: AsyncClient, email: str, password: str = "testpassword"
) -> str:
    """Register a user and get their auth token."""
    reg_response = await http_client.post(
        "/v1/auth/register", json={"email": email, "password": password}
    )
    assert reg_response.status_code == 201, f"Register failed: {reg_response.text}"
    login_response = await http_client.post(
        "/v1/auth/jwt/login",
        data={"username": email, "password": password},
    )
    assert login_response.status_code == 200, f"Login failed: {login_response.text}"
    return login_response.json()["access_token"]


async def _create_room(http_client: AsyncClient, token: str) -> str:
    """Create a room and return its ID."""
    room_id = str(uuid.uuid4())
    response = await http_client.post(
        "/v1/rooms",
        json={"room_id": room_id},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 201
    RoomCreateResponse.model_validate(response.json())
    return room_id


@pytest.mark.asyncio
async def test_socketio_join_room(server: str, http_client: AsyncClient) -> None:
    """Test joining a Socket.IO room."""
    token = await _get_user_token(http_client, "roomuser@example.com")
    room_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    assert sio_client.connected

    tsio_client = wrap(sio_client)
    result = await tsio_client.call(
        RoomJoin(room_id=room_id), response_model=RoomJoinResponse
    )
    assert result.room_id == room_id

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_join_nonexistent_room(
    server: str, http_client: AsyncClient
) -> None:
    """Test joining a non-existent room returns error."""
    token = await _get_user_token(http_client, "roomuser2@example.com")

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})

    raw_result = await sio_client.call(
        "room_join", {"room_id": "nonexistent-room-99999"}
    )
    assert isinstance(raw_result, dict)
    assert raw_result["type"] == RoomNotFound.type_uri()

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_leave_room(server: str, http_client: AsyncClient) -> None:
    """Test leaving a Socket.IO room."""
    token = await _get_user_token(http_client, "leaveuser@example.com")
    room_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)
    await tsio_client.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    result = await tsio_client.call(
        RoomLeave(room_id=room_id), response_model=RoomLeaveResponse
    )
    assert result.room_id == room_id

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_leave_room_not_in_room(
    server: str, http_client: AsyncClient
) -> None:
    """Test leaving a room when not in it returns error."""
    token = await _get_user_token(http_client, "notinroom@example.com")

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})

    raw_result = await sio_client.call("room_leave", {"room_id": "room-1"})
    assert isinstance(raw_result, dict)
    assert raw_result["type"] == NotInRoom.type_uri()

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_leave_room_after_switch_is_idempotent(
    server: str, http_client: AsyncClient
) -> None:
    """Test leaving a room after switching to another succeeds (idempotent).

    This simulates the race condition where frontend cleanup runs after
    backend has already switched to a new room via room_join.
    """
    token = await _get_user_token(http_client, "roomswitch@example.com")
    room1_id = await _create_room(http_client, token)
    room2_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)

    # Join room1
    await tsio_client.call(RoomJoin(room_id=room1_id), response_model=RoomJoinResponse)

    # Switch to room2 (backend automatically leaves room1)
    await tsio_client.call(RoomJoin(room_id=room2_id), response_model=RoomJoinResponse)

    # Try to leave room1 again (simulates frontend cleanup race)
    # Should succeed idempotently since we're already in room2
    result = await tsio_client.call(
        RoomLeave(room_id=room1_id), response_model=RoomLeaveResponse
    )
    assert result.room_id == room1_id

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_join_system_room_overview(
    server: str, http_client: AsyncClient
) -> None:
    """Test joining @overview system room (no DB validation, no presence)."""
    token = await _get_user_token(http_client, "overview@example.com")

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)

    # Join @overview system room (does not require DB record)
    result = await tsio_client.call(
        RoomJoin(room_id="@overview", client_type="frontend"),
        response_model=RoomJoinResponse,
    )
    assert result.room_id == "@overview"
    assert result.session_id is not None  # Server assigns session ID
    assert result.step == 0  # System rooms have no frames
    assert result.frame_count == 0  # System rooms have no frames
    assert result.locked is False  # System rooms are never locked

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_heartbeat(server: str, http_client: AsyncClient) -> None:
    """Test heartbeat updates presence."""
    token = await _get_user_token(http_client, "heartbeatuser@example.com")
    room_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)
    await tsio_client.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    result = await tsio_client.call(
        Heartbeat(room_id=room_id), response_model=HeartbeatResponse
    )
    assert result.status == "ok"

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_heartbeat_wrong_room(
    server: str, http_client: AsyncClient
) -> None:
    """Test heartbeat for wrong room returns error."""
    token = await _get_user_token(http_client, "wronghb@example.com")
    room_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)
    await tsio_client.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    raw_result = await sio_client.call(
        "heartbeat", {"room_id": "nonexistent-room-99999"}
    )
    assert isinstance(raw_result, dict)
    assert raw_result["type"] == NotInRoom.type_uri()

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_typing_events(server: str, http_client: AsyncClient) -> None:
    """Test typing start/stop events."""
    token = await _get_user_token(http_client, "typinguser@example.com")
    room_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)
    await tsio_client.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    result = await tsio_client.call(
        TypingStart(room_id=room_id), response_model=TypingResponse
    )
    assert result.status == "ok"

    result = await tsio_client.call(
        TypingStop(room_id=room_id), response_model=TypingResponse
    )
    assert result.status == "ok"

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_session_joined_broadcast(
    server: str, http_client: AsyncClient
) -> None:
    """Test that session_joined is broadcast when someone joins."""
    token1 = await _get_user_token(http_client, "first@example.com")
    token2 = await _get_user_token(http_client, "second@example.com")

    room_id = await _create_room(http_client, token1)

    # User1 connects and joins Socket.IO room
    sio_client1 = socketio.AsyncClient()
    received_events: list[SessionJoined] = []

    @sio_client1.on("session_joined")
    async def on_session_joined(data: dict) -> None:
        received_events.append(SessionJoined.model_validate(data))

    await sio_client1.connect(server, auth={"token": token1})
    tsio_client1 = wrap(sio_client1)
    await tsio_client1.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    # User2 connects and joins Socket.IO room
    sio_client2 = socketio.AsyncClient()
    await sio_client2.connect(server, auth={"token": token2})
    tsio_client2 = wrap(sio_client2)

    # Get user2's user_id via UserGet
    user2_info = await tsio_client2.call(UserGet(), response_model=UserGetResponse)

    await tsio_client2.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    await asyncio.sleep(0.5)

    assert len(received_events) == 1
    assert received_events[0].room_id == room_id
    assert received_events[0].user_id == user2_info.id
    assert received_events[0].sid  # sid should be present

    await sio_client1.disconnect()
    await sio_client2.disconnect()


@pytest.mark.asyncio
async def test_socketio_presence(server: str, http_client: AsyncClient) -> None:
    """Test presence endpoint reflects online users."""
    token = await _get_user_token(http_client, "presenceuser@example.com")
    room_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)
    await tsio_client.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    # Check presence via REST
    presence_response = await http_client.get(
        f"/v1/rooms/{room_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert presence_response.status_code == 200
    presence = PresenceResponse.model_validate(presence_response.json())
    assert len(presence.items) == 1
    session_presence = presence.items[0]
    assert session_presence.email == "presenceuser@example.com"
    assert session_presence.sid  # SID should be present

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_session_left_broadcast(
    server: str, http_client: AsyncClient
) -> None:
    """Test that session_left is broadcast when someone leaves."""
    token1 = await _get_user_token(http_client, "stayer@example.com")
    token2 = await _get_user_token(http_client, "leaver2@example.com")

    room_id = await _create_room(http_client, token1)

    # User1 connects and joins Socket.IO room
    sio_client1 = socketio.AsyncClient()
    received_events: list[SessionLeft] = []

    @sio_client1.on("session_left")
    async def on_session_left(data: dict) -> None:
        received_events.append(SessionLeft.model_validate(data))

    await sio_client1.connect(server, auth={"token": token1})
    tsio_client1 = wrap(sio_client1)
    await tsio_client1.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    # User2 connects, joins, then leaves Socket.IO room
    sio_client2 = socketio.AsyncClient()
    await sio_client2.connect(server, auth={"token": token2})
    tsio_client2 = wrap(sio_client2)

    user2_info = await tsio_client2.call(UserGet(), response_model=UserGetResponse)

    await tsio_client2.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)
    await tsio_client2.call(
        RoomLeave(room_id=room_id), response_model=RoomLeaveResponse
    )

    await asyncio.sleep(0.5)

    assert len(received_events) == 1
    assert received_events[0].room_id == room_id
    assert received_events[0].user_id == user2_info.id
    assert received_events[0].sid  # sid should be present

    await sio_client1.disconnect()
    await sio_client2.disconnect()


@pytest.mark.asyncio
async def test_socketio_session_left_on_disconnect(
    server: str, http_client: AsyncClient
) -> None:
    """Test that session_left is broadcast when someone disconnects."""
    token1 = await _get_user_token(http_client, "observer@example.com")
    token2 = await _get_user_token(http_client, "disconnector@example.com")

    room_id = await _create_room(http_client, token1)

    # User1 connects and joins Socket.IO room
    sio_client1 = socketio.AsyncClient()
    received_events: list[SessionLeft] = []

    @sio_client1.on("session_left")
    async def on_session_left(data: dict) -> None:
        received_events.append(SessionLeft.model_validate(data))

    await sio_client1.connect(server, auth={"token": token1})
    tsio_client1 = wrap(sio_client1)
    await tsio_client1.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    # User2 connects, joins, then disconnects (abruptly)
    sio_client2 = socketio.AsyncClient()
    await sio_client2.connect(server, auth={"token": token2})
    tsio_client2 = wrap(sio_client2)

    user2_info = await tsio_client2.call(UserGet(), response_model=UserGetResponse)

    await tsio_client2.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)
    await sio_client2.disconnect()

    await asyncio.sleep(0.5)

    assert len(received_events) == 1
    assert received_events[0].room_id == room_id
    assert received_events[0].user_id == user2_info.id
    assert received_events[0].sid  # sid should be present

    await sio_client1.disconnect()


# =============================================================================
# Multi-Session Tests (Session-Centric Presence)
# =============================================================================


@pytest.mark.asyncio
async def test_same_user_two_sessions_same_room(
    server: str, http_client: AsyncClient
) -> None:
    """User in same room from two sessions - both tracked independently."""
    token = await _get_user_token(http_client, "multitab@example.com")
    room_id = await _create_room(http_client, token)

    # Connect session 1 (first tab)
    sio_client1 = socketio.AsyncClient()
    session1_joined_events: list[SessionJoined] = []
    session1_left_events: list[SessionLeft] = []

    @sio_client1.on("session_joined")
    async def on_session1_joined(data: dict) -> None:
        session1_joined_events.append(SessionJoined.model_validate(data))

    @sio_client1.on("session_left")
    async def on_session1_left(data: dict) -> None:
        session1_left_events.append(SessionLeft.model_validate(data))

    await sio_client1.connect(server, auth={"token": token})
    tsio_client1 = wrap(sio_client1)

    user_info = await tsio_client1.call(UserGet(), response_model=UserGetResponse)

    await tsio_client1.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    # Connect session 2 (second tab, same user)
    sio_client2 = socketio.AsyncClient()
    await sio_client2.connect(server, auth={"token": token})
    tsio_client2 = wrap(sio_client2)
    await tsio_client2.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    await asyncio.sleep(0.3)

    # Session 1 should have received session_joined for session 2
    assert len(session1_joined_events) == 1
    assert session1_joined_events[0].user_id == user_info.id
    assert session1_joined_events[0].sid  # Different SID than session 1

    # Check presence - should show 2 sessions (same user, different SIDs)
    presence_response = await http_client.get(
        f"/v1/rooms/{room_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert presence_response.status_code == 200
    presence = PresenceResponse.model_validate(presence_response.json())
    assert len(presence.items) == 2  # Two sessions for same user
    sids = {s.sid for s in presence.items}
    assert len(sids) == 2  # Different SIDs
    assert all(s.user_id == user_info.id for s in presence.items)

    # Disconnect session 2 - session 1 should receive session_left
    await sio_client2.disconnect()
    await asyncio.sleep(0.3)

    assert len(session1_left_events) == 1
    assert session1_left_events[0].user_id == user_info.id
    assert session1_left_events[0].sid  # SID of session 2

    # Presence should show 1 session (session 1 is active)
    presence_response = await http_client.get(
        f"/v1/rooms/{room_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    presence = PresenceResponse.model_validate(presence_response.json())
    assert len(presence.items) == 1
    assert presence.items[0].user_id == user_info.id

    # Disconnect session 1 - presence should be empty
    await sio_client1.disconnect()
    await asyncio.sleep(0.3)

    presence_response = await http_client.get(
        f"/v1/rooms/{room_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    presence = PresenceResponse.model_validate(presence_response.json())
    assert len(presence.items) == 0


@pytest.mark.asyncio
async def test_same_user_two_sessions_different_rooms(
    server: str, http_client: AsyncClient
) -> None:
    """User in different rooms from two sessions."""
    token = await _get_user_token(http_client, "multiroom@example.com")

    room_a_id = await _create_room(http_client, token)
    room_b_id = await _create_room(http_client, token)

    # Session 1 joins Room A
    sio_client1 = socketio.AsyncClient()
    await sio_client1.connect(server, auth={"token": token})
    tsio_client1 = wrap(sio_client1)

    user_info = await tsio_client1.call(UserGet(), response_model=UserGetResponse)

    await tsio_client1.call(
        RoomJoin(room_id=room_a_id), response_model=RoomJoinResponse
    )

    # Session 2 joins Room B
    sio_client2 = socketio.AsyncClient()
    await sio_client2.connect(server, auth={"token": token})
    tsio_client2 = wrap(sio_client2)
    await tsio_client2.call(
        RoomJoin(room_id=room_b_id), response_model=RoomJoinResponse
    )

    await asyncio.sleep(0.3)

    # Check presence in Room A - should show 1 session
    presence_a = await http_client.get(
        f"/v1/rooms/{room_a_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    presence_a_data = PresenceResponse.model_validate(presence_a.json())
    assert len(presence_a_data.items) == 1
    assert presence_a_data.items[0].user_id == user_info.id

    # Check presence in Room B - should show 1 session
    presence_b = await http_client.get(
        f"/v1/rooms/{room_b_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    presence_b_data = PresenceResponse.model_validate(presence_b.json())
    assert len(presence_b_data.items) == 1
    assert presence_b_data.items[0].user_id == user_info.id

    # Disconnect session 1 (Room A)
    await sio_client1.disconnect()
    await asyncio.sleep(0.3)

    # Room A should be empty now
    presence_a = await http_client.get(
        f"/v1/rooms/{room_a_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    presence_a_data = PresenceResponse.model_validate(presence_a.json())
    assert len(presence_a_data.items) == 0

    # Room B should still show session
    presence_b = await http_client.get(
        f"/v1/rooms/{room_b_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    presence_b_data = PresenceResponse.model_validate(presence_b.json())
    assert len(presence_b_data.items) == 1
    assert presence_b_data.items[0].user_id == user_info.id

    await sio_client2.disconnect()


@pytest.mark.asyncio
async def test_session_broadcast_includes_sid(
    server: str, http_client: AsyncClient
) -> None:
    """Session broadcasts include SID for callback targeting."""
    token1 = await _get_user_token(http_client, "listener@example.com")
    token2 = await _get_user_token(http_client, "joiner@example.com")

    room_id = await _create_room(http_client, token1)

    # User1 connects and joins Socket.IO room
    sio_client1 = socketio.AsyncClient()
    received_events: list[SessionJoined] = []

    @sio_client1.on("session_joined")
    async def on_session_joined(data: dict) -> None:
        received_events.append(SessionJoined.model_validate(data))

    await sio_client1.connect(server, auth={"token": token1})
    tsio_client1 = wrap(sio_client1)
    await tsio_client1.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    # User2 connects and joins Socket.IO room
    sio_client2 = socketio.AsyncClient()
    await sio_client2.connect(server, auth={"token": token2})
    tsio_client2 = wrap(sio_client2)

    user2_info = await tsio_client2.call(UserGet(), response_model=UserGetResponse)

    await tsio_client2.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    await asyncio.sleep(0.3)

    # Verify broadcast includes all expected fields
    assert len(received_events) == 1
    event = received_events[0]
    assert event.room_id == room_id
    assert event.user_id == user2_info.id
    assert event.email == "joiner@example.com"
    assert event.sid is not None
    assert len(event.sid) > 0  # SID should be a non-empty string

    await sio_client1.disconnect()
    await sio_client2.disconnect()


# =============================================================================
# TTL Expiration Tests
# =============================================================================


@pytest.mark.asyncio
async def test_presence_expires_without_heartbeat(
    server_short_ttl: str, http_client_short_ttl: AsyncClient
) -> None:
    """Test that presence expires after TTL without heartbeat.

    Uses server_short_ttl fixture which configures a 2-second presence TTL.
    """
    token = await _get_user_token(http_client_short_ttl, "ttluser@example.com")
    room_id = await _create_room(http_client_short_ttl, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server_short_ttl, auth={"token": token})
    tsio_client = wrap(sio_client)
    await tsio_client.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    # Verify presence is set
    presence_response = await http_client_short_ttl.get(
        f"/v1/rooms/{room_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert presence_response.status_code == 200
    presence = PresenceResponse.model_validate(presence_response.json())
    assert len(presence.items) == 1

    # Wait for TTL to expire (2s TTL + buffer)
    await asyncio.sleep(2.5)

    # Presence should be empty now (no heartbeat was sent)
    presence_response = await http_client_short_ttl.get(
        f"/v1/rooms/{room_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    presence = PresenceResponse.model_validate(presence_response.json())
    assert len(presence.items) == 0

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_heartbeat_extends_presence_ttl(
    server_short_ttl: str, http_client_short_ttl: AsyncClient
) -> None:
    """Test that heartbeat keeps presence alive past initial TTL.

    Uses server_short_ttl fixture which configures a 2-second presence TTL.
    """
    token = await _get_user_token(http_client_short_ttl, "heartbeatttl@example.com")
    room_id = await _create_room(http_client_short_ttl, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server_short_ttl, auth={"token": token})
    tsio_client = wrap(sio_client)
    await tsio_client.call(RoomJoin(room_id=room_id), response_model=RoomJoinResponse)

    # Wait 1.5s (less than 2s TTL), then send heartbeat
    await asyncio.sleep(1.5)
    result = await tsio_client.call(
        Heartbeat(room_id=room_id), response_model=HeartbeatResponse
    )
    assert result.status == "ok"

    # Wait another 1.5s (total 3s from join, but only 1.5s from heartbeat)
    await asyncio.sleep(1.5)

    # Presence should still be set (heartbeat refreshed TTL)
    presence_response = await http_client_short_ttl.get(
        f"/v1/rooms/{room_id}/presence",
        headers={"Authorization": f"Bearer {token}"},
    )
    presence = PresenceResponse.model_validate(presence_response.json())
    assert len(presence.items) == 1

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_frontend_join_creates_session_camera(
    server: str, http_client: AsyncClient
) -> None:
    """Frontend clients get a session camera on join."""
    token = await _get_user_token(http_client, "cam-frontend@example.com")
    room_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)

    result = await tsio_client.call(
        RoomJoin(room_id=room_id, client_type="frontend"),
        response_model=RoomJoinResponse,
    )
    assert result.camera_key is not None
    assert result.camera_key.startswith("cam:")

    # Verify camera geometry exists via REST
    response = await http_client.get(
        f"/v1/rooms/{room_id}/geometries/{result.camera_key}",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert response.status_code == 200

    await sio_client.disconnect()


@pytest.mark.asyncio
async def test_pyclient_join_does_not_create_session_camera(
    server: str, http_client: AsyncClient
) -> None:
    """Python clients must not get a session camera on join."""
    token = await _get_user_token(http_client, "cam-pyclient@example.com")
    room_id = await _create_room(http_client, token)

    sio_client = socketio.AsyncClient()
    await sio_client.connect(server, auth={"token": token})
    tsio_client = wrap(sio_client)

    result = await tsio_client.call(
        RoomJoin(room_id=room_id, client_type="pyclient"),
        response_model=RoomJoinResponse,
    )
    assert result.camera_key is None

    await sio_client.disconnect()


# =============================================================================
# Session Camera Ownership Tests
# =============================================================================


@pytest.mark.asyncio
async def test_rest_rejects_updating_other_users_session_camera(
    server_auth: str, http_client_auth: AsyncClient
) -> None:
    """REST PUT to another user's session camera returns 403."""
    token_a = await _get_user_token(http_client_auth, "cam-rest-owner@example.com")
    token_b = await _get_user_token(http_client_auth, "cam-rest-intruder@example.com")
    room_id = await _create_room(http_client_auth, token_a)

    # User A joins as frontend (creates session camera)
    sio_a = socketio.AsyncClient()
    await sio_a.connect(server_auth, auth={"token": token_a})
    tsio_a = wrap(sio_a)
    result_a = await tsio_a.call(
        RoomJoin(room_id=room_id, client_type="frontend"),
        response_model=RoomJoinResponse,
    )
    assert result_a.camera_key is not None

    # User B tries to update User A's camera via REST → 403
    response = await http_client_auth.put(
        f"/v1/rooms/{room_id}/geometries/{result_a.camera_key}",
        json={
            "type": "Camera",
            "data": {"position": [99.0, 99.0, 99.0]},
        },
        headers={"Authorization": f"Bearer {token_b}"},
    )
    assert response.status_code == 403
    assert "forbidden" in response.json()["type"]

    await sio_a.disconnect()
