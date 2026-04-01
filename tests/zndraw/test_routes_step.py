"""Tests for Step REST API endpoints."""

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
    make_raw_frame,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.schemas import StepResponse, StepUpdateResponse
from zndraw.storage import FrameStorage

# =============================================================================
# GET Step Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_step_returns_zero_initially(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test GET returns step=0 for new room with no step set."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add some frames to the room
    await frame_storage[room.id].extend(
        [make_raw_frame({"a": 1}), make_raw_frame({"b": 2}), make_raw_frame({"c": 3})]
    )

    response = await client.get(
        f"/v1/rooms/{room.id}/step",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = StepResponse.model_validate(response.json())
    assert result.step == 0
    assert result.total_frames == 3


@pytest.mark.asyncio
async def test_get_step_returns_current_step(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test GET returns previously set step."""
    from zndraw.models import MemberRole, Room, RoomMembership

    user, token = await create_test_user_in_db(session)
    # Create room with step=2 manually (create_test_room doesn't support step param)
    room = Room(
        description="Test Room",
        created_by_id=user.id,  # type: ignore[arg-type]
        is_public=True,
        step=2,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)
    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=user.id,  # type: ignore[arg-type]
        role=MemberRole.OWNER,
    )
    session.add(membership)
    await session.commit()

    # Add frames
    await frame_storage[room.id].extend(
        [make_raw_frame({"a": 1}), make_raw_frame({"b": 2}), make_raw_frame({"c": 3})]
    )

    response = await client.get(
        f"/v1/rooms/{room.id}/step",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = StepResponse.model_validate(response.json())
    assert result.step == 2
    assert result.total_frames == 3


# =============================================================================
# PUT Step Tests
# =============================================================================


@pytest.mark.asyncio
async def test_set_step_updates_and_returns(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
    mock_sio: MockSioServer,
) -> None:
    """Test PUT updates step and returns new value."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add frames
    await frame_storage[room.id].extend(
        [make_raw_frame({"a": 1}), make_raw_frame({"b": 2}), make_raw_frame({"c": 3})]
    )

    response = await client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": 1},
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = StepUpdateResponse.model_validate(response.json())
    assert result.success is True
    assert result.step == 1

    # Verify step persisted in DB
    await session.refresh(room)
    assert room.step == 1

    # Verify Socket.IO broadcast was sent
    assert len(mock_sio.emitted) == 1
    assert mock_sio.emitted[0]["event"] == "frame_update"
    assert mock_sio.emitted[0]["room"] == f"room:{room.id}"


@pytest.mark.asyncio
async def test_set_step_out_of_bounds_returns_422(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test PUT with step > total_frames returns 422."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add 3 frames (indices 0, 1, 2)
    await frame_storage[room.id].extend(
        [make_raw_frame({"a": 1}), make_raw_frame({"b": 2}), make_raw_frame({"c": 3})]
    )

    # Request step=100 — should return 422
    response = await client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": 100},
        headers=auth_header(token),
    )
    assert response.status_code == 422
    body = response.json()
    assert body["type"] == "/v1/problems/step-out-of-bounds"
    assert "out of bounds" in body["detail"]


@pytest.mark.asyncio
async def test_set_step_empty_room_rejects_nonzero(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test PUT to room with no frames rejects non-zero step."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Room has no frames — step=5 should be rejected
    response = await client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": 5},
        headers=auth_header(token),
    )
    assert response.status_code == 422
    body = response.json()
    assert body["type"] == "/v1/problems/step-out-of-bounds"


@pytest.mark.asyncio
async def test_set_step_negative_returns_422(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test PUT with negative step returns 422 (Pydantic ge=0 validation)."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await frame_storage[room.id].extend([make_raw_frame({"a": 1})])

    response = await client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": -1},
        headers=auth_header(token),
    )
    assert response.status_code == 422
    body = response.json()
    assert body["type"] == "/v1/problems/unprocessable-content"
    assert "body.step" in body["detail"]


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_step_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(f"/v1/rooms/{room.id}/step")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_set_step_requires_auth(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/step",
        json={"step": 1},
    )
    assert response.status_code == 401


# =============================================================================
# Room Not Found Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_step_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test GET for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/99999/step",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_set_step_returns_404_for_nonexistent_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test PUT for non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.put(
        "/v1/rooms/99999/step",
        json={"step": 1},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
