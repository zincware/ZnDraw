"""Tests for provider-aware frame dispatch in GET /v1/rooms/{room_id}/frames/{index}."""

import asyncio

import msgpack
import pytest
from helpers import (
    InMemoryResultBackend,
    auth_header,
    create_test_room,
    create_test_user_in_db,
    decode_msgpack_response,
    make_raw_frame,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.exceptions import FrameNotFound, ProblemDetail
from zndraw.storage import FrameStorage, RawFrame
from zndraw_auth import User
from zndraw_joblib.dependencies import request_hash
from zndraw_joblib.exceptions import ProviderTimeout
from zndraw_joblib.models import ProviderRecord, Worker

# =============================================================================
# Helpers
# =============================================================================


async def _create_provider(
    session: AsyncSession, room_id: str, user: User
) -> ProviderRecord:
    """Create a frames provider for the given room."""
    worker = Worker(user_id=user.id)
    session.add(worker)
    await session.flush()
    provider = ProviderRecord(
        room_id=room_id,
        category="frames",
        name=f"mount-{room_id}",
        content_type="application/x-msgpack",
        user_id=user.id,
        worker_id=worker.id,
    )
    session.add(provider)
    await session.commit()
    await session.refresh(provider)
    return provider


# =============================================================================
# Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_frame_storage_hit_ignores_provider(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Frame in storage, provider exists -- returns frame (no provider touch)."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _create_provider(session, room.id, user)

    await frame_storage[room.id].extend([make_raw_frame({"a": 1})])

    response = await client.get(
        f"/v1/rooms/{room.id}/frames/0", headers=auth_header(token)
    )
    assert response.status_code == 200
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert msgpack.unpackb(frames[0][b"a"]) == 1


@pytest.mark.asyncio
async def test_get_frame_provider_cache_hit(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
    result_backend: InMemoryResultBackend,
) -> None:
    """Frame in provider cache, storage slot is None -- returns 200 with frame."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    provider = await _create_provider(session, room.id, user)

    # Reserve slots (provider has 3 frames), slot 0 is None
    await frame_storage[room.id].reserve(3)

    # Pre-populate provider cache
    positions_packed = msgpack.packb([1, 2, 3])
    assert isinstance(positions_packed, bytes)
    frame_data: RawFrame = {b"arrays.positions": positions_packed}
    packed = msgpack.packb(frame_data, use_bin_type=True)
    assert isinstance(packed, bytes)
    params = {"index": "0"}
    rhash = request_hash(params)
    cache_key = f"provider-result:{provider.full_name}:{rhash}"
    await result_backend.store(cache_key, packed, 300)

    response = await client.get(
        f"/v1/rooms/{room.id}/frames/0", headers=auth_header(token)
    )
    assert response.status_code == 200
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert b"arrays.positions" in frames[0]


@pytest.mark.asyncio
async def test_get_frame_provider_timeout(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Frame not cached, provider exists -- long-poll times out → 504."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _create_provider(session, room.id, user)

    # Reserve slots, leave them empty
    await frame_storage[room.id].reserve(5)

    response = await client.get(
        f"/v1/rooms/{room.id}/frames/2", headers=auth_header(token)
    )
    assert response.status_code == 504

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == ProviderTimeout.type_uri()
    assert response.headers.get("retry-after") == "2"


@pytest.mark.asyncio
async def test_get_frame_no_provider_returns_404(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Frame missing, no provider registered -- returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Reserve slots but no provider registered
    await frame_storage[room.id].reserve(3)

    response = await client.get(
        f"/v1/rooms/{room.id}/frames/1", headers=auth_header(token)
    )
    assert response.status_code == 404
    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


@pytest.mark.asyncio
async def test_get_frame_dispatch_acquires_inflight(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
    result_backend: InMemoryResultBackend,
) -> None:
    """After dispatch, inflight lock is acquired."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    provider = await _create_provider(session, room.id, user)

    await frame_storage[room.id].reserve(3)

    response = await client.get(
        f"/v1/rooms/{room.id}/frames/0", headers=auth_header(token)
    )
    assert response.status_code == 504  # timeout, but inflight lock was set

    # Verify inflight lock was acquired
    params = {"index": "0"}
    rhash = request_hash(params)
    inflight_key = f"provider-inflight:{provider.full_name}:{rhash}"
    assert inflight_key in result_backend._inflight


@pytest.mark.asyncio
async def test_get_frame_notify_wakes_long_poll(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
    result_backend: InMemoryResultBackend,
) -> None:
    """Provider uploads result mid-poll — long-poll wakes up and returns 200."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    provider = await _create_provider(session, room.id, user)

    await frame_storage[room.id].reserve(3)

    # Build the cache key that _dispatch_provider_frame will wait on
    params = {"index": "1"}
    rhash = request_hash(params)
    cache_key = f"provider-result:{provider.full_name}:{rhash}"

    # Prepare the frame data that the "provider" will upload
    positions_packed = msgpack.packb([10, 20, 30])
    assert isinstance(positions_packed, bytes)
    frame_data: RawFrame = {b"arrays.positions": positions_packed}
    packed = msgpack.packb(frame_data, use_bin_type=True)
    assert isinstance(packed, bytes)

    async def _simulate_provider_upload() -> None:
        """Simulate provider uploading result after a short delay."""
        await asyncio.sleep(0.1)
        await result_backend.store(cache_key, packed, 300)
        await result_backend.notify_key(cache_key)

    # Start the simulated provider upload concurrently with the GET request
    upload_task = asyncio.create_task(_simulate_provider_upload())

    response = await client.get(
        f"/v1/rooms/{room.id}/frames/1", headers=auth_header(token)
    )

    await upload_task

    # The long-poll should have woken up and returned 200
    assert response.status_code == 200
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert b"arrays.positions" in frames[0]
    assert msgpack.unpackb(frames[0][b"arrays.positions"]) == [10, 20, 30]


@pytest.mark.asyncio
async def test_list_frames_notify_wakes_concurrent_dispatch(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
    result_backend: InMemoryResultBackend,
) -> None:
    """Multiple missing frames dispatched concurrently — all wake on notify."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    provider = await _create_provider(session, room.id, user)

    # Reserve 3 slots, fill only index 1
    await frame_storage[room.id].reserve(3)
    await frame_storage[room.id][1].set(make_raw_frame({"existing": 1}))

    # Build cache keys for the missing frames (index 0 and 2)
    def _cache_key(idx: int) -> str:
        return (
            f"provider-result:{provider.full_name}:{request_hash({'index': str(idx)})}"
        )

    async def _simulate_provider_upload() -> None:
        """Upload results for both missing frames after a short delay."""
        await asyncio.sleep(0.1)
        for idx in (0, 2):
            idx_packed = msgpack.packb(idx)
            assert isinstance(idx_packed, bytes)
            frame: RawFrame = {b"idx": idx_packed}
            packed = msgpack.packb(frame, use_bin_type=True)
            assert isinstance(packed, bytes)
            key = _cache_key(idx)
            await result_backend.store(key, packed, 300)
            await result_backend.notify_key(key)

    upload_task = asyncio.create_task(_simulate_provider_upload())

    response = await client.get(
        f"/v1/rooms/{room.id}/frames?indices=0,1,2", headers=auth_header(token)
    )

    await upload_task

    assert response.status_code == 200
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 3
    # Index 0 and 2 came from provider, index 1 from storage
    assert msgpack.unpackb(frames[0][b"idx"]) == 0
    assert msgpack.unpackb(frames[1][b"existing"]) == 1
    assert msgpack.unpackb(frames[2][b"idx"]) == 2
