"""Tests for Frame REST API endpoints."""

from collections.abc import AsyncIterator
from typing import Any
from unittest.mock import AsyncMock

import ase
import msgpack
import pytest
import pytest_asyncio
from conftest import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
    decode_msgpack_response,
    make_raw_frame,
)
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel

from zndraw.client import atoms_to_json_dict
from zndraw.exceptions import FrameNotFound, ProblemDetail, RoomNotFound
from zndraw.schemas import (
    FrameBulkResponse,
    FrameMergeResponse,
    FrameResponse,
    StatusResponse,
)
from zndraw.storage import InMemoryStorage
from zndraw.storage.base import RawFrame


def _make_json_frame(formula: str = "H2") -> dict[str, Any]:
    """Create a JSON-serializable frame dict from an ASE Atoms formula."""
    atoms = ase.Atoms(
        formula,
        positions=[
            [i, 0, 0] for i in range(ase.Atoms(formula).get_global_number_of_atoms())
        ],
    )
    return atoms_to_json_dict(atoms)


def raw_frame_to_dict(frame: RawFrame) -> dict[str, Any]:
    """Convert a raw frame to a string-keyed dict for easier assertion.

    Keys are decoded from bytes, values are unpacked from msgpack.
    """
    return {k.decode(): msgpack.unpackb(v) for k, v in frame.items()}


# =============================================================================
# Test-specific Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="frame_session")
async def frame_session_fixture() -> AsyncIterator[AsyncSession]:
    """Create a fresh in-memory async database session for each test."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )

    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)

        async_session_factory = async_sessionmaker(
            bind=engine,
            class_=AsyncSession,
            expire_on_commit=False,
        )

        async with async_session_factory() as session:
            yield session
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="frame_storage")
async def frame_storage_fixture() -> AsyncIterator[InMemoryStorage]:
    """Create a fresh InMemoryStorage instance for each test."""
    storage = InMemoryStorage()
    yield storage
    await storage.close()


@pytest_asyncio.fixture(name="frame_client")
async def frame_client_fixture(
    frame_session: AsyncSession, frame_storage: InMemoryStorage
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with session and storage dependencies overridden."""
    from contextlib import asynccontextmanager

    from zndraw_auth import get_session
    from zndraw_auth.settings import AuthSettings
    from zndraw_joblib.settings import JobLibSettings

    from zndraw.app import app
    from zndraw.dependencies import (
        get_joblib_settings,
        get_redis,
        get_result_backend,
        get_storage,
        get_tsio,
    )

    mock_sio = MockSioServer()

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield frame_session

    @asynccontextmanager
    async def test_session_maker():
        yield frame_session

    def get_storage_override() -> InMemoryStorage:
        return frame_storage

    def get_sio_override() -> MockSioServer:
        return mock_sio

    # Mock Redis for WritableRoomDep (returns None = no edit lock)
    mock_redis = AsyncMock()
    mock_redis.get = AsyncMock(return_value=None)

    app.state.auth_settings = AuthSettings()
    app.state.session_maker = test_session_maker
    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_storage] = get_storage_override
    app.dependency_overrides[get_tsio] = get_sio_override
    app.dependency_overrides[get_redis] = lambda: mock_redis
    app.dependency_overrides[get_result_backend] = lambda: AsyncMock()
    app.dependency_overrides[get_joblib_settings] = lambda: JobLibSettings()

    async with AsyncClient(
        transport=ASGITransport(app=app),
        base_url="http://test",
    ) as client:
        yield client

    app.dependency_overrides.clear()


_create_user = create_test_user_in_db
_create_room = create_test_room
_auth_header = auth_header


# =============================================================================
# List Frames Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_frames_empty_room(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test listing frames from an empty room returns empty list."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/x-msgpack"

    frames = decode_msgpack_response(response.content)
    assert frames == []


@pytest.mark.asyncio
async def test_list_frames_with_data(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test listing frames with data returns all frames."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add frames to storage
    await frame_storage.extend(room.id, [{"a": 1}, {"b": 2}, {"c": 3}])  # type: ignore

    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 3
    # Verify first frame content
    assert raw_frame_to_dict(frames[0]) == {"a": 1}


@pytest.mark.asyncio
async def test_list_frames_with_range(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test listing frames with range query params."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add frames to storage
    await frame_storage.extend(room.id, [{"a": 1}, {"b": 2}, {"c": 3}, {"d": 4}])  # type: ignore

    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames?start=1&stop=3",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 2
    assert raw_frame_to_dict(frames[0]) == {"b": 2}
    assert raw_frame_to_dict(frames[1]) == {"c": 3}


@pytest.mark.asyncio
async def test_list_frames_room_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test listing frames from non-existent room returns 404."""
    _, token = await _create_user(frame_session)

    response = await frame_client.get(
        "/v1/rooms/99999/frames",
        headers=_auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


@pytest.mark.asyncio
async def test_list_frames_with_indices(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test listing specific frames by indices parameter."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add 5 frames
    await frame_storage.extend(
        room.id,
        [{"a": 0}, {"b": 1}, {"c": 2}, {"d": 3}, {"e": 4}],  # type: ignore
    )

    # Request specific indices
    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames?indices=1,3",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 2
    assert raw_frame_to_dict(frames[0]) == {"b": 1}
    assert raw_frame_to_dict(frames[1]) == {"d": 3}


@pytest.mark.asyncio
async def test_list_frames_with_keys_filter(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test listing frames with keys parameter to filter frame data."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add frames with multiple keys
    await frame_storage.extend(
        room.id,
        [
            {"x": 1, "y": 2, "z": 3},
            {"x": 4, "y": 5, "z": 6},
        ],
    )  # type: ignore

    # Request only x and z keys
    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames?keys=x,z",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 2
    assert raw_frame_to_dict(frames[0]) == {"x": 1, "z": 3}
    assert raw_frame_to_dict(frames[1]) == {"x": 4, "z": 6}


@pytest.mark.asyncio
async def test_list_frames_with_indices_and_keys(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test listing specific indices with keys filter."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add frames
    await frame_storage.extend(
        room.id,
        [
            {"a": 1, "b": 2},
            {"a": 3, "b": 4},
            {"a": 5, "b": 6},
        ],
    )  # type: ignore

    # Request index 2 with only key 'a'
    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames?indices=2&keys=a",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert raw_frame_to_dict(frames[0]) == {"a": 5}


# =============================================================================
# Get Single Frame Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_frame(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test getting a single frame by index."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add frames to storage
    await frame_storage.extend(room.id, [{"a": 1}, {"b": 2}])  # type: ignore

    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames/1",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/x-msgpack"

    # GET frame returns a list with single element
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert raw_frame_to_dict(frames[0]) == {"b": 2}


@pytest.mark.asyncio
async def test_get_frame_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test getting non-existent frame returns 404."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames/99",
        headers=_auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


@pytest.mark.asyncio
async def test_get_frame_room_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test getting frame from non-existent room returns 404."""
    _, token = await _create_user(frame_session)

    response = await frame_client.get(
        "/v1/rooms/99999/frames/0",
        headers=_auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


# =============================================================================
# Frame Metadata Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_frame_metadata(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test getting metadata for a frame with mixed scalar and array data."""
    from ase import Atoms
    from asebytes import encode

    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Create an Atoms object with calc results
    atoms = Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    atoms.info["energy"] = -42.5

    raw = encode(atoms)
    await frame_storage.extend(room.id, [raw])  # type: ignore

    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames/0/metadata",
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    data = response.json()
    assert data["frame_id"] == 0
    assert data["source_room"] == room.id

    meta = data["metadata"]
    # Check that key frame keys exist
    assert "arrays.positions" in meta
    assert meta["arrays.positions"]["type"] == "array"
    assert meta["arrays.positions"]["shape"] == [3, 3]  # 3 atoms, 3 coords

    assert "arrays.numbers" in meta
    assert meta["arrays.numbers"]["type"] == "array"
    assert meta["arrays.numbers"]["shape"] == [3]

    assert "info.energy" in meta
    assert meta["info.energy"]["type"] == "scalar"
    assert meta["info.energy"]["dtype"] == "float64"


@pytest.mark.asyncio
async def test_get_frame_metadata_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test getting metadata for non-existent frame returns 404."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    response = await frame_client.get(
        f"/v1/rooms/{room.id}/frames/99/metadata",
        headers=_auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


@pytest.mark.asyncio
async def test_get_frame_metadata_room_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test getting metadata for a frame in a non-existent room returns 404."""
    _, token = await _create_user(frame_session)

    response = await frame_client.get(
        "/v1/rooms/nonexistent-room/frames/0/metadata",
        headers=_auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


# =============================================================================
# Append Frames Tests
# =============================================================================


@pytest.mark.asyncio
async def test_append_frames(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test appending frames to storage."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    frame_a = _make_json_frame("H2")
    frame_b = _make_json_frame("H2O")

    response = await frame_client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [frame_a, frame_b]},
        headers=_auth_header(token),
    )
    assert response.status_code == 201

    result = FrameBulkResponse.model_validate(response.json())
    assert len(result.frames) == 2
    assert result.total == 2
    assert result.start == 0
    assert result.stop == 2


@pytest.mark.asyncio
async def test_append_frames_multiple_times(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test appending frames multiple times."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # First append
    response = await frame_client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [_make_json_frame("H2")]},
        headers=_auth_header(token),
    )
    assert response.status_code == 201
    result = FrameBulkResponse.model_validate(response.json())
    assert result.total == 1
    assert result.start == 0
    assert result.stop == 1

    # Second append
    response = await frame_client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [_make_json_frame("H2O"), _make_json_frame("CH4")]},
        headers=_auth_header(token),
    )
    assert response.status_code == 201
    result = FrameBulkResponse.model_validate(response.json())
    assert result.total == 3
    assert result.start == 1
    assert result.stop == 3


@pytest.mark.asyncio
async def test_append_frames_room_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test appending frames to non-existent room returns 404."""
    _, token = await _create_user(frame_session)

    response = await frame_client.post(
        "/v1/rooms/99999/frames",
        json={"frames": [_make_json_frame("H2")]},
        headers=_auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


@pytest.mark.asyncio
async def test_append_frames_empty_list_rejected(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test appending empty frames list is rejected (422)."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    response = await frame_client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": []},
        headers=_auth_header(token),
    )
    # FastAPI returns 422 for validation errors
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_append_frames_exceeds_max_length(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test appending more than 1000 frames is rejected (422)."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    frame = _make_json_frame("H2")
    response = await frame_client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [frame] * 1001},
        headers=_auth_header(token),
    )
    assert response.status_code == 422


# =============================================================================
# Update Frame Tests
# =============================================================================


@pytest.mark.asyncio
async def test_update_frame(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test updating a frame at specific index."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add frames to storage
    await frame_storage.extend(room.id, [{"a": 1}, {"b": 2}])  # type: ignore

    new_frame = _make_json_frame("He")

    response = await frame_client.put(
        f"/v1/rooms/{room.id}/frames/1",
        json={"data": new_frame},
        headers=_auth_header(token),
    )
    assert response.status_code == 200

    result = FrameResponse.model_validate(response.json())
    assert result.index == 1


@pytest.mark.asyncio
async def test_update_frame_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test updating non-existent frame returns 404."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    response = await frame_client.put(
        f"/v1/rooms/{room.id}/frames/99",
        json={"data": _make_json_frame("H2")},
        headers=_auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


# =============================================================================
# Merge Frame (PATCH) Tests
# =============================================================================


@pytest.mark.asyncio
async def test_merge_frame(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test partial update merges new keys into existing frame."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    await frame_storage.extend(room.id, [{"a": 1, "b": 2}])  # type: ignore

    # Send PATCH with msgpack body updating key "a" and adding key "c"
    patch_data = msgpack.packb({"a": 99, "c": 3})
    response = await frame_client.patch(
        f"/v1/rooms/{room.id}/frames/0",
        content=patch_data,
        headers={**_auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 200

    result = FrameMergeResponse.model_validate(response.json())
    assert result.index == 0
    assert set(result.updated_keys) == {"a", "c"}

    # Verify merged data in storage
    stored = await frame_storage.get(room.id, 0)  # type: ignore
    assert stored is not None
    assert raw_frame_to_dict(stored) == {"a": 99, "b": 2, "c": 3}


@pytest.mark.asyncio
async def test_merge_frame_preserves_untouched_keys(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test partial update does not remove keys not in the patch."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    await frame_storage.extend(room.id, [{"x": 10, "y": 20, "z": 30}])  # type: ignore

    # Only update "y"
    patch_data = msgpack.packb({"y": 99})
    response = await frame_client.patch(
        f"/v1/rooms/{room.id}/frames/0",
        content=patch_data,
        headers={**_auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 200

    stored = await frame_storage.get(room.id, 0)  # type: ignore
    assert stored is not None
    assert raw_frame_to_dict(stored) == {"x": 10, "y": 99, "z": 30}


@pytest.mark.asyncio
async def test_merge_frame_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test merging non-existent frame returns 404."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    patch_data = msgpack.packb({"a": 1})
    response = await frame_client.patch(
        f"/v1/rooms/{room.id}/frames/99",
        content=patch_data,
        headers={**_auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


@pytest.mark.asyncio
async def test_merge_frame_room_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test merging frame in non-existent room returns 404."""
    _, token = await _create_user(frame_session)

    patch_data = msgpack.packb({"a": 1})
    response = await frame_client.patch(
        "/v1/rooms/99999/frames/0",
        content=patch_data,
        headers={**_auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


@pytest.mark.asyncio
async def test_merge_frame_preserves_msgpack_str_type(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test that PATCH preserves msgpack str/bin distinction for numpy arrays.

    The frontend's msgpack-numpy format uses msgpack str type for dtype
    (e.g., "type": "<f4"). If the backend re-packs these as bin type,
    the frontend decoder fails because it checks `typeof typeStr === "string"`.

    This reproduces the bug: after PATCH, reading the frame back should
    produce identical msgpack encoding to what the frontend originally sent.
    """
    import struct

    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Create initial frame with a numpy-format position array (float64)
    # This mimics what asebytes.encode produces
    original_positions = struct.pack("<9d", *range(9))  # 9 float64 values
    original_numpy = msgpack.packb(
        {"nd": True, "type": "<f8", "shape": [3, 3], "data": original_positions},
        use_bin_type=True,
    )
    await frame_storage.extend(
        room.id,
        [  # type: ignore[arg-type]  # raw bytes frames bypass extend's str API
            {
                b"arrays.positions": original_numpy,
                b"arrays.numbers": msgpack.packb([1, 1, 1]),
            }
        ],
    )

    # PATCH with float32 positions (what the frontend sends after editing)
    # Frontend's packBinary produces: {str"arrays.positions": {str"nd": true, str"type": str"<f4", ...}}
    edited_positions = struct.pack("<9f", *[float(x) for x in range(9)])
    patch_body = msgpack.packb(
        {
            "arrays.positions": {
                "nd": True,
                "type": "<f4",
                "shape": [3, 3],
                "data": edited_positions,
            }
        },
        use_bin_type=True,
    )
    response = await frame_client.patch(
        f"/v1/rooms/{room.id}/frames/0",
        content=patch_body,
        headers={**_auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 200

    # Read the merged frame back and decode the stored positions value
    stored_frame = await frame_storage.get(room.id, 0)  # type: ignore
    assert stored_frame is not None

    # The stored value must preserve msgpack str/bin distinction.
    # Frontend decodes dtype as JS string only if it's msgpack str type.
    # With raw=True re-packing, strings become bytes (msgpack bin) — breaking the frontend.
    stored_positions_bytes = stored_frame[b"arrays.positions"]

    # Decode with raw=False: msgpack str → Python str, msgpack bin → Python bytes
    decoded = msgpack.unpackb(stored_positions_bytes, raw=False)

    # Keys must be strings (msgpack str type), not bytes (msgpack bin type)
    assert "type" in decoded, (
        f"Inner keys must be msgpack str type, got bytes keys: {list(decoded.keys())!r}"
    )
    # dtype value must be a string, not bytes
    assert isinstance(decoded["type"], str), (
        f"dtype field must be msgpack str type for frontend compatibility, "
        f"got {type(decoded['type']).__name__}: {decoded['type']!r}"
    )
    assert decoded["type"] == "<f4"
    assert decoded["nd"] is True
    assert decoded["shape"] == [3, 3]

    # Verify untouched keys are preserved
    assert b"arrays.numbers" in stored_frame


# =============================================================================
# Delete Frame Tests
# =============================================================================


@pytest.mark.asyncio
async def test_delete_frame(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """Test deleting a frame at specific index."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add frames to storage
    await frame_storage.extend(room.id, [{"a": 1}, {"b": 2}, {"c": 3}])  # type: ignore

    response = await frame_client.delete(
        f"/v1/rooms/{room.id}/frames/1",
        headers=_auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify frame was deleted and others shifted
    assert await frame_storage.get_length(room.id) == 2  # type: ignore
    assert await frame_storage.get(room.id, 0) == make_raw_frame({"a": 1})  # type: ignore
    assert await frame_storage.get(room.id, 1) == make_raw_frame({"c": 3})  # type: ignore


@pytest.mark.asyncio
async def test_delete_frame_not_found(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test deleting non-existent frame returns 404."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    response = await frame_client.delete(
        f"/v1/rooms/{room.id}/frames/99",
        headers=_auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_frames_require_authentication(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """Test that all frame endpoints require authentication."""
    user, _ = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # All endpoints should return 401 without auth
    endpoints = [
        ("GET", f"/v1/rooms/{room.id}/frames"),
        ("GET", f"/v1/rooms/{room.id}/frames/0"),
        ("POST", f"/v1/rooms/{room.id}/frames"),
        ("PUT", f"/v1/rooms/{room.id}/frames/0"),
        ("PATCH", f"/v1/rooms/{room.id}/frames/0"),
        ("DELETE", f"/v1/rooms/{room.id}/frames/0"),
    ]

    for method, url in endpoints:
        if method == "GET":
            response = await frame_client.get(url)
        elif method == "POST":
            response = await frame_client.post(
                url, json={"frames": [_make_json_frame("H2")]}
            )
        elif method == "PUT":
            response = await frame_client.put(
                url, json={"data": _make_json_frame("H2")}
            )
        elif method == "PATCH":
            response = await frame_client.patch(
                url,
                content=msgpack.packb({"a": 1}),
                headers={"Content-Type": "application/msgpack"},
            )
        else:  # DELETE
            response = await frame_client.delete(url)

        assert response.status_code == 401, f"{method} {url} should require auth"


# =============================================================================
# Frame Validation Tests
# =============================================================================


def _make_bare_json_frame(formula: str = "H2") -> dict[str, Any]:
    """Create a JSON-serializable frame dict WITHOUT colors/radii."""
    from asebytes import encode

    atoms = ase.Atoms(
        formula,
        positions=[
            [i, 0, 0] for i in range(ase.Atoms(formula).get_global_number_of_atoms())
        ],
    )
    # Encode without enrichment → no arrays.colors / arrays.radii
    encoded = encode(atoms)
    result: dict[str, Any] = {}
    for key, value in encoded.items():
        import base64

        key_str = "b64:" + base64.b64encode(key).decode("ascii")
        value_str = base64.b64encode(value).decode("ascii")
        result[key_str] = value_str
    return result


@pytest.mark.asyncio
async def test_append_rejects_frames_without_colors_radii(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """POST /frames rejects frames missing arrays.colors and arrays.radii."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    bare_frame = _make_bare_json_frame("H2")

    response = await frame_client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [bare_frame]},
        headers=_auth_header(token),
    )
    assert response.status_code == 422

    problem = ProblemDetail.model_validate(response.json())
    assert "colors" in (problem.detail or "").lower()
    assert "radii" in (problem.detail or "").lower()


@pytest.mark.asyncio
async def test_update_rejects_frame_without_colors_radii(
    frame_client: AsyncClient,
    frame_session: AsyncSession,
    frame_storage: InMemoryStorage,
) -> None:
    """PUT /frames/{index} rejects frames missing arrays.colors and arrays.radii."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    # Add a valid frame so index 0 exists
    await frame_storage.extend(room.id, [{"a": 1}])  # type: ignore

    bare_frame = _make_bare_json_frame("H2")

    response = await frame_client.put(
        f"/v1/rooms/{room.id}/frames/0",
        json={"data": bare_frame},
        headers=_auth_header(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_append_accepts_enriched_frames(
    frame_client: AsyncClient, frame_session: AsyncSession
) -> None:
    """POST /frames accepts frames that already have colors and radii."""
    user, token = await _create_user(frame_session)
    room = await _create_room(frame_session, user)

    enriched_frame = _make_json_frame("H2")

    response = await frame_client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [enriched_frame]},
        headers=_auth_header(token),
    )
    assert response.status_code == 201
