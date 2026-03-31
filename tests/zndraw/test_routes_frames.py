"""Tests for Frame REST API endpoints."""

from typing import Any

import ase
import msgpack
import pytest
from helpers import (
    auth_header,
    create_test_room,
    create_test_user_in_db,
    decode_msgpack_response,
    make_raw_frame,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.client import atoms_to_json_dict
from zndraw.exceptions import FrameNotFound, ProblemDetail, RoomNotFound
from zndraw.schemas import (
    FrameBulkResponse,
    FrameMergeResponse,
    FrameResponse,
    StatusResponse,
)
from zndraw.storage import FrameStorage, RawFrame


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
# List Frames Tests
# =============================================================================


@pytest.mark.asyncio
async def test_list_frames_empty_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test listing frames from an empty room returns empty list."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/frames",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/x-msgpack"

    frames = decode_msgpack_response(response.content)
    assert frames == []


@pytest.mark.asyncio
async def test_list_frames_with_data(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test listing frames with data returns all frames."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add frames to storage
    await frame_storage[room.id].extend(
        [make_raw_frame({"a": 1}), make_raw_frame({"b": 2}), make_raw_frame({"c": 3})]
    )
    response = await client.get(
        f"/v1/rooms/{room.id}/frames",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 3
    # Verify first frame content
    assert raw_frame_to_dict(frames[0]) == {"a": 1}


@pytest.mark.asyncio
async def test_list_frames_with_range(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test listing frames with range query params."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add frames to storage
    await frame_storage[room.id].extend(
        [
            make_raw_frame({"a": 1}),
            make_raw_frame({"b": 2}),
            make_raw_frame({"c": 3}),
            make_raw_frame({"d": 4}),
        ]
    )
    response = await client.get(
        f"/v1/rooms/{room.id}/frames?start=1&stop=3",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 2
    assert raw_frame_to_dict(frames[0]) == {"b": 2}
    assert raw_frame_to_dict(frames[1]) == {"c": 3}


@pytest.mark.asyncio
async def test_list_frames_room_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test listing frames from non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/99999/frames",
        headers=auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


@pytest.mark.asyncio
async def test_list_frames_with_indices(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test listing specific frames by indices parameter."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add 5 frames
    await frame_storage[room.id].extend(
        [
            make_raw_frame({"a": 0}),
            make_raw_frame({"b": 1}),
            make_raw_frame({"c": 2}),
            make_raw_frame({"d": 3}),
            make_raw_frame({"e": 4}),
        ]
    )

    # Request specific indices
    response = await client.get(
        f"/v1/rooms/{room.id}/frames?indices=1,3",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 2
    assert raw_frame_to_dict(frames[0]) == {"b": 1}
    assert raw_frame_to_dict(frames[1]) == {"d": 3}


@pytest.mark.asyncio
async def test_list_frames_with_keys_filter(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test listing frames with keys parameter to filter frame data."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add frames with multiple keys
    await frame_storage[room.id].extend(
        [
            make_raw_frame({"x": 1, "y": 2, "z": 3}),
            make_raw_frame({"x": 4, "y": 5, "z": 6}),
        ]
    )
    # Request only x and z keys
    response = await client.get(
        f"/v1/rooms/{room.id}/frames?keys=x,z",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    frames = decode_msgpack_response(response.content)
    assert len(frames) == 2
    assert raw_frame_to_dict(frames[0]) == {"x": 1, "z": 3}
    assert raw_frame_to_dict(frames[1]) == {"x": 4, "z": 6}


@pytest.mark.asyncio
async def test_list_frames_with_indices_and_keys(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test listing specific indices with keys filter."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add frames
    await frame_storage[room.id].extend(
        [
            make_raw_frame({"a": 1, "b": 2}),
            make_raw_frame({"a": 3, "b": 4}),
            make_raw_frame({"a": 5, "b": 6}),
        ]
    )
    # Request index 2 with only key 'a'
    response = await client.get(
        f"/v1/rooms/{room.id}/frames?indices=2&keys=a",
        headers=auth_header(token),
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
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test getting a single frame by index."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add frames to storage
    await frame_storage[room.id].extend(
        [make_raw_frame({"a": 1}), make_raw_frame({"b": 2})]
    )
    response = await client.get(
        f"/v1/rooms/{room.id}/frames/1",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/x-msgpack"

    # GET frame returns a list with single element
    frames = decode_msgpack_response(response.content)
    assert len(frames) == 1
    assert raw_frame_to_dict(frames[0]) == {"b": 2}


@pytest.mark.asyncio
async def test_get_frame_not_found(client: AsyncClient, session: AsyncSession) -> None:
    """Test getting non-existent frame returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/frames/99",
        headers=auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


@pytest.mark.asyncio
async def test_get_frame_room_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test getting frame from non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/99999/frames/0",
        headers=auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


# =============================================================================
# Frame Metadata Tests
# =============================================================================


@pytest.mark.asyncio
async def test_get_frame_metadata(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test getting metadata for a frame with mixed scalar and array data."""
    from ase import Atoms
    from asebytes import encode

    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Create an Atoms object with calc results
    atoms = Atoms("H2O", positions=[[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    atoms.info["energy"] = -42.5

    raw = encode(atoms)
    await frame_storage[room.id].extend([raw])
    response = await client.get(
        f"/v1/rooms/{room.id}/frames/0/metadata",
        headers=auth_header(token),
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
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test getting metadata for non-existent frame returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/frames/99/metadata",
        headers=auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


@pytest.mark.asyncio
async def test_get_frame_metadata_room_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test getting metadata for a frame in a non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.get(
        "/v1/rooms/nonexistent-room/frames/0/metadata",
        headers=auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


# =============================================================================
# Append Frames Tests
# =============================================================================


@pytest.mark.asyncio
async def test_append_frames(client: AsyncClient, session: AsyncSession) -> None:
    """Test appending frames to storage."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    frame_a = _make_json_frame("H2")
    frame_b = _make_json_frame("H2O")

    response = await client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [frame_a, frame_b]},
        headers=auth_header(token),
    )
    assert response.status_code == 201

    result = FrameBulkResponse.model_validate(response.json())
    assert len(result.frames) == 2
    assert result.total == 2
    assert result.start == 0
    assert result.stop == 2


@pytest.mark.asyncio
async def test_append_frames_multiple_times(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test appending frames multiple times."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # First append
    response = await client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [_make_json_frame("H2")]},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    result = FrameBulkResponse.model_validate(response.json())
    assert result.total == 1
    assert result.start == 0
    assert result.stop == 1

    # Second append
    response = await client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [_make_json_frame("H2O"), _make_json_frame("CH4")]},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    result = FrameBulkResponse.model_validate(response.json())
    assert result.total == 3
    assert result.start == 1
    assert result.stop == 3


@pytest.mark.asyncio
async def test_append_frames_room_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test appending frames to non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    response = await client.post(
        "/v1/rooms/99999/frames",
        json={"frames": [_make_json_frame("H2")]},
        headers=auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


@pytest.mark.asyncio
async def test_append_frames_empty_list_rejected(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test appending empty frames list is rejected (422)."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": []},
        headers=auth_header(token),
    )
    # FastAPI returns 422 for validation errors
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_append_frames_exceeds_max_length(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test appending more than 1000 frames is rejected (422)."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    frame = _make_json_frame("H2")
    response = await client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [frame] * 1001},
        headers=auth_header(token),
    )
    assert response.status_code == 422


# =============================================================================
# Update Frame Tests
# =============================================================================


@pytest.mark.asyncio
async def test_update_frame(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test updating a frame at specific index."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add frames to storage
    await frame_storage[room.id].extend(
        [make_raw_frame({"a": 1}), make_raw_frame({"b": 2})]
    )
    new_frame = _make_json_frame("He")

    response = await client.put(
        f"/v1/rooms/{room.id}/frames/1",
        json={"data": new_frame},
        headers=auth_header(token),
    )
    assert response.status_code == 200

    result = FrameResponse.model_validate(response.json())
    assert result.index == 1


@pytest.mark.asyncio
async def test_update_frame_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test updating non-existent frame returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/rooms/{room.id}/frames/99",
        json={"data": _make_json_frame("H2")},
        headers=auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


# =============================================================================
# Merge Frame (PATCH) Tests
# =============================================================================


@pytest.mark.asyncio
async def test_merge_frame(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test partial update merges new keys into existing frame."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await frame_storage[room.id].extend([make_raw_frame({"a": 1, "b": 2})])
    # Send PATCH with msgpack body updating key "a" and adding key "c"
    patch_data = msgpack.packb({"a": 99, "c": 3})
    response = await client.patch(
        f"/v1/rooms/{room.id}/frames/0",
        content=patch_data,
        headers={**auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 200

    result = FrameMergeResponse.model_validate(response.json())
    assert result.index == 0
    assert set(result.updated_keys) == {"a", "c"}

    # Verify merged data in storage
    stored = await frame_storage[room.id][0]
    assert stored is not None
    assert raw_frame_to_dict(stored) == {"a": 99, "b": 2, "c": 3}


@pytest.mark.asyncio
async def test_merge_frame_preserves_untouched_keys(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test partial update does not remove keys not in the patch."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await frame_storage[room.id].extend([make_raw_frame({"x": 10, "y": 20, "z": 30})])
    # Only update "y"
    patch_data = msgpack.packb({"y": 99})
    response = await client.patch(
        f"/v1/rooms/{room.id}/frames/0",
        content=patch_data,
        headers={**auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 200

    stored = await frame_storage[room.id][0]
    assert stored is not None
    assert raw_frame_to_dict(stored) == {"x": 10, "y": 99, "z": 30}


@pytest.mark.asyncio
async def test_merge_frame_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test merging non-existent frame returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    patch_data = msgpack.packb({"a": 1})
    response = await client.patch(
        f"/v1/rooms/{room.id}/frames/99",
        content=patch_data,
        headers={**auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


@pytest.mark.asyncio
async def test_merge_frame_room_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test merging frame in non-existent room returns 404."""
    _, token = await create_test_user_in_db(session)

    patch_data = msgpack.packb({"a": 1})
    response = await client.patch(
        "/v1/rooms/99999/frames/0",
        content=patch_data,
        headers={**auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomNotFound.type_uri()


@pytest.mark.asyncio
async def test_merge_frame_preserves_msgpack_str_type(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test that PATCH preserves msgpack str/bin distinction for numpy arrays.

    The frontend's msgpack-numpy format uses msgpack str type for dtype
    (e.g., "type": "<f4"). If the backend re-packs these as bin type,
    the frontend decoder fails because it checks `typeof typeStr === "string"`.

    This reproduces the bug: after PATCH, reading the frame back should
    produce identical msgpack encoding to what the frontend originally sent.
    """
    import struct

    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Create initial frame with a numpy-format position array (float64)
    # This mimics what asebytes.encode produces
    original_positions = struct.pack("<9d", *range(9))  # 9 float64 values
    original_numpy = msgpack.packb(
        {"nd": True, "type": "<f8", "shape": [3, 3], "data": original_positions},
        use_bin_type=True,
    )
    await frame_storage[room.id].extend(
        [
            {
                b"arrays.positions": original_numpy,
                b"arrays.numbers": msgpack.packb([1, 1, 1]),
            }
        ]
    )

    # PATCH with float32 positions (what the frontend sends after editing)
    # Frontend's packBinary produces:
    # {str"arrays.positions": {str"nd": true, str"type": str"<f4", ...}}
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
    response = await client.patch(
        f"/v1/rooms/{room.id}/frames/0",
        content=patch_body,
        headers={**auth_header(token), "Content-Type": "application/msgpack"},
    )
    assert response.status_code == 200

    # Read the merged frame back and decode the stored positions value
    stored_frame = await frame_storage[room.id][0]
    assert stored_frame is not None

    # The stored value must preserve msgpack str/bin distinction.
    # Frontend decodes dtype as JS string only if it's msgpack str type.
    # With raw=True re-packing, strings become bytes
    # (msgpack bin) — breaking the frontend.
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
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test deleting a frame at specific index."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add frames to storage
    await frame_storage[room.id].extend(
        [make_raw_frame({"a": 1}), make_raw_frame({"b": 2}), make_raw_frame({"c": 3})]
    )
    response = await client.delete(
        f"/v1/rooms/{room.id}/frames/1",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())

    # Verify frame was deleted and others shifted
    assert await frame_storage.get_length(room.id) == 2
    assert await frame_storage[room.id][0] == make_raw_frame({"a": 1})
    assert await frame_storage[room.id][1] == make_raw_frame({"c": 3})


@pytest.mark.asyncio
async def test_delete_frame_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test deleting non-existent frame returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.delete(
        f"/v1/rooms/{room.id}/frames/99",
        headers=auth_header(token),
    )
    assert response.status_code == 404

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == FrameNotFound.type_uri()


# =============================================================================
# Authentication Tests
# =============================================================================


@pytest.mark.asyncio
async def test_frames_require_authentication(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Test that all frame endpoints require authentication."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

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
            response = await client.get(url)
        elif method == "POST":
            response = await client.post(url, json={"frames": [_make_json_frame("H2")]})
        elif method == "PUT":
            response = await client.put(url, json={"data": _make_json_frame("H2")})
        elif method == "PATCH":
            response = await client.patch(
                url,
                content=msgpack.packb({"a": 1}),
                headers={"Content-Type": "application/msgpack"},
            )
        else:  # DELETE
            response = await client.delete(url)

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
    client: AsyncClient, session: AsyncSession
) -> None:
    """POST /frames rejects frames missing arrays.colors and arrays.radii."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    bare_frame = _make_bare_json_frame("H2")

    response = await client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [bare_frame]},
        headers=auth_header(token),
    )
    assert response.status_code == 422

    problem = ProblemDetail.model_validate(response.json())
    assert "colors" in (problem.detail or "").lower()
    assert "radii" in (problem.detail or "").lower()


@pytest.mark.asyncio
async def test_update_rejects_frame_without_colors_radii(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """PUT /frames/{index} rejects frames missing arrays.colors and arrays.radii."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Add a valid frame so index 0 exists
    await frame_storage[room.id].extend([make_raw_frame({"a": 1})])
    bare_frame = _make_bare_json_frame("H2")

    response = await client.put(
        f"/v1/rooms/{room.id}/frames/0",
        json={"data": bare_frame},
        headers=auth_header(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_append_accepts_enriched_frames(
    client: AsyncClient, session: AsyncSession
) -> None:
    """POST /frames accepts frames that already have colors and radii."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    enriched_frame = _make_json_frame("H2")

    response = await client.post(
        f"/v1/rooms/{room.id}/frames",
        json={"frames": [enriched_frame]},
        headers=auth_header(token),
    )
    assert response.status_code == 201
