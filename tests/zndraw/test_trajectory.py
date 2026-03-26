"""Tests for Trajectory REST API endpoints (download/upload + download tokens)."""

import io
from typing import Any

import ase
import ase.io
import numpy as np
import pytest
from asebytes import decode, encode
from helpers import auth_header, create_test_room, create_test_user_in_db
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.exceptions import InvalidPayload, ProblemDetail, RoomLocked, RoomReadOnly
from zndraw.schemas import FrameBulkResponse
from zndraw.storage import FrameStorage

# =============================================================================
# Test Helpers
# =============================================================================


def _make_atoms(
    symbols: str = "H2",
    positions: list[list[float]] | None = None,
    info: dict[str, Any] | None = None,
) -> ase.Atoms:
    """Create an ase.Atoms object with optional positions and info."""
    if positions is None:
        positions = [[0, 0, 0], [1, 0, 0]]
    atoms = ase.Atoms(symbols, positions=positions)
    if info:
        atoms.info.update(info)
    return atoms


def _atoms_to_file_bytes(atoms_list: list[ase.Atoms], fmt: str = "extxyz") -> bytes:
    """Write a list of Atoms to a bytes buffer in the given format."""
    buf = io.StringIO()
    ase.io.write(buf, atoms_list, format=fmt)
    return buf.getvalue().encode("utf-8")


def _parse_trajectory(text: str, fmt: str = "extxyz") -> list[ase.Atoms]:
    """Parse trajectory text into a list of Atoms (with correct type for pyright)."""
    result = ase.io.read(io.StringIO(text), index=":", format=fmt)
    if not isinstance(result, list):
        return [result]  # type: ignore[list-item]
    return result  # type: ignore[return-value]


async def _add_atoms_to_storage(
    storage: FrameStorage, room_id: str, atoms_list: list[ase.Atoms]
) -> None:
    """Encode and store a list of Atoms objects in storage."""
    frames = [encode(a) for a in atoms_list]
    await storage[room_id].extend(frames)


# =============================================================================
# Download Tests
# =============================================================================


@pytest.mark.asyncio
async def test_download_single_frame(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test downloading a single frame by index."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms("H2", [[0, 0, 0], [1, 0, 0]])
    await _add_atoms_to_storage(frame_storage, room.id, [atoms])

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?indices=0",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    parsed = _parse_trajectory(response.text)
    assert len(parsed) == 1
    assert len(parsed[0]) == 2


@pytest.mark.asyncio
async def test_download_all_frames(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test downloading all frames when no indices specified."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms_list = [
        _make_atoms("H2", [[0, 0, 0], [1, 0, 0]]),
        _make_atoms("H3", [[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
        _make_atoms("H", [[0, 0, 0]]),
    ]
    await _add_atoms_to_storage(frame_storage, room.id, atoms_list)

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    parsed = _parse_trajectory(response.text)
    assert len(parsed) == 3
    assert len(parsed[0]) == 2
    assert len(parsed[1]) == 3
    assert len(parsed[2]) == 1


@pytest.mark.asyncio
async def test_download_specific_indices(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test downloading specific frame indices."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms_list = [
        _make_atoms("H", [[0, 0, 0]]),
        _make_atoms("H2", [[0, 0, 0], [1, 0, 0]]),
        _make_atoms("H3", [[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
        _make_atoms("He2", [[0, 0, 0], [5, 0, 0]]),
        _make_atoms("O", [[3, 0, 0]]),
    ]
    await _add_atoms_to_storage(frame_storage, room.id, atoms_list)

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?indices=0,2,4",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    parsed = _parse_trajectory(response.text)
    assert len(parsed) == 3
    # Frame 0: H (1 atom), Frame 2: H3 (3 atoms), Frame 4: O (1 atom)
    assert len(parsed[0]) == 1
    assert len(parsed[1]) == 3
    assert len(parsed[2]) == 1


@pytest.mark.asyncio
async def test_download_with_atom_selection(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test downloading with atom selection filter."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms("H2O", [[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    await _add_atoms_to_storage(frame_storage, room.id, [atoms])

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?selection=0,2",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    parsed = _parse_trajectory(response.text)
    assert len(parsed) == 1
    assert len(parsed[0]) == 2
    # Atom 0 is H, atom 2 is O
    np.testing.assert_array_equal(parsed[0].get_atomic_numbers(), [1, 8])


@pytest.mark.asyncio
async def test_download_preserves_info(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test that atoms.info is preserved through download."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms("H2", [[0, 0, 0], [1, 0, 0]], info={"key": "value"})
    await _add_atoms_to_storage(frame_storage, room.id, [atoms])

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    parsed = _parse_trajectory(response.text)
    assert parsed[0].info["key"] == "value"


@pytest.mark.asyncio
async def test_download_empty_room(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test downloading from empty room returns 400."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory",
        headers=auth_header(token),
    )
    assert response.status_code == 400

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == InvalidPayload.type_uri()


@pytest.mark.asyncio
async def test_download_invalid_index(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test downloading with out-of-range index returns 400."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms_list = [_make_atoms() for _ in range(3)]
    await _add_atoms_to_storage(frame_storage, room.id, atoms_list)

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?indices=99",
        headers=auth_header(token),
    )
    assert response.status_code == 400

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == InvalidPayload.type_uri()


@pytest.mark.asyncio
async def test_download_custom_filename(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test that custom filename appears in Content-Disposition header."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms()
    await _add_atoms_to_storage(frame_storage, room.id, [atoms])

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?filename=my_traj.extxyz",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert 'filename="my_traj.extxyz"' in response.headers["content-disposition"]


@pytest.mark.asyncio
@pytest.mark.parametrize("fmt", ["extxyz", "xyz"])
async def test_download_formats(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
    fmt: str,
) -> None:
    """Test downloading in different supported formats."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms("H2", [[0, 0, 0], [1, 0, 0]])
    await _add_atoms_to_storage(frame_storage, room.id, [atoms])

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?format={fmt}",
        headers=auth_header(token),
    )
    assert response.status_code == 200

    parsed = _parse_trajectory(response.text, fmt)
    assert len(parsed) == 1
    assert len(parsed[0]) == 2


@pytest.mark.asyncio
async def test_download_requires_auth(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test downloading without authentication returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(f"/v1/rooms/{room.id}/trajectory")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_download_unsupported_format(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test downloading with unsupported format returns 400."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms()
    await _add_atoms_to_storage(frame_storage, room.id, [atoms])

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?format=invalid",
        headers=auth_header(token),
    )
    assert response.status_code == 400

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == InvalidPayload.type_uri()


# =============================================================================
# Upload Tests
# =============================================================================


@pytest.mark.asyncio
async def test_upload_extxyz(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test uploading an extxyz trajectory file stores frames."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms_list = [
        _make_atoms("H2", [[0, 0, 0], [1, 0, 0]]),
        _make_atoms("H3", [[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
    ]
    content = _atoms_to_file_bytes(atoms_list, "extxyz")

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("traj.extxyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 201

    result = FrameBulkResponse.model_validate(response.json())
    assert result.total == 2
    assert result.start == 0
    assert result.stop == 2

    # Verify frames stored
    assert await frame_storage.get_length(room.id) == 2


@pytest.mark.asyncio
async def test_upload_xyz(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test uploading an xyz format trajectory file."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms("H2", [[0, 0, 0], [1, 0, 0]])
    content = _atoms_to_file_bytes([atoms], "xyz")

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory?format=xyz",
        files={"file": ("traj.xyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 201

    result = FrameBulkResponse.model_validate(response.json())
    assert result.total == 1
    assert await frame_storage.get_length(room.id) == 1


@pytest.mark.asyncio
async def test_upload_format_from_extension(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test that format is inferred from file extension when not specified."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms("H2", [[0, 0, 0], [1, 0, 0]])
    content = _atoms_to_file_bytes([atoms], "xyz")

    # No format query param, but filename ends with .xyz
    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("trajectory.xyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    assert await frame_storage.get_length(room.id) == 1


@pytest.mark.asyncio
async def test_upload_explicit_format(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test that explicit format param overrides extension inference."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms("H2", [[0, 0, 0], [1, 0, 0]])
    content = _atoms_to_file_bytes([atoms], "extxyz")

    # Filename says .xyz but format param says extxyz
    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory?format=extxyz",
        files={"file": ("traj.xyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    assert await frame_storage.get_length(room.id) == 1


@pytest.mark.asyncio
async def test_upload_preserves_positions(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test that positions survive the upload roundtrip."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    original_positions = [[0.0, 0.0, 0.0], [1.5, 2.5, 3.5]]
    atoms = _make_atoms("H2", original_positions)
    content = _atoms_to_file_bytes([atoms], "extxyz")

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("traj.extxyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 201

    # Read back from storage and decode
    raw_frame = await frame_storage[room.id][0]
    assert raw_frame is not None
    decoded_atoms = decode(raw_frame)
    np.testing.assert_allclose(
        decoded_atoms.get_positions(), original_positions, atol=1e-6
    )


@pytest.mark.asyncio
async def test_upload_appends_to_nonempty(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Test that uploading to a room with existing frames appends."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Pre-populate with 2 frames
    existing = [_make_atoms(), _make_atoms()]
    await _add_atoms_to_storage(frame_storage, room.id, existing)
    assert await frame_storage.get_length(room.id) == 2

    # Upload 1 more frame
    new_atoms = _make_atoms("He", [[0, 0, 0]])
    content = _atoms_to_file_bytes([new_atoms], "extxyz")

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("traj.extxyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 201

    result = FrameBulkResponse.model_validate(response.json())
    assert result.total == 3
    assert result.start == 2
    assert result.stop == 3

    assert await frame_storage.get_length(room.id) == 3


@pytest.mark.asyncio
async def test_upload_empty_file(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test uploading an empty file returns 400."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("traj.extxyz", b"", "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 400

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == InvalidPayload.type_uri()


@pytest.mark.asyncio
async def test_upload_requires_auth(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test uploading without authentication returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    atoms = _make_atoms()
    content = _atoms_to_file_bytes([atoms], "extxyz")

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("traj.extxyz", content, "application/octet-stream")},
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_upload_locked_room(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Test uploading to a locked room returns 423."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    room.locked = True
    session.add(room)
    await session.commit()

    atoms = _make_atoms()
    content = _atoms_to_file_bytes([atoms], "extxyz")

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("traj.extxyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 423

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomLocked.type_uri()


@pytest.mark.asyncio
async def test_upload_provider_backed_readonly(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Upload to a provider-backed room returns 409 RoomReadOnly."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    await frame_storage.set_frame_count(room.id, 10)

    atoms = _make_atoms("H2", [[0, 0, 0], [1, 0, 0]])
    content = _atoms_to_file_bytes([atoms], "extxyz")

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("traj.extxyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 409

    problem = ProblemDetail.model_validate(response.json())
    assert problem.type == RoomReadOnly.type_uri()


# =============================================================================
# Download Token Tests
# =============================================================================


@pytest.mark.asyncio
async def test_create_download_token(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """POST creates a download token with default TTL."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _add_atoms_to_storage(frame_storage, room.id, [_make_atoms()])

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory/download-tokens",
        headers=auth_header(token),
    )
    assert response.status_code == 201

    data = response.json()
    assert "token" in data
    assert "url" in data
    assert "expires_in" in data
    assert data["expires_in"] == 300  # default TTL
    assert data["token"] in data["url"]
    assert f"/v1/rooms/{room.id}/trajectory" in data["url"]


@pytest.mark.asyncio
async def test_create_download_token_custom_ttl(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """POST with custom TTL sets that TTL."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _add_atoms_to_storage(frame_storage, room.id, [_make_atoms()])

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory/download-tokens",
        headers=auth_header(token),
        json={"ttl": 60},
    )
    assert response.status_code == 201
    assert response.json()["expires_in"] == 60


@pytest.mark.asyncio
async def test_create_download_token_ttl_exceeds_max_rejected(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """TTL above server max is rejected by Pydantic validation."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _add_atoms_to_storage(frame_storage, room.id, [_make_atoms()])

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory/download-tokens",
        headers=auth_header(token),
        json={"ttl": 999999},
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_create_download_token_requires_auth(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """POST without auth returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory/download-tokens",
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_download_with_token_no_auth_header(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """GET with valid download token works without Authorization header."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _add_atoms_to_storage(frame_storage, room.id, [_make_atoms()])

    # Create download token
    create_resp = await client.post(
        f"/v1/rooms/{room.id}/trajectory/download-tokens",
        headers=auth_header(token),
    )
    download_token = create_resp.json()["token"]

    # Download WITHOUT auth header, using token param
    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?token={download_token}",
    )
    assert response.status_code == 200

    parsed = _parse_trajectory(response.text)
    assert len(parsed) == 1


@pytest.mark.asyncio
async def test_download_with_invalid_token(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """GET with invalid token and no auth header returns 401."""
    user, _ = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _add_atoms_to_storage(frame_storage, room.id, [_make_atoms()])

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?token=bogus-token",
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_download_token_wrong_room(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Token for room A cannot download from room B."""
    user, token = await create_test_user_in_db(session)
    room_a = await create_test_room(session, user, description="Room A")
    room_b = await create_test_room(session, user, description="Room B")
    await _add_atoms_to_storage(frame_storage, room_a.id, [_make_atoms()])
    await _add_atoms_to_storage(frame_storage, room_b.id, [_make_atoms()])

    # Create token for room A
    create_resp = await client.post(
        f"/v1/rooms/{room_a.id}/trajectory/download-tokens",
        headers=auth_header(token),
    )
    download_token = create_resp.json()["token"]

    # Try to use it on room B
    response = await client.get(
        f"/v1/rooms/{room_b.id}/trajectory?token={download_token}",
    )
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_download_token_single_use(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Token is consumed after first use."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    await _add_atoms_to_storage(frame_storage, room.id, [_make_atoms()])

    create_resp = await client.post(
        f"/v1/rooms/{room.id}/trajectory/download-tokens",
        headers=auth_header(token),
    )
    download_token = create_resp.json()["token"]

    # First use succeeds
    resp1 = await client.get(
        f"/v1/rooms/{room.id}/trajectory?token={download_token}",
    )
    assert resp1.status_code == 200

    # Second use fails — token was consumed
    resp2 = await client.get(
        f"/v1/rooms/{room.id}/trajectory?token={download_token}",
    )
    assert resp2.status_code == 401


@pytest.mark.asyncio
async def test_upload_enriches_frames(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Uploaded bare atoms get colors, radii, and connectivity added."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Bare atoms — no colors, radii, or connectivity
    atoms = _make_atoms("H2", [[0, 0, 0], [1, 0, 0]])
    content = _atoms_to_file_bytes([atoms], "extxyz")

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory",
        files={"file": ("traj.extxyz", content, "application/octet-stream")},
        headers=auth_header(token),
    )
    assert response.status_code == 201

    raw_frame = await frame_storage[room.id][0]
    assert raw_frame is not None
    assert b"arrays.colors" in raw_frame
    assert b"arrays.radii" in raw_frame
    assert b"info.connectivity" in raw_frame


@pytest.mark.asyncio
async def test_upload_malformed_file(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Uploading a non-trajectory file returns 400."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/trajectory?format=extxyz",
        headers=auth_header(token),
        files={"file": ("garbage.xyz", b"this is not a trajectory", "text/plain")},
    )
    assert response.status_code == 400

    problem = ProblemDetail.model_validate(response.json())
    assert problem.detail is not None
    assert "parse" in problem.detail.lower() or "failed" in problem.detail.lower()


@pytest.mark.asyncio
async def test_download_atom_selection_out_of_range(
    client: AsyncClient,
    session: AsyncSession,
    frame_storage: FrameStorage,
) -> None:
    """Atom selection with out-of-range index returns 400."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    # H2 has 2 atoms (indices 0, 1)
    await _add_atoms_to_storage(frame_storage, room.id, [_make_atoms()])

    response = await client.get(
        f"/v1/rooms/{room.id}/trajectory?selection=0,99",
        headers=auth_header(token),
    )
    assert response.status_code == 400


# =============================================================================
# Section 4: Full integration tests (real server + real Redis)
# =============================================================================


def test_download_token_full_roundtrip(server: str) -> None:
    """Upload frames via ZnDraw client, get a download token, download, verify."""
    import uuid

    import httpx

    from zndraw import ZnDraw

    room_id = uuid.uuid4().hex
    client = ZnDraw(url=server, room=room_id)

    # Upload 3 frames
    uploaded = [
        ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]]),
        ase.Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
        ase.Atoms("O", positions=[[5, 5, 5]]),
    ]
    client.extend(uploaded)
    assert len(client) == 3

    # Request a download token using the client's JWT
    token_resp = httpx.post(
        f"{server}/v1/rooms/{room_id}/trajectory/download-tokens",
        headers={"Authorization": f"Bearer {client.api.token}"},
    )
    assert token_resp.status_code == 201
    token_data = token_resp.json()
    assert "token" in token_data
    assert "url" in token_data

    # Download using the token (no auth header) — like a browser would
    download_resp = httpx.get(
        f"{server}/v1/rooms/{room_id}/trajectory",
        params={"token": token_data["token"], "format": "extxyz"},
    )
    assert download_resp.status_code == 200

    # Parse the downloaded trajectory
    buf = io.StringIO(download_resp.text)
    downloaded = ase.io.read(buf, index=":", format="extxyz")
    if not isinstance(downloaded, list):
        downloaded = [downloaded]

    assert len(downloaded) == 3

    # Verify each frame
    for orig, dl in zip(uploaded, downloaded, strict=True):
        assert list(orig.symbols) == list(dl.symbols)
        np.testing.assert_allclose(orig.positions, dl.positions, atol=1e-6)

    client.disconnect()
