"""E2E test: Guest can authenticate, create a room, and write/read frames via REST.

Uses a real uvicorn server via the http_client fixture (which starts the server).
All interaction is through the REST API only (no Socket.IO).
"""

import uuid

import ase
import msgpack
import pytest
from httpx import AsyncClient

from zndraw.client import atoms_to_json_dict


def _make_atoms(x: float, formula: str = "H") -> ase.Atoms:
    """Create an Atoms object for testing."""
    atoms = ase.Atoms(formula, positions=[[x, 0, 0]])
    atoms.info["x_pos"] = x
    return atoms


def _decode_msgpack_frames(content: bytes) -> list:
    """Decode a msgpack-encoded list of frames from a response body."""
    return msgpack.unpackb(content, raw=True)


# =============================================================================
# Guest Auth + Room Flow
# =============================================================================


@pytest.mark.asyncio
async def test_guest_auth_returns_token(http_client: AsyncClient):
    """POST /v1/auth/guest returns a valid bearer token."""
    response = await http_client.post("/v1/auth/guest")
    assert response.status_code == 200

    body = response.json()
    assert body["token_type"] == "bearer"
    assert body["access_token"]
    assert isinstance(body["access_token"], str)
    assert len(body["access_token"]) > 10


@pytest.mark.asyncio
async def test_guest_can_create_room(http_client: AsyncClient):
    """A guest can create a room after authenticating."""
    # Get guest token
    auth_resp = await http_client.post("/v1/auth/guest")
    assert auth_resp.status_code == 200
    token = auth_resp.json()["access_token"]

    room_id = uuid.uuid4().hex
    create_resp = await http_client.post(
        "/v1/rooms",
        json={"room_id": room_id},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert create_resp.status_code == 201

    body = create_resp.json()
    assert body["room_id"] == room_id
    assert body["created"] is True


@pytest.mark.asyncio
async def test_guest_write_then_read_frame(http_client: AsyncClient):
    """Guest creates a room, writes a frame via REST, then reads it back."""
    # Step 1: Get guest token
    auth_resp = await http_client.post("/v1/auth/guest")
    assert auth_resp.status_code == 200
    token = auth_resp.json()["access_token"]
    headers = {"Authorization": f"Bearer {token}"}

    # Step 2: Create room (with @none so it starts empty)
    room_id = uuid.uuid4().hex
    create_resp = await http_client.post(
        "/v1/rooms",
        json={"room_id": room_id, "copy_from": "@none"},
        headers=headers,
    )
    assert create_resp.status_code == 201

    # Step 3: Write a frame via REST
    frame_data = atoms_to_json_dict(_make_atoms(7.5))
    post_resp = await http_client.post(
        f"/v1/rooms/{room_id}/frames",
        json={"frames": [frame_data]},
        headers=headers,
    )
    assert post_resp.status_code == 201

    result = post_resp.json()
    assert result["total"] == 1
    assert result["start"] == 0
    assert result["stop"] == 1

    # Step 4: Read the frame back via REST
    get_resp = await http_client.get(
        f"/v1/rooms/{room_id}/frames",
        headers=headers,
    )
    assert get_resp.status_code == 200
    assert get_resp.headers["content-type"] == "application/x-msgpack"

    frames = _decode_msgpack_frames(get_resp.content)
    assert len(frames) == 1


@pytest.mark.asyncio
async def test_guest_write_multiple_frames(http_client: AsyncClient):
    """Guest writes multiple frames and reads them back by index."""
    auth_resp = await http_client.post("/v1/auth/guest")
    token = auth_resp.json()["access_token"]
    headers = {"Authorization": f"Bearer {token}"}

    room_id = uuid.uuid4().hex
    await http_client.post(
        "/v1/rooms",
        json={"room_id": room_id, "copy_from": "@none"},
        headers=headers,
    )

    # Write 3 frames
    frames = [atoms_to_json_dict(_make_atoms(float(i))) for i in range(3)]
    post_resp = await http_client.post(
        f"/v1/rooms/{room_id}/frames",
        json={"frames": frames},
        headers=headers,
    )
    assert post_resp.status_code == 201
    assert post_resp.json()["total"] == 3

    # Read all frames back
    get_resp = await http_client.get(
        f"/v1/rooms/{room_id}/frames",
        headers=headers,
    )
    assert get_resp.status_code == 200
    stored_frames = _decode_msgpack_frames(get_resp.content)
    assert len(stored_frames) == 3


@pytest.mark.asyncio
async def test_guest_cannot_access_other_room_without_auth(http_client: AsyncClient):
    """Unauthenticated requests to frame endpoints return 401."""
    # Create a room with auth
    auth_resp = await http_client.post("/v1/auth/guest")
    token = auth_resp.json()["access_token"]
    headers = {"Authorization": f"Bearer {token}"}

    room_id = uuid.uuid4().hex
    await http_client.post(
        "/v1/rooms",
        json={"room_id": room_id, "copy_from": "@none"},
        headers=headers,
    )

    # Try to access frames without authorization
    get_resp = await http_client.get(f"/v1/rooms/{room_id}/frames")
    assert get_resp.status_code == 401
