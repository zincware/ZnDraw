"""Tests for Screenshot REST API endpoints."""

import json
from pathlib import Path

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.redis import RedisKey
from zndraw.schemas import StatusResponse

# =============================================================================
# Fixture: adds media_path override on top of the shared client
# =============================================================================


@pytest.fixture(name="media_path")
def media_path_fixture(tmp_path: Path):
    """Provide a temporary media directory and install the dependency override.

    The shared ``client`` fixture (from conftest) sets up the base overrides.
    This fixture layers ``get_media_path`` on top — tests that need it request
    both ``client`` and ``media_path``.
    """
    from zndraw.app import app
    from zndraw.dependencies import get_media_path

    media = tmp_path / "media"
    app.dependency_overrides[get_media_path] = lambda: media
    yield media
    app.dependency_overrides.pop(get_media_path, None)


# =============================================================================
# Helpers unique to this test file
# =============================================================================


def _png_bytes(size: int = 100) -> bytes:
    """Create minimal valid-ish PNG bytes for testing."""
    return b"\x89PNG\r\n\x1a\n" + b"\x00" * size


def _make_camera_entry(sid: str, owner_id: str, email: str = "u@test") -> str:
    """Build a JSON camera entry as stored in room_cameras hash."""
    from zndraw.geometries.camera import Camera

    camera = Camera(owner=owner_id)
    return json.dumps({"sid": sid, "email": email, "data": camera.model_dump()})


# =============================================================================
# POST /upload — Upload screenshot
# =============================================================================


@pytest.mark.asyncio
async def test_upload_screenshot(
    client: AsyncClient,
    session: AsyncSession,
    media_path: Path,
) -> None:
    """Upload a PNG screenshot, verify 201 and file on disk."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    data = _png_bytes()
    response = await client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", data, "image/png")},
        data={"format": "png"},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    body = response.json()
    assert body["format"] == "png"
    assert body["size"] == len(data)
    assert body["status"] == "completed"
    assert body["room_id"] == room.id

    # Verify file exists on disk
    file_path = media_path / room.id / "screenshots" / f"{body['id']}.png"
    assert file_path.exists()
    assert file_path.read_bytes() == data


@pytest.mark.asyncio
async def test_upload_invalid_format(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Upload with unsupported format returns 422."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.bmp", b"data", "image/bmp")},
        data={"format": "bmp"},
        headers=auth_header(token),
    )
    assert response.status_code == 422
    assert "invalid-screenshot-format" in response.json()["type"]


@pytest.mark.asyncio
async def test_upload_too_large(client: AsyncClient, session: AsyncSession) -> None:
    """Upload exceeding 10 MB returns 413."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    big_data = b"\x00" * (10 * 1024 * 1024 + 1)
    response = await client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", big_data, "image/png")},
        data={"format": "png"},
        headers=auth_header(token),
    )
    assert response.status_code == 413
    assert "screenshot-too-large" in response.json()["type"]


# =============================================================================
# GET — List screenshots
# =============================================================================


@pytest.mark.asyncio
async def test_list_screenshots(client: AsyncClient, session: AsyncSession) -> None:
    """Upload 3 screenshots, verify list with pagination."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    for i in range(3):
        await client.post(
            f"/v1/rooms/{room.id}/screenshots/upload",
            files={"file": (f"shot{i}.png", _png_bytes(50 + i), "image/png")},
            data={"format": "png"},
            headers=auth_header(token),
        )

    response = await client.get(
        f"/v1/rooms/{room.id}/screenshots?limit=2&offset=0",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    body = response.json()
    assert body["total"] == 3
    assert len(body["items"]) == 2
    assert body["limit"] == 2
    assert body["offset"] == 0

    # Second page
    response2 = await client.get(
        f"/v1/rooms/{room.id}/screenshots?limit=2&offset=2",
        headers=auth_header(token),
    )
    assert response2.status_code == 200
    body2 = response2.json()
    assert len(body2["items"]) == 1


# =============================================================================
# GET — Single screenshot
# =============================================================================


@pytest.mark.asyncio
async def test_get_screenshot(client: AsyncClient, session: AsyncSession) -> None:
    """Upload then GET, verify data field is base64-encoded."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)
    data = _png_bytes()

    upload = await client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", data, "image/png")},
        data={"format": "png"},
        headers=auth_header(token),
    )
    screenshot_id = upload.json()["id"]

    response = await client.get(
        f"/v1/rooms/{room.id}/screenshots/{screenshot_id}",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    body = response.json()
    assert body["status"] == "completed"
    assert body["data"] is not None

    import base64

    assert base64.b64decode(body["data"]) == data


@pytest.mark.asyncio
async def test_get_screenshot_not_found(
    client: AsyncClient, session: AsyncSession
) -> None:
    """GET nonexistent screenshot returns 404."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.get(
        f"/v1/rooms/{room.id}/screenshots/99999",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "screenshot-not-found" in response.json()["type"]


# =============================================================================
# DELETE — Delete screenshot
# =============================================================================


@pytest.mark.asyncio
async def test_delete_screenshot(
    client: AsyncClient,
    session: AsyncSession,
    media_path: Path,
) -> None:
    """Upload, DELETE, verify 204 and file removed."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    upload = await client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", _png_bytes(), "image/png")},
        data={"format": "png"},
        headers=auth_header(token),
    )
    screenshot_id = upload.json()["id"]
    file_path = media_path / room.id / "screenshots" / f"{screenshot_id}.png"
    assert file_path.exists()

    response = await client.delete(
        f"/v1/rooms/{room.id}/screenshots/{screenshot_id}",
        headers=auth_header(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())
    assert not file_path.exists()


# =============================================================================
# POST — Request capture (JSON)
# =============================================================================


@pytest.mark.asyncio
async def test_request_capture(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
    redis_client,
) -> None:
    """JSON POST with live session owned by requesting user returns 202."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    target_sid = "frontend-sid-123"
    cam_key = f"cam:{user.email}:abcd1234"

    # Pre-populate Redis: mark session as active and add camera entry
    await redis_client.hset(RedisKey.active_cameras(room.id), target_sid, cam_key)
    await redis_client.hset(
        RedisKey.room_cameras(room.id),
        cam_key,
        _make_camera_entry(target_sid, str(user.id), user.email),
    )

    response = await client.post(
        f"/v1/rooms/{room.id}/screenshots",
        json={"session_id": target_sid},
        headers=auth_header(token),
    )
    assert response.status_code == 202
    body = response.json()
    assert body["status"] == "pending"
    assert body["data"] is None
    assert "Location" in response.headers

    # Verify socket emit
    assert len(mock_sio.emitted) == 1
    emit = mock_sio.emitted[0]
    assert emit["event"] == "screenshot_request"
    assert emit["data"]["screenshot_id"] == body["id"]
    assert emit["to"] == target_sid


@pytest.mark.asyncio
async def test_request_capture_rejects_other_users_session(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
    redis_client,
) -> None:
    """Requesting a screenshot from another user's session returns 409."""
    user_a, token_a = await create_test_user_in_db(session, email="alice@test")
    room = await create_test_room(session, user_a)
    user_b, _ = await create_test_user_in_db(session, email="bob@test")

    target_sid = "bobs-browser-sid"
    cam_key = "cam:bob@test:abcd1234"

    # Session exists but is owned by user_b, not user_a
    await redis_client.hset(RedisKey.active_cameras(room.id), target_sid, cam_key)
    await redis_client.hset(
        RedisKey.room_cameras(room.id),
        cam_key,
        _make_camera_entry(target_sid, str(user_b.id), user_b.email),
    )

    response = await client.post(
        f"/v1/rooms/{room.id}/screenshots",
        json={"session_id": target_sid},
        headers=auth_header(token_a),
    )
    assert response.status_code == 409
    assert "no-frontend-session" in response.json()["type"]

    # No socket event should have been emitted
    assert len(mock_sio.emitted) == 0


@pytest.mark.asyncio
async def test_request_capture_invalid_session(
    client: AsyncClient,
    session: AsyncSession,
) -> None:
    """Request capture with non-active session returns 409."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Session not active (Redis has no entry for this sid)
    response = await client.post(
        f"/v1/rooms/{room.id}/screenshots",
        json={"session_id": "nonexistent-sid"},
        headers=auth_header(token),
    )
    assert response.status_code == 409
    assert "no-frontend-session" in response.json()["type"]


# =============================================================================
# PATCH — Complete pending screenshot
# =============================================================================


@pytest.mark.asyncio
async def test_patch_pending_screenshot(
    client: AsyncClient,
    session: AsyncSession,
    redis_client,
    media_path: Path,
) -> None:
    """Create pending via capture request, then PATCH with file."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    target_sid = "frontend-sid"
    cam_key = f"cam:{user.email}:abcd1234"

    await redis_client.hset(RedisKey.active_cameras(room.id), target_sid, cam_key)
    await redis_client.hset(
        RedisKey.room_cameras(room.id),
        cam_key,
        _make_camera_entry(target_sid, str(user.id), user.email),
    )

    # Create pending screenshot
    capture_resp = await client.post(
        f"/v1/rooms/{room.id}/screenshots",
        json={"session_id": target_sid},
        headers=auth_header(token),
    )
    screenshot_id = capture_resp.json()["id"]

    # PATCH to complete
    data = _png_bytes(200)
    response = await client.patch(
        f"/v1/rooms/{room.id}/screenshots/{screenshot_id}",
        files={"file": ("shot.png", data, "image/png")},
        data={"format": "png", "width": "1920", "height": "1080"},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    body = response.json()
    assert body["status"] == "completed"
    assert body["size"] == len(data)
    assert body["width"] == 1920
    assert body["height"] == 1080

    # Verify file on disk
    file_path = media_path / room.id / "screenshots" / f"{screenshot_id}.png"
    assert file_path.exists()


@pytest.mark.asyncio
async def test_patch_completed_screenshot(
    client: AsyncClient, session: AsyncSession
) -> None:
    """PATCH on already-completed screenshot returns 409."""
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    # Upload a completed screenshot
    upload = await client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", _png_bytes(), "image/png")},
        data={"format": "png"},
        headers=auth_header(token),
    )
    screenshot_id = upload.json()["id"]

    # Try to PATCH (should fail — already completed)
    response = await client.patch(
        f"/v1/rooms/{room.id}/screenshots/{screenshot_id}",
        files={"file": ("shot2.png", _png_bytes(), "image/png")},
        data={"format": "png"},
        headers=auth_header(token),
    )
    assert response.status_code == 409
    assert "screenshot-not-pending" in response.json()["type"]
