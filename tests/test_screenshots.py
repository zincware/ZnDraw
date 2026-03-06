"""Tests for Screenshot REST API endpoints."""

from collections.abc import AsyncIterator
from pathlib import Path
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
from zndraw.models import MemberRole, Room, RoomMembership
from zndraw.schemas import StatusResponse

# =============================================================================
# Fixtures
# =============================================================================


@pytest_asyncio.fixture(name="ss_session")
async def ss_session_fixture() -> AsyncIterator[AsyncSession]:
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


@pytest.fixture(name="media_path")
def media_path_fixture(tmp_path: Path) -> Path:
    """Provide a temporary media directory for screenshots."""
    return tmp_path / "media"


@pytest.fixture(name="mock_redis")
def mock_redis_fixture() -> AsyncMock:
    """Provide a mock Redis with configurable hexists."""
    mock = AsyncMock()
    mock.get = AsyncMock(return_value=None)
    mock.hexists = AsyncMock(return_value=False)
    return mock


@pytest_asyncio.fixture(name="ss_client")
async def ss_client_fixture(
    ss_session: AsyncSession,
    mock_sio: MockSioServer,
    media_path: Path,
    mock_redis: AsyncMock,
) -> AsyncIterator[AsyncClient]:
    from zndraw_auth import get_session

    from zndraw.app import app
    from zndraw.dependencies import get_media_path, get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield ss_session

    def get_sio_override() -> MockSioServer:
        return mock_sio

    settings = Settings()
    settings.media_path = media_path  # type: ignore[assignment]
    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_tsio] = get_sio_override
    app.dependency_overrides[get_redis] = lambda: mock_redis
    app.dependency_overrides[get_media_path] = lambda: media_path
    app.state.settings = settings
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


def _auth(token: str) -> dict[str, str]:
    return {"Authorization": f"Bearer {token}"}


def _png_bytes(size: int = 100) -> bytes:
    """Create minimal valid-ish PNG bytes for testing."""
    return b"\x89PNG\r\n\x1a\n" + b"\x00" * size


# =============================================================================
# POST /upload — Upload screenshot
# =============================================================================


@pytest.mark.asyncio
async def test_upload_screenshot(
    ss_client: AsyncClient,
    ss_session: AsyncSession,
    media_path: Path,
) -> None:
    """Upload a PNG screenshot, verify 201 and file on disk."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    data = _png_bytes()
    response = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", data, "image/png")},
        data={"format": "png"},
        headers=_auth(token),
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
    ss_client: AsyncClient, ss_session: AsyncSession
) -> None:
    """Upload with unsupported format returns 422."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    response = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.bmp", b"data", "image/bmp")},
        data={"format": "bmp"},
        headers=_auth(token),
    )
    assert response.status_code == 422
    assert "invalid-screenshot-format" in response.json()["type"]


@pytest.mark.asyncio
async def test_upload_too_large(
    ss_client: AsyncClient, ss_session: AsyncSession
) -> None:
    """Upload exceeding 10 MB returns 413."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    big_data = b"\x00" * (10 * 1024 * 1024 + 1)
    response = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", big_data, "image/png")},
        data={"format": "png"},
        headers=_auth(token),
    )
    assert response.status_code == 413
    assert "screenshot-too-large" in response.json()["type"]


# =============================================================================
# GET — List screenshots
# =============================================================================


@pytest.mark.asyncio
async def test_list_screenshots(
    ss_client: AsyncClient, ss_session: AsyncSession
) -> None:
    """Upload 3 screenshots, verify list with pagination."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    for i in range(3):
        await ss_client.post(
            f"/v1/rooms/{room.id}/screenshots/upload",
            files={"file": (f"shot{i}.png", _png_bytes(50 + i), "image/png")},
            data={"format": "png"},
            headers=_auth(token),
        )

    response = await ss_client.get(
        f"/v1/rooms/{room.id}/screenshots?limit=2&offset=0",
        headers=_auth(token),
    )
    assert response.status_code == 200
    body = response.json()
    assert body["total"] == 3
    assert len(body["items"]) == 2
    assert body["limit"] == 2
    assert body["offset"] == 0

    # Second page
    response2 = await ss_client.get(
        f"/v1/rooms/{room.id}/screenshots?limit=2&offset=2",
        headers=_auth(token),
    )
    assert response2.status_code == 200
    body2 = response2.json()
    assert len(body2["items"]) == 1


# =============================================================================
# GET — Single screenshot
# =============================================================================


@pytest.mark.asyncio
async def test_get_screenshot(ss_client: AsyncClient, ss_session: AsyncSession) -> None:
    """Upload then GET, verify data field is base64-encoded."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)
    data = _png_bytes()

    upload = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", data, "image/png")},
        data={"format": "png"},
        headers=_auth(token),
    )
    screenshot_id = upload.json()["id"]

    response = await ss_client.get(
        f"/v1/rooms/{room.id}/screenshots/{screenshot_id}",
        headers=_auth(token),
    )
    assert response.status_code == 200
    body = response.json()
    assert body["status"] == "completed"
    assert body["data"] is not None

    import base64

    assert base64.b64decode(body["data"]) == data


@pytest.mark.asyncio
async def test_get_screenshot_not_found(
    ss_client: AsyncClient, ss_session: AsyncSession
) -> None:
    """GET nonexistent screenshot returns 404."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    response = await ss_client.get(
        f"/v1/rooms/{room.id}/screenshots/99999",
        headers=_auth(token),
    )
    assert response.status_code == 404
    assert "screenshot-not-found" in response.json()["type"]


# =============================================================================
# DELETE — Delete screenshot
# =============================================================================


@pytest.mark.asyncio
async def test_delete_screenshot(
    ss_client: AsyncClient,
    ss_session: AsyncSession,
    media_path: Path,
) -> None:
    """Upload, DELETE, verify 204 and file removed."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    upload = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", _png_bytes(), "image/png")},
        data={"format": "png"},
        headers=_auth(token),
    )
    screenshot_id = upload.json()["id"]
    file_path = media_path / room.id / "screenshots" / f"{screenshot_id}.png"
    assert file_path.exists()

    response = await ss_client.delete(
        f"/v1/rooms/{room.id}/screenshots/{screenshot_id}",
        headers=_auth(token),
    )
    assert response.status_code == 200
    StatusResponse.model_validate(response.json())
    assert not file_path.exists()


# =============================================================================
# POST — Request capture (JSON)
# =============================================================================


@pytest.mark.asyncio
async def test_request_capture(
    ss_client: AsyncClient,
    ss_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> None:
    """JSON POST with live session returns 202 and emits socket event."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    # Mark session as active frontend
    mock_redis.hexists = AsyncMock(return_value=True)

    response = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots",
        json={"session_id": "frontend-sid-123"},
        headers=_auth(token),
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
    assert emit["to"] == "frontend-sid-123"


@pytest.mark.asyncio
async def test_request_capture_invalid_session(
    ss_client: AsyncClient,
    ss_session: AsyncSession,
    mock_redis: AsyncMock,
) -> None:
    """Request capture with non-active session returns 409."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    # Session not active (default mock_redis.hexists returns False)
    response = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots",
        json={"session_id": "nonexistent-sid"},
        headers=_auth(token),
    )
    assert response.status_code == 409
    assert "no-frontend-session" in response.json()["type"]


# =============================================================================
# PATCH — Complete pending screenshot
# =============================================================================


@pytest.mark.asyncio
async def test_patch_pending_screenshot(
    ss_client: AsyncClient,
    ss_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
    media_path: Path,
) -> None:
    """Create pending via capture request, then PATCH with file."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)
    mock_redis.hexists = AsyncMock(return_value=True)

    # Create pending screenshot
    capture_resp = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots",
        json={"session_id": "frontend-sid"},
        headers=_auth(token),
    )
    screenshot_id = capture_resp.json()["id"]

    # PATCH to complete
    data = _png_bytes(200)
    response = await ss_client.patch(
        f"/v1/rooms/{room.id}/screenshots/{screenshot_id}",
        files={"file": ("shot.png", data, "image/png")},
        data={"format": "png", "width": "1920", "height": "1080"},
        headers=_auth(token),
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
    ss_client: AsyncClient, ss_session: AsyncSession, media_path: Path
) -> None:
    """PATCH on already-completed screenshot returns 409."""
    user, token = await _create_user(ss_session)
    room = await _create_room(ss_session, user)

    # Upload a completed screenshot
    upload = await ss_client.post(
        f"/v1/rooms/{room.id}/screenshots/upload",
        files={"file": ("shot.png", _png_bytes(), "image/png")},
        data={"format": "png"},
        headers=_auth(token),
    )
    screenshot_id = upload.json()["id"]

    # Try to PATCH (should fail — already completed)
    response = await ss_client.patch(
        f"/v1/rooms/{room.id}/screenshots/{screenshot_id}",
        files={"file": ("shot2.png", _png_bytes(), "image/png")},
        data={"format": "png"},
        headers=_auth(token),
    )
    assert response.status_code == 409
    assert "screenshot-not-pending" in response.json()["type"]
