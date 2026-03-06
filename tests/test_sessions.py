"""Integration tests for vis.sessions (user-scoped session control)."""

import json
import uuid

import pytest
from fastapi_users.jwt import decode_jwt
from redis import Redis
from zndraw_auth import AuthSettings

from zndraw import ZnDraw
from zndraw.client import Sessions
from zndraw.geometries.camera import Camera
from zndraw.redis import RedisKey
from zndraw.settings import RoomConfig


def _get_user_id(vis: ZnDraw) -> str:
    """Extract the user UUID from the client's JWT token."""
    auth = AuthSettings()
    assert vis.api.token is not None
    payload = decode_jwt(
        vis.api.token,
        auth.secret_key.get_secret_value(),
        audience=["fastapi-users:auth"],
    )
    return payload["sub"]


def _seed_frontend_session(
    r: Redis, room_id: str, user_id: str, fake_sid: str, camera_key: str
) -> None:
    """Seed Redis with a fake frontend session for testing.

    Creates the presence key, active-cameras entry, and a camera geometry
    so that ``list_sessions`` returns the fake SID and camera operations work.
    """
    # Presence key (maps SID → user_id)
    r.set(RedisKey.presence_sid(room_id, fake_sid), user_id, ex=300)

    # Active cameras hash (marks SID as frontend)
    r.hset(RedisKey.active_cameras(room_id), fake_sid, camera_key)

    # Camera geometry in Redis hash
    camera = Camera(owner=user_id)
    camera_value = json.dumps(
        {"sid": fake_sid, "email": "test@local.test", "data": camera.model_dump()}
    )
    r.hset(RedisKey.room_cameras(room_id), camera_key, camera_value)


# =============================================================================
# Mapping Interface Tests
# =============================================================================


def test_sessions_returns_sessions_type(server: str):
    """vis.sessions returns a Sessions instance."""
    vis = ZnDraw(url=server, room=uuid.uuid4().hex)
    assert isinstance(vis.sessions, Sessions)
    vis.disconnect()


def test_sessions_empty_without_frontend(server: str):
    """Python client (pyclient) sees no sessions — only frontend counts."""
    vis = ZnDraw(url=server, room=uuid.uuid4().hex)
    assert len(vis.sessions) == 0
    assert list(vis.sessions) == []
    vis.disconnect()


def test_sessions_repr_empty(server: str):
    """repr shows empty list when no sessions."""
    vis = ZnDraw(url=server, room=uuid.uuid4().hex)
    assert repr(vis.sessions) == "Sessions([])"
    vis.disconnect()


def test_sessions_str_empty(server: str):
    """str shows count when no sessions."""
    vis = ZnDraw(url=server, room=uuid.uuid4().hex)
    assert str(vis.sessions) == "Sessions(n=0)"
    vis.disconnect()


def test_sessions_getitem_missing_raises_key_error(server: str):
    """Accessing non-existent session raises KeyError (validated in Session.__post_init__)."""
    vis = ZnDraw(url=server, room=uuid.uuid4().hex)
    with pytest.raises(KeyError):
        vis.sessions["nonexistent-sid"]
    vis.disconnect()


def test_sessions_contains_false_for_missing(server: str):
    """Membership test returns False for non-existent session."""
    vis = ZnDraw(url=server, room=uuid.uuid4().hex)
    assert "nonexistent-sid" not in vis.sessions
    vis.disconnect()


def test_session_repr_and_str(server: str):
    """Session has clean repr/str from dataclass (uses seeded session)."""
    room_id = uuid.uuid4().hex
    vis = ZnDraw(url=server, room=room_id)
    user_id = _get_user_id(vis)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-repr"
    camera_key = f"cam:test@local.test:{fake_sid[:8]}"
    _seed_frontend_session(r, room_id, user_id, fake_sid, camera_key)

    session = vis.sessions[fake_sid]
    assert repr(session) == f"Session(sid={fake_sid!r})"
    assert str(session) == f"Session({fake_sid!r})"

    r.close()
    vis.disconnect()


def test_sessions_is_frozen(server: str):
    """Sessions dataclass is frozen — cannot set attributes."""
    vis = ZnDraw(url=server, room=uuid.uuid4().hex)
    sessions = vis.sessions
    with pytest.raises(AttributeError):
        sessions._api = None  # type: ignore[misc]
    vis.disconnect()


# =============================================================================
# Roundtrip Tests (seeded frontend session)
# =============================================================================


def test_sessions_lists_seeded_frontend(server: str):
    """Seeded frontend session appears in vis.sessions."""
    room_id = uuid.uuid4().hex
    vis = ZnDraw(url=server, room=room_id)
    user_id = _get_user_id(vis)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-001"
    camera_key = f"cam:test@local.test:{fake_sid[:8]}"
    _seed_frontend_session(r, room_id, user_id, fake_sid, camera_key)

    assert fake_sid in vis.sessions
    assert len(vis.sessions) == 1
    assert list(vis.sessions) == [fake_sid]

    r.close()
    vis.disconnect()


def test_session_active_camera_roundtrip(server: str):
    """Get/set active_camera through Session proxy."""
    room_id = uuid.uuid4().hex
    vis = ZnDraw(url=server, room=room_id)
    user_id = _get_user_id(vis)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-002"
    camera_key = f"cam:test@local.test:{fake_sid[:8]}"
    _seed_frontend_session(r, room_id, user_id, fake_sid, camera_key)

    session = vis.sessions[fake_sid]
    assert session.active_camera == camera_key

    # Set a different active camera
    new_key = "custom_camera"
    session.active_camera = new_key
    assert session.active_camera == new_key

    r.close()
    vis.disconnect()


def test_session_settings_roundtrip(server: str):
    """Get/set settings through Session proxy."""
    room_id = uuid.uuid4().hex
    vis = ZnDraw(url=server, room=room_id)
    user_id = _get_user_id(vis)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-003"
    camera_key = f"cam:test@local.test:{fake_sid[:8]}"
    _seed_frontend_session(r, room_id, user_id, fake_sid, camera_key)

    session = vis.sessions[fake_sid]

    # Default settings
    settings = session.settings
    assert isinstance(settings, RoomConfig)
    assert settings.studio_lighting.key_light == pytest.approx(0.7)

    # Update settings
    new_settings = settings.model_copy(
        update={
            "studio_lighting": settings.studio_lighting.model_copy(
                update={"key_light": 1.5}
            )
        }
    )
    session.settings = new_settings

    # Read back
    updated = session.settings
    assert updated.studio_lighting.key_light == pytest.approx(1.5)

    r.close()
    vis.disconnect()


def test_session_camera_roundtrip(server: str):
    """Get/set camera through Session proxy (resolves via active_camera)."""
    room_id = uuid.uuid4().hex
    vis = ZnDraw(url=server, room=room_id)
    user_id = _get_user_id(vis)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-004"
    camera_key = f"cam:test@local.test:{fake_sid[:8]}"
    _seed_frontend_session(r, room_id, user_id, fake_sid, camera_key)

    session = vis.sessions[fake_sid]

    # Read default camera
    cam = session.camera
    assert isinstance(cam, Camera)
    assert cam.position == pytest.approx((-10.0, 10.0, 30.0))

    # Update camera
    new_cam = cam.model_copy(update={"position": (5.0, 5.0, 5.0)})
    session.camera = new_cam

    # Read back
    updated_cam = session.camera
    assert updated_cam.position == pytest.approx((5.0, 5.0, 5.0))

    r.close()
    vis.disconnect()


def test_session_settings_not_accessible_by_other_user(server: str):
    """Another user cannot access a session's settings (404)."""
    room_id = uuid.uuid4().hex
    vis1 = ZnDraw(url=server, room=room_id)
    user1_id = _get_user_id(vis1)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-005"
    camera_key = f"cam:test@local.test:{fake_sid[:8]}"
    _seed_frontend_session(r, room_id, user1_id, fake_sid, camera_key)

    # Create a second user client in the same room
    vis2 = ZnDraw(url=server, room=room_id)

    # vis2 should not see vis1's sessions
    assert fake_sid not in vis2.sessions
    assert len(vis2.sessions) == 0

    r.close()
    vis2.disconnect()
    vis1.disconnect()
