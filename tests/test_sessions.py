"""Integration tests for vis.sessions (room-scoped session visibility)."""

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
    r: Redis,
    room_id: str,
    user_id: str,
    fake_sid: str,
    camera_key: str,
    email: str = "test@local.test",
) -> None:
    """Seed Redis with a fake frontend session for testing.

    Creates active-cameras entry and a camera geometry so that
    ``list_sessions`` returns the fake SID and camera operations work.
    """
    # Active cameras hash (marks SID as frontend)
    r.hset(RedisKey.active_cameras(room_id), fake_sid, camera_key)

    # Camera geometry in Redis hash
    camera = Camera(owner=user_id)
    camera_value = json.dumps(
        {"sid": fake_sid, "email": email, "data": camera.model_dump()}
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

    # Seed a second camera and switch to it
    second_key = "cam:other@local.test:second"
    second_cam = Camera(owner=user_id)
    r.hset(
        RedisKey.room_cameras(room_id),
        second_key,
        json.dumps(
            {
                "sid": "other-sid",
                "email": "other@local.test",
                "data": second_cam.model_dump(),
            }
        ),
    )
    session.active_camera = second_key
    assert session.active_camera == second_key

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


# =============================================================================
# Cross-User Visibility Tests (room-scoped)
# =============================================================================


def test_cross_user_sees_other_users_sessions(server: str):
    """User2 can see user1's frontend session in the listing."""
    room_id = uuid.uuid4().hex
    vis1 = ZnDraw(url=server, room=room_id)
    user1_id = _get_user_id(vis1)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-cross-list"
    camera_key = f"cam:user1@local.test:{fake_sid[:8]}"
    _seed_frontend_session(
        r, room_id, user1_id, fake_sid, camera_key, email="user1@local.test"
    )

    vis2 = ZnDraw(url=server, room=room_id)

    # User2 can see user1's session
    assert fake_sid in vis2.sessions
    assert len(vis2.sessions) == 1

    # The listing includes email and camera_key
    items = vis2.api.list_sessions()
    assert len(items) == 1
    assert items[0].sid == fake_sid
    assert items[0].email == "user1@local.test"
    assert items[0].camera_key == camera_key

    r.close()
    vis2.disconnect()
    vis1.disconnect()


def test_cross_user_can_read_active_camera(server: str):
    """User2 can read user1's session active camera (shared state)."""
    room_id = uuid.uuid4().hex
    vis1 = ZnDraw(url=server, room=room_id)
    user1_id = _get_user_id(vis1)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-cross-read"
    camera_key = f"cam:user1@local.test:{fake_sid[:8]}"
    _seed_frontend_session(r, room_id, user1_id, fake_sid, camera_key)

    vis2 = ZnDraw(url=server, room=room_id)
    active = vis2.api.get_active_camera(fake_sid)
    assert active == camera_key

    r.close()
    vis2.disconnect()
    vis1.disconnect()


def test_cross_user_cannot_set_active_camera(server: str):
    """User2 cannot set active camera on user1's session (KeyError from SessionNotFound)."""
    room_id = uuid.uuid4().hex
    vis1 = ZnDraw(url=server, room=room_id)
    user1_id = _get_user_id(vis1)

    r = Redis.from_url("redis://localhost", decode_responses=True)
    fake_sid = "fake-frontend-cross-write"
    camera_key = f"cam:user1@local.test:{fake_sid[:8]}"
    _seed_frontend_session(r, room_id, user1_id, fake_sid, camera_key)

    vis2 = ZnDraw(url=server, room=room_id)

    # Can see the session and read it
    session = vis2.sessions[fake_sid]
    assert session.active_camera == camera_key

    # Cannot set active camera (VerifiedSessionDep rejects — raises KeyError)
    with pytest.raises(KeyError):
        session.active_camera = camera_key

    r.close()
    vis2.disconnect()
    vis1.disconnect()
