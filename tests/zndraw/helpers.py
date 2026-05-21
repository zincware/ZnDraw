"""Shared test utility functions and classes.

These are importable helpers (not pytest fixtures).
Fixtures live in conftest.py.
"""

import asyncio
from typing import Any

import msgpack
from fastapi_users.password import PasswordHelper
from httpx import AsyncClient
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.models import Room
from zndraw.storage import RawFrame
from zndraw_auth import User
from zndraw_auth.settings import AuthSettings

_password_helper = PasswordHelper()


def make_raw_frame(data: dict) -> RawFrame:
    """Convert a simple dict to RawFrame format for test assertions.

    The storage converts input dicts to raw bytes format.
    This helper creates the expected output for simple test dicts.
    """
    result: RawFrame = {}
    for k, v in data.items():
        key: bytes = k.encode() if isinstance(k, str) else k
        packed = msgpack.packb(v)
        val: bytes = packed if packed is not None else b""
        result[key] = val
    return result


def create_test_user_model(
    email: str = "testuser@local.test",
    password: str = "testpassword",  # noqa: S107
    is_superuser: bool = False,
) -> User:
    """Create a User model instance with hashed password for tests."""
    return User(
        email=email,
        hashed_password=_password_helper.hash(password),
        is_active=True,
        is_superuser=is_superuser,
        is_verified=True,
    )


def create_test_token(user: User) -> str:
    """Create a JWT token for a test user."""
    from fastapi_users.jwt import generate_jwt

    settings = AuthSettings()
    data = {"sub": str(user.id), "aud": "fastapi-users:auth"}
    return generate_jwt(
        data,
        settings.secret_key.get_secret_value(),
        settings.token_lifetime_seconds,
    )


def decode_msgpack_response(content: bytes) -> list[dict[bytes, bytes]]:
    """Decode a MessagePack response body to a list of raw frames."""
    return msgpack.unpackb(content, raw=True)


async def create_test_user_in_db(
    session: AsyncSession,
    email: str = "testuser@local.test",
    *,
    is_superuser: bool = False,
) -> tuple[User, str]:
    """Create a user in the DB and return (user, token)."""
    user = create_test_user_model(email=email, is_superuser=is_superuser)
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user, create_test_token(user)


async def create_test_room(
    session: AsyncSession,
    user: User,
    description: str = "Test Room",
    room_name: str = "test-room",
) -> Room:
    """Create a PUBLIC user-owned room and return it."""
    from zndraw.access import Visibility

    room = Room(
        room_name=room_name,
        description=description,
        created_by_id=user.id,  # type: ignore[arg-type]
        owner_user_id=user.id,  # type: ignore[arg-type]
        visibility=Visibility.PUBLIC,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)
    return room


def auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


async def _register_and_login(
    client: AsyncClient, email: str, password: str = "test12345"
) -> str:
    """Register a new user and return their JWT access token."""
    await client.post("/v1/auth/register", json={"email": email, "password": password})
    r = await client.post(
        "/v1/auth/jwt/login",
        data={"username": email, "password": password},
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    r.raise_for_status()
    return r.json()["access_token"]


async def get_user_id(client: AsyncClient, token: str) -> str:
    """Return the UUID string of the authenticated user."""
    r = await client.get("/v1/auth/users/me", headers={"Authorization": f"Bearer {token}"})
    r.raise_for_status()
    return r.json()["id"]


async def create_room_via_api(
    client: AsyncClient,
    token: str,
    name: str,
    owner_id: str | None = None,
    visibility: str = "public",
    copy_from: str | None = None,
    description: str | None = None,
) -> str:
    """Create a room via POST /v1/rooms and return its composed room_id.

    If owner_id is omitted, the authenticated user's id is used.
    """
    if owner_id is None:
        owner_id = await get_user_id(client, token)
    body: dict = {"owner_id": owner_id, "name": name, "visibility": visibility}
    if copy_from is not None:
        body["copy_from"] = copy_from
    if description is not None:
        body["description"] = description
    r = await client.post(
        "/v1/rooms",
        json=body,
        headers={"Authorization": f"Bearer {token}"},
    )
    r.raise_for_status()
    return r.json()["room_id"]


class MockSioServer:
    """Mock Socket.IO server for testing broadcasts.

    Compatible with zndraw-socketio's AsyncServerWrapper emit pattern:
    - sio.emit(PydanticModel(), room=...) - model as first arg
    - sio.emit("event", data, room=...) - classic pattern
    """

    def __init__(self) -> None:
        self.emitted: list[dict[str, Any]] = []
        self.rooms: dict[str, set[str]] = {}
        self.sessions: dict[str, dict[str, Any]] = {}

    async def emit(
        self,
        event_or_model: str | BaseModel,
        data: Any = None,
        *,
        room: str | None = None,
        _skip_sid: str | None = None,
        to: str | None = None,
        **_kwargs: Any,
    ) -> None:
        if isinstance(event_or_model, BaseModel):
            # zndraw-socketio pattern: model class name -> snake_case event
            cls_name = type(event_or_model).__name__
            event = "".join(
                f"_{c.lower()}" if c.isupper() else c for c in cls_name
            ).lstrip("_")
            data = event_or_model.model_dump()
        else:
            event = event_or_model
        self.emitted.append({"event": event, "data": data, "room": room, "to": to})

    async def enter_room(self, sid: str, room: str) -> None:
        if room not in self.rooms:
            self.rooms[room] = set()
        self.rooms[room].add(sid)

    async def get_session(self, sid: str) -> dict[str, Any]:
        return self.sessions.get(sid, {})

    async def save_session(self, sid: str, session: dict[str, Any]) -> None:
        self.sessions[sid] = session


class InMemoryResultBackend:
    """In-memory result backend for testing.

    Drop-in replacement for the real Redis-based ResultBackend.
    Extracted from zndraw_joblib/conftest.py for shared use.
    """

    def __init__(self) -> None:
        self._store: dict[str, bytes] = {}
        self._inflight: set[str] = set()
        self._waiters: dict[str, list[asyncio.Event]] = {}

    async def store(self, key: str, data: bytes, _ttl: int) -> None:
        self._store[key] = data
        await self.notify_key(key)

    async def get(self, key: str) -> bytes | None:
        return self._store.get(key)

    async def delete(self, key: str) -> None:
        self._store.pop(key, None)

    async def acquire_inflight(self, key: str, _ttl: int) -> bool:
        if key in self._inflight:
            return False
        self._inflight.add(key)
        return True

    async def release_inflight(self, key: str) -> None:
        self._inflight.discard(key)

    async def wait_for_key(self, key: str, wait_timeout: float) -> bytes | None:
        cached = self._store.get(key)
        if cached is not None:
            return cached
        event = asyncio.Event()
        self._waiters.setdefault(key, []).append(event)
        try:
            await asyncio.wait_for(event.wait(), timeout=wait_timeout)
            return self._store.get(key)
        except TimeoutError:
            return None
        finally:
            waiters = self._waiters.get(key, [])
            if event in waiters:
                waiters.remove(event)

    async def notify_key(self, key: str) -> None:
        for event in self._waiters.pop(key, []):
            event.set()
