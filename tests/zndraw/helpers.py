"""Shared test utility functions and classes.

These are importable helpers (not pytest fixtures).
Fixtures live in conftest.py.
"""

from typing import Any

import msgpack
from fastapi_users.password import PasswordHelper
from pydantic import BaseModel
from sqlalchemy.ext.asyncio import AsyncSession
from zndraw_auth import User
from zndraw_auth.settings import AuthSettings

from zndraw.models import Room
from zndraw.storage import RawFrame

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
    session: AsyncSession, email: str = "testuser@local.test"
) -> tuple[User, str]:
    """Create a user in the DB and return (user, token)."""
    user = create_test_user_model(email=email)
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user, create_test_token(user)


async def create_test_room(
    session: AsyncSession, user: User, description: str = "Test Room"
) -> Room:
    """Create a room with user as owner and return it."""
    from zndraw.models import MemberRole, RoomMembership

    room = Room(
        description=description,
        created_by_id=user.id,  # type: ignore[arg-type]
        is_public=True,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)

    membership = RoomMembership(
        room_id=room.id,  # type: ignore[arg-type]
        user_id=user.id,  # type: ignore[arg-type]
        role=MemberRole.OWNER,
    )
    session.add(membership)
    await session.commit()
    return room


def auth_header(token: str) -> dict[str, str]:
    """Return Authorization header dict."""
    return {"Authorization": f"Bearer {token}"}


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
