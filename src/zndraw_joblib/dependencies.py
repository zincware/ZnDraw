# src/zndraw_joblib/dependencies.py
import hashlib
import json
from typing import Annotated, Any, Protocol, runtime_checkable

from fastapi import Depends, Path, Request
from fastapi_users.authentication import JWTStrategy
from sqlmodel import select
from zndraw_socketio import AsyncServerWrapper

from zndraw_auth import User
from zndraw_auth.db import SessionDep
from zndraw_joblib.exceptions import InvalidRoomId
from zndraw_joblib.registry import InternalRegistry
from zndraw_joblib.settings import JobLibSettings


def get_joblib_settings(request: Request) -> JobLibSettings:
    """Return joblib settings from app.state."""
    return request.app.state.joblib_settings


JobLibSettingsDep = Annotated[JobLibSettings, Depends(get_joblib_settings)]


async def get_internal_registry(request: Request) -> InternalRegistry | None:
    """Return the internal registry from app.state, or None if not configured."""
    return getattr(request.app.state, "internal_registry", None)


def get_tsio(request: Request) -> AsyncServerWrapper | None:
    """Return the Socket.IO server wrapper from app.state.

    Returns None if tsio is not configured, which disables
    all real-time event emissions.
    """
    return getattr(request.app.state, "tsio", None)


def validate_room_id(room_id: str) -> None:
    """Validate room_id doesn't contain @ or : (except @global and @internal)."""
    if room_id in ("@global", "@internal"):
        return
    if not room_id or "@" in room_id or ":" in room_id:
        raise InvalidRoomId.exception(
            detail=f"Room ID '{room_id}' contains invalid characters (@ or :)"
        )


async def verify_writable_room(room_id: str = Path()) -> str:
    """Verify the room is writable. Override in host app for lock checks.

    Default implementation only validates the room_id format.
    Host apps override this via ``app.dependency_overrides[verify_writable_room]``
    to add lock checks (e.g. admin lock, edit lock).
    """
    validate_room_id(room_id)
    return room_id


WritableRoomDep = Annotated[str, Depends(verify_writable_room)]


@runtime_checkable
class ResultBackend(Protocol):
    """Protocol for storing and retrieving cached provider results."""

    async def store(self, key: str, data: bytes, ttl: int) -> None: ...

    async def get(self, key: str) -> bytes | None: ...

    async def delete(self, key: str) -> None: ...

    async def acquire_inflight(self, key: str, ttl: int) -> bool:
        """Return True if lock acquired (SET NX semantics)."""
        ...

    async def release_inflight(self, key: str) -> None: ...

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:  # noqa: ASYNC109
        """Wait for a cache key to be populated.

        Subscribes to a notification channel, checks cache (race-safe),
        then awaits notification or timeout.

        Returns
        -------
        bytes | None
            Cached data if it arrives within timeout, None otherwise.
        """
        ...

    async def notify_key(self, key: str) -> None:
        """Notify waiters that a cache key has been populated."""
        ...


async def get_result_backend() -> ResultBackend:
    """Return the configured ResultBackend.

    Host apps must override this dependency via
    ``app.dependency_overrides[get_result_backend]``.
    """
    raise NotImplementedError(
        "ResultBackend not configured — host app must override get_result_backend"
    )


ResultBackendDep = Annotated[ResultBackend, Depends(get_result_backend)]


@runtime_checkable
class FrameRoomCleanup(Protocol):
    """Protocol for cleaning up frame state when providers are removed."""

    async def __call__(self, room_ids: set[str]) -> None:
        """Clean up external frame state (e.g. Redis keys) for given room IDs."""
        ...


async def _noop_frame_cleanup(room_ids: set[str]) -> None:
    """Default no-op frame cleanup."""


async def get_frame_room_cleanup() -> FrameRoomCleanup:
    """Return the configured FrameRoomCleanup callback.

    Host apps should override this dependency via
    ``app.dependency_overrides[get_frame_room_cleanup]``
    to add Redis key cleanup and FramesInvalidate emission.
    """
    return _noop_frame_cleanup  # type: ignore[return-value]


FrameRoomCleanupDep = Annotated[FrameRoomCleanup, Depends(get_frame_room_cleanup)]


def request_hash(params: dict[str, Any]) -> str:
    """Return a SHA-256 hex digest of the canonicalized JSON representation."""
    canonical = json.dumps(params, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


async def get_worker_token(request: Request, session: SessionDep) -> str:
    """Return a fresh JWT for the internal worker user.

    Reads ``settings`` and ``auth_settings`` from ``request.app.state``
    (set in the host app lifespan).  Accepts ``session`` via DI so
    FastAPI reuses the request-scoped session (avoids SQLite deadlock).
    """
    settings = request.app.state.settings
    auth_settings = request.app.state.auth_settings
    result = await session.exec(
        select(User).where(User.email == settings.internal_worker_email)  # type: ignore[arg-type]
    )
    user = result.one_or_none()
    if user is None:
        raise RuntimeError(
            f"Internal worker user '{settings.internal_worker_email}' not found. "
            "Has the database been initialized?"
        )
    strategy = JWTStrategy(
        secret=auth_settings.secret_key.get_secret_value(),
        lifetime_seconds=auth_settings.token_lifetime_seconds,
    )
    return await strategy.write_token(user)


WorkerTokenDep = Annotated[str, Depends(get_worker_token)]
