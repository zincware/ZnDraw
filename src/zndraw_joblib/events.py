# src/zndraw_joblib/events.py
"""Socket.IO event models for real-time notifications.

All models are frozen for hashability, enabling set-based deduplication
of emissions via the Emission NamedTuple.
"""

from __future__ import annotations

import json
import uuid as _uuid
from datetime import datetime  # noqa: TC003
from typing import TYPE_CHECKING, Any, NamedTuple
from uuid import UUID

from pydantic import BaseModel, ConfigDict

from zndraw.socket_events import RoomScopedEvent

if TYPE_CHECKING:
    from sqlmodel.ext.asyncio.session import AsyncSession
    from zndraw_socketio import AsyncServerWrapper

from zndraw_joblib.models import Task, TaskStatus  # noqa: TC001

NIL_ROOM_UUID = UUID(int=0)

# Stable namespace for deriving a placeholder room UUID from a composed
# address when no persisted Room is available (joblib-only test envs).
_ROOM_ADDRESS_NS = UUID("c4a4f5fd-8b8a-5e7a-9c8e-1f4a2b3c4d5e")


def event_room_uuid(room_id: str) -> UUID:
    """Map a joblib room_id (surrogate UUID, sigil, or composed address) to a UUID.

    For composed ``<owner>/<name>`` addresses with no persisted Room
    (joblib-only test envs), derives a stable uuid5 from the address so
    the ``RoomScopedEvent`` validator accepts the payload.
    """
    if room_id in ("@global", "@internal"):
        return NIL_ROOM_UUID
    try:
        return UUID(room_id)
    except ValueError:
        pass
    if "/" in room_id:
        return _uuid.uuid5(_ROOM_ADDRESS_NS, room_id)
    return NIL_ROOM_UUID


class FrozenEvent(BaseModel):
    """Base class for all frozen event models.

    Provides frozen=True config for hashability, required by Emission sets.
    """

    model_config = ConfigDict(frozen=True)


class JobsInvalidate(RoomScopedEvent):
    """Frontend should refetch the job list."""

    model_config = ConfigDict(frozen=True)


class TaskAvailable(FrozenEvent):
    """A new task is available for claiming."""

    job_name: str
    room_id: str
    task_id: str


class TaskStatusEvent(RoomScopedEvent):
    """A task's status changed."""

    model_config = ConfigDict(frozen=True)

    id: str
    name: str
    status: TaskStatus
    created_at: datetime
    started_at: datetime | None = None
    completed_at: datetime | None = None
    queue_position: int | None = None
    worker_id: str | None = None
    error: str | None = None


class JoinJobRoom(FrozenEvent):
    """Worker requests to join a job's notification room.

    Sent by the client after REST job registration. The host app's
    socketio handler should call ``tsio.enter_room(sid, f"jobs:{job_name}")``
    and store the ``worker_id`` in the SIO session for disconnect cleanup.
    """

    job_name: str
    worker_id: str


class LeaveJobRoom(FrozenEvent):
    """Worker requests to leave a job's notification room.

    Sent by the client on graceful disconnect or job unregistration.
    The host app's handler should call ``tsio.leave_room(sid, f"jobs:{job_name}")``.
    """

    job_name: str
    worker_id: str


class ProvidersInvalidate(RoomScopedEvent):
    """Frontend should refetch the provider list."""

    model_config = ConfigDict(frozen=True)


class ProviderRequest(FrozenEvent):
    """Server dispatches a read request to a provider client.

    ``params`` is stored as a canonical JSON string (sorted keys, compact
    separators) so it is both hashable (frozen model) and directly usable
    by Socket.IO clients without tuple-to-dict conversion.
    """

    request_id: str
    provider_name: str  # full_name: room_id:category:name
    params: str  # canonical JSON string

    @classmethod
    def from_dict_params(
        cls,
        *,
        request_id: str,
        provider_name: str,
        params: dict[str, Any],
    ) -> ProviderRequest:
        """Create from a dict, converting params to canonical JSON."""
        return cls(
            request_id=request_id,
            provider_name=provider_name,
            params=json.dumps(params, sort_keys=True, separators=(",", ":")),
        )


class ProviderResultReady(RoomScopedEvent):
    """Server notifies frontend that a provider result is cached."""

    model_config = ConfigDict(frozen=True)

    provider_name: str  # full_name: room_id:category:name
    request_hash: str


class JoinProviderRoom(FrozenEvent):
    """Client joins a provider dispatch room.

    Sent by the client after REST provider registration. The host app's
    socketio handler should call ``tsio.enter_room(sid, f"providers:{provider_name}")``
    and store the ``worker_id`` in the SIO session for disconnect cleanup.
    """

    provider_name: str  # full_name: room_id:category:name
    worker_id: str


class LeaveProviderRoom(FrozenEvent):
    """Client leaves a provider dispatch room.

    Sent by the client on graceful disconnect or provider unregistration.
    The host app's handler should call
    ``tsio.leave_room(sid, f"providers:{provider_name}")``.
    """

    provider_name: str  # full_name: room_id:category:name
    worker_id: str


class Emission(NamedTuple):
    """Hashable (event, room) pair for set-based deduplication."""

    event: BaseModel
    room: str


def build_task_status_emission(
    task: Task,
    job_full_name: str,
    room_address: str,
    queue_position: int | None = None,
) -> Emission:
    """Build a TaskStatusEvent emission from task data."""
    return Emission(
        TaskStatusEvent(
            room_id=event_room_uuid(task.room_id),
            room_address=room_address,
            id=str(task.id),
            name=job_full_name,
            status=task.status,
            created_at=task.created_at,
            started_at=task.started_at,
            completed_at=task.completed_at,
            queue_position=queue_position,
            worker_id=str(task.worker_id) if task.worker_id else None,
            error=task.error,
        ),
        f"room:{task.room_id}",
    )


async def build_room_scoped_emission(
    session: AsyncSession,
    event_cls: type[RoomScopedEvent],
    room_id: str,
    **fields: Any,
) -> Emission:
    """Construct an Emission for a joblib room-scoped event.

    Handles the sigil/real-room/unknown split centrally:
      * for persisted rooms, defers to ``event_cls.for_room(room, **fields)``;
      * for sigils (``@global``/``@internal``), uses ``NIL_ROOM_UUID`` —
        the ``RoomScopedEvent`` validator allows this because sigil
        addresses do not look composed;
      * for a bare UUID with no persisted row (e.g. soft-deleted room),
        uses the UUID directly as ``room_id``;
      * for a composed ``<owner>/<name>`` address with no persisted row
        (joblib-only test environments), derives a stable
        uuid5 from the address so the validator passes.

    The channel is always ``f"room:{room_id}"``.
    """
    from zndraw_joblib.room_lookup import fetch_room

    room = await fetch_room(session, room_id)
    if room is not None:
        from zndraw.models import build_public_address

        room_address = await build_public_address(session, room)
        return Emission(
            event_cls.for_room(room, room_address=room_address, **fields),
            f"room:{room.id}",
        )
    return Emission(
        event_cls(room_id=event_room_uuid(room_id), room_address=room_id, **fields),
        f"room:{room_id}",
    )


async def emit(tsio: AsyncServerWrapper | None, emissions: set[Emission]) -> None:
    """Emit a set of events via the Socket.IO server wrapper.

    Uses the zndraw-socketio API (passing Pydantic models directly).
    No-op if tsio is None.
    """
    if not tsio:
        return
    for emission in emissions:
        await tsio.emit(emission.event, room=emission.room)
