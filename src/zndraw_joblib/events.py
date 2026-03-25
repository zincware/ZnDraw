# src/zndraw_joblib/events.py
"""Socket.IO event models for real-time notifications.

All models are frozen for hashability, enabling set-based deduplication
of emissions via the Emission NamedTuple.
"""

from __future__ import annotations

import json
from datetime import datetime
from typing import Any, NamedTuple

from pydantic import BaseModel, ConfigDict
from zndraw_socketio import AsyncServerWrapper

from zndraw_joblib.models import Task, TaskStatus


class FrozenEvent(BaseModel):
    """Base class for all frozen event models.

    Provides frozen=True config for hashability, required by Emission sets.
    """

    model_config = ConfigDict(frozen=True)


class JobsInvalidate(FrozenEvent):
    """Frontend should refetch the job list."""


class TaskAvailable(FrozenEvent):
    """A new task is available for claiming."""

    job_name: str
    room_id: str
    task_id: str


class TaskStatusEvent(FrozenEvent):
    """A task's status changed."""

    id: str
    name: str
    room_id: str
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


class ProvidersInvalidate(FrozenEvent):
    """Frontend should refetch the provider list."""


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
    ) -> "ProviderRequest":
        """Create from a dict, converting params to canonical JSON."""
        return cls(
            request_id=request_id,
            provider_name=provider_name,
            params=json.dumps(params, sort_keys=True, separators=(",", ":")),
        )


class ProviderResultReady(FrozenEvent):
    """Server notifies frontend that a provider result is cached."""

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
    The host app's handler should call ``tsio.leave_room(sid, f"providers:{provider_name}")``.
    """

    provider_name: str  # full_name: room_id:category:name
    worker_id: str


class Emission(NamedTuple):
    """Hashable (event, room) pair for set-based deduplication."""

    event: FrozenEvent
    room: str


def build_task_status_emission(
    task: Task,
    job_full_name: str,
    queue_position: int | None = None,
) -> Emission:
    """Build a TaskStatusEvent emission from task data."""
    return Emission(
        TaskStatusEvent(
            id=str(task.id),
            name=job_full_name,
            room_id=task.room_id,
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


async def emit(tsio: AsyncServerWrapper | None, emissions: set[Emission]) -> None:
    """Emit a set of events via the Socket.IO server wrapper.

    Uses the zndraw-socketio API (passing Pydantic models directly).
    No-op if tsio is None.
    """
    if not tsio:
        return
    for emission in emissions:
        await tsio.emit(emission.event, room=emission.room)
