"""Progress tracking REST API endpoints.

Ephemeral progress trackers stored in Redis (hash), broadcast via Socket.IO.
"""

import json
from uuid import UUID

from fastapi import APIRouter, Response, status

from zndraw.dependencies import (
    AccessEditDep,
    RedisDep,
    SioDep,
    room_channel,
)
from zndraw.exceptions import (
    NotAuthenticated,
    ProgressNotFound,
    RoomNotFound,
    problem_responses,
)
from zndraw.redis import RedisKey
from zndraw.schemas import ProgressCreate, ProgressPatch, ProgressResponse
from zndraw.socket_events import ProgressComplete, ProgressStart, ProgressUpdate

router = APIRouter(prefix="/v1/rooms/{owner_id}/{room_name}/progress", tags=["progress"])

PROGRESS_TTL = 3600  # 1 hour — auto-cleanup for orphaned trackers


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def create_progress(
    sio: SioDep,
    redis: RedisDep,
    access: AccessEditDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
    request: ProgressCreate,
) -> ProgressResponse:
    """Create a new progress tracker in the room."""
    room_id = access.room.id

    tracker = ProgressResponse(
        progress_id=request.progress_id,
        description=request.description,
        unit=request.unit,
    )
    key = RedisKey.room_progress(room_id)
    await redis.hset(key, request.progress_id, tracker.model_dump_json())  # type: ignore[misc]
    await redis.expire(key, PROGRESS_TTL)  # type: ignore[misc]

    await sio.emit(
        ProgressStart(
            progress_id=request.progress_id,
            description=request.description,
            unit=request.unit,
        ),
        room=room_channel(room_id),
    )

    return tracker


@router.patch(
    "/{progress_id}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, ProgressNotFound),
)
async def update_progress(
    sio: SioDep,
    redis: RedisDep,
    access: AccessEditDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
    progress_id: str,
    request: ProgressPatch,
) -> ProgressResponse:
    """Update an existing progress tracker."""
    room_id = access.room.id

    key = RedisKey.room_progress(room_id)
    raw = await redis.hget(key, progress_id)  # type: ignore[misc]
    if raw is None:
        raise ProgressNotFound.exception(f"Progress tracker {progress_id} not found")

    current = json.loads(raw)
    for field in ("description", "n", "total", "elapsed", "unit"):
        value = getattr(request, field)
        if value is not None:
            current[field] = value

    await redis.hset(key, progress_id, json.dumps(current))  # type: ignore[misc]
    await redis.expire(key, PROGRESS_TTL)  # type: ignore[misc]

    await sio.emit(
        ProgressUpdate(
            progress_id=progress_id,
            description=request.description,
            n=request.n,
            total=request.total,
            elapsed=request.elapsed,
            unit=request.unit,
        ),
        room=room_channel(room_id),
    )

    return ProgressResponse(**current)


@router.delete(
    "/{progress_id}",
    status_code=status.HTTP_204_NO_CONTENT,
    responses=problem_responses(NotAuthenticated, RoomNotFound, ProgressNotFound),
)
async def delete_progress(
    sio: SioDep,
    redis: RedisDep,
    access: AccessEditDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
    progress_id: str,
) -> Response:
    """Complete and remove a progress tracker."""
    room_id = access.room.id

    deleted = await redis.hdel(RedisKey.room_progress(room_id), progress_id)  # type: ignore[misc]
    if not deleted:
        raise ProgressNotFound.exception(f"Progress tracker {progress_id} not found")

    await sio.emit(
        ProgressComplete(progress_id=progress_id),
        room=room_channel(room_id),
    )

    return Response(status_code=status.HTTP_204_NO_CONTENT)
