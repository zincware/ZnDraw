"""Progress tracking REST API endpoints.

Ephemeral progress trackers stored in Redis (hash), broadcast via Socket.IO.
"""

import json

from fastapi import APIRouter, Response, status

from zndraw.dependencies import (
    CurrentUserDep,
    RedisDep,
    SessionDep,
    SioDep,
    room_channel,
    verify_room,
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

router = APIRouter(prefix="/v1/rooms/{room_id}/progress", tags=["progress"])


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def create_progress(
    session: SessionDep,
    sio: SioDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
    request: ProgressCreate,
) -> ProgressResponse:
    """Create a new progress tracker in the room."""
    await verify_room(session, room_id)

    tracker = ProgressResponse(
        progress_id=request.progress_id,
        description=request.description,
    )
    await redis.hset(  # type: ignore[misc]
        RedisKey.room_progress(room_id),
        request.progress_id,
        tracker.model_dump_json(),
    )

    await sio.emit(
        ProgressStart(
            progress_id=request.progress_id,
            description=request.description,
        ),
        room=room_channel(room_id),
    )

    return tracker


@router.patch(
    "/{progress_id}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, ProgressNotFound),
)
async def update_progress(
    session: SessionDep,
    sio: SioDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
    progress_id: str,
    request: ProgressPatch,
) -> ProgressResponse:
    """Update an existing progress tracker."""
    await verify_room(session, room_id)

    raw = await redis.hget(RedisKey.room_progress(room_id), progress_id)  # type: ignore[misc]
    if raw is None:
        raise ProgressNotFound.exception(f"Progress tracker {progress_id} not found")

    current = json.loads(raw)
    if request.description is not None:
        current["description"] = request.description
    if request.progress is not None:
        current["progress"] = request.progress

    await redis.hset(  # type: ignore[misc]
        RedisKey.room_progress(room_id),
        progress_id,
        json.dumps(current),
    )

    await sio.emit(
        ProgressUpdate(
            progress_id=progress_id,
            description=request.description,
            progress=request.progress,
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
    session: SessionDep,
    sio: SioDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
    progress_id: str,
) -> Response:
    """Complete and remove a progress tracker."""
    await verify_room(session, room_id)

    deleted = await redis.hdel(RedisKey.room_progress(room_id), progress_id)  # type: ignore[misc]
    if not deleted:
        raise ProgressNotFound.exception(f"Progress tracker {progress_id} not found")

    await sio.emit(
        ProgressComplete(progress_id=progress_id),
        room=room_channel(room_id),
    )

    return Response(status_code=status.HTTP_204_NO_CONTENT)
