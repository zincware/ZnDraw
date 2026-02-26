"""Step REST API endpoints for room frame navigation."""

from fastapi import APIRouter

from zndraw.dependencies import (
    CurrentUserDep,
    SessionDep,
    SioDep,
    StorageDep,
    WritableRoomDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    StepOutOfBounds,
    problem_responses,
)
from zndraw.schemas import StepResponse, StepUpdateRequest, StepUpdateResponse
from zndraw.socket_events import FrameUpdate

router = APIRouter(prefix="/v1/rooms/{room_id}/step", tags=["step"])


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def get_step(
    session: SessionDep,
    storage: StorageDep,
    current_user: CurrentUserDep,
    room_id: str,
) -> StepResponse:
    """Get current step (frame index) for a room.

    Returns current step and total frame count. Clamps step to valid range
    if frames were deleted.
    """
    room = await verify_room(session, room_id)
    step = room.step
    total = await storage.get_length(room_id)

    # Clamp step to valid range (frames may have been deleted)
    if total > 0 and step >= total:
        step = max(0, total - 1)
        room.step = step
        await session.commit()

    return StepResponse(step=step, total_frames=total)


@router.put(
    "",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, StepOutOfBounds
    ),
)
async def set_step(
    session: SessionDep,
    storage: StorageDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    request: StepUpdateRequest,
) -> StepUpdateResponse:
    """Set current step (frame index) for a room.

    Broadcasts frame:update to the room.
    """
    total = await storage.get_length(room_id)

    if request.step >= total:
        raise StepOutOfBounds.exception(
            f"Step {request.step} out of bounds (0-{total - 1})"
            if total > 0
            else f"Step {request.step} out of bounds (no frames)"
        )

    room.step = request.step
    await session.commit()

    # Broadcast frame update
    await sio.emit(
        FrameUpdate(room_id=room_id, frame=request.step),
        room=room_channel(room_id),
    )

    return StepUpdateResponse(success=True, step=request.step)
