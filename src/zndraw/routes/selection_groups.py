"""Selection Groups REST API endpoints."""

import json

from fastapi import APIRouter

from zndraw.dependencies import (
    CurrentUserDep,
    SessionDep,
    SioDep,
    WritableRoomDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    SelectionGroupNotFound,
    problem_responses,
)
from zndraw.models import SelectionGroup
from zndraw.schemas import (
    SelectionGroupResponse,
    SelectionGroupsListResponse,
    SelectionGroupUpdateRequest,
    StatusResponse,
)
from zndraw.socket_events import SelectionGroupsInvalidate

router = APIRouter(
    prefix="/v1/rooms/{room_id}/selection-groups", tags=["selection-groups"]
)


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_selection_groups(
    session: SessionDep,
    current_user: CurrentUserDep,
    room_id: str,
) -> SelectionGroupsListResponse:
    """List all selection groups for a room."""
    await verify_room(session, room_id)

    from sqlmodel import select

    result = await session.execute(
        select(SelectionGroup).where(SelectionGroup.room_id == room_id)
    )
    groups: dict[str, dict[str, list[int]]] = {}
    for row in result.scalars().all():
        groups[row.name] = json.loads(row.selections)

    return SelectionGroupsListResponse(items=groups)


@router.get(
    "/{group_name}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, SelectionGroupNotFound),
)
async def get_selection_group(
    session: SessionDep,
    current_user: CurrentUserDep,
    room_id: str,
    group_name: str,
) -> SelectionGroupResponse:
    """Get a selection group by name."""
    await verify_room(session, room_id)
    row = await session.get(SelectionGroup, (room_id, group_name))
    if row is None:
        raise SelectionGroupNotFound.exception(
            f"Selection group '{group_name}' not found"
        )
    return SelectionGroupResponse(group=json.loads(row.selections))


@router.put(
    "/{group_name}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, RoomLocked),
)
async def update_selection_group(
    session: SessionDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    group_name: str,
    request: SelectionGroupUpdateRequest,
) -> StatusResponse:
    """Create or update a selection group."""
    row = await session.get(SelectionGroup, (room_id, group_name))
    if row is None:
        row = SelectionGroup(
            room_id=room_id,
            name=group_name,
            selections=json.dumps(request.selections),
        )
        session.add(row)
    else:
        row.selections = json.dumps(request.selections)
    await session.commit()

    await sio.emit(
        SelectionGroupsInvalidate(room_id=room_id), room=room_channel(room_id)
    )
    return StatusResponse()


@router.delete(
    "/{group_name}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, SelectionGroupNotFound
    ),
)
async def delete_selection_group(
    session: SessionDep,
    sio: SioDep,
    room: WritableRoomDep,
    room_id: str,
    group_name: str,
) -> StatusResponse:
    """Delete a selection group."""
    row = await session.get(SelectionGroup, (room_id, group_name))
    if row is None:
        raise SelectionGroupNotFound.exception(
            f"Selection group '{group_name}' not found"
        )
    await session.delete(row)
    await session.commit()

    await sio.emit(
        SelectionGroupsInvalidate(room_id=room_id), room=room_channel(room_id)
    )

    return StatusResponse()
