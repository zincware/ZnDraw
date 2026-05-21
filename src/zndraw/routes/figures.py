"""Figures REST API endpoints for room Plotly figures."""

from fastapi import APIRouter, status
from sqlmodel import select

from zndraw.dependencies import (
    AccessReadDep,
    SessionDep,
    SioDep,
    WritableRoomDep,
    room_channel,
)
from zndraw.exceptions import (
    FigureNotFound,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    problem_responses,
)
from zndraw.models import RoomFigure
from zndraw.schemas import (
    CollectionResponse,
    FigureCreateRequest,
    FigureCreateResponse,
    FigureData,
    FigureResponse,
    StatusResponse,
)
from zndraw.socket_events import FigureInvalidate

router = APIRouter(prefix="/v1/rooms/{owner_id}/{room_name}/figures", tags=["figures"])


@router.get(
    "",
    responses=problem_responses(RoomNotFound),
)
async def list_figures(
    session: SessionDep,
    access: AccessReadDep,
) -> CollectionResponse[str]:
    """List all figure keys in a room."""
    room_id = access.room.id
    result = await session.exec(
        select(RoomFigure.key).where(RoomFigure.room_id == room_id)
    )
    keys = list(result.all())
    return CollectionResponse(items=keys)


@router.get(
    "/{key}",
    responses=problem_responses(RoomNotFound, FigureNotFound),
)
async def get_figure(
    session: SessionDep,
    access: AccessReadDep,
    key: str,
) -> FigureResponse:
    """Get a single figure by key."""
    room_id = access.room.id
    row = await session.get(RoomFigure, (room_id, key))
    if row is None:
        raise FigureNotFound.exception(f"Figure '{key}' not found")
    return FigureResponse(key=key, figure=FigureData(type=row.type, data=row.data))


@router.post(
    "/{key}",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomNotFound, RoomLocked),
)
async def create_figure(
    session: SessionDep,
    sio: SioDep,
    _room: WritableRoomDep,
    key: str,
    request: FigureCreateRequest,
) -> FigureCreateResponse:
    """Create or update a figure."""
    room_id = _room.id

    row = await session.get(RoomFigure, (room_id, key))
    created = row is None
    if created:
        row = RoomFigure(
            room_id=room_id, key=key, type=request.figure.type, data=request.figure.data
        )
        session.add(row)
    else:
        row.type = request.figure.type
        row.data = request.figure.data
    await session.commit()

    await sio.emit(
        FigureInvalidate(room_id=room_id, key=key, operation="set"),
        room=room_channel(room_id),
    )

    return FigureCreateResponse(key=key, created=created)


@router.delete(
    "/{key}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, FigureNotFound
    ),
)
async def delete_figure(
    session: SessionDep,
    sio: SioDep,
    _room: WritableRoomDep,
    key: str,
) -> StatusResponse:
    """Delete a figure."""
    room_id = _room.id

    row = await session.get(RoomFigure, (room_id, key))
    if row is None:
        raise FigureNotFound.exception(f"Figure '{key}' not found")
    await session.delete(row)
    await session.commit()

    await sio.emit(
        FigureInvalidate(room_id=room_id, key=key, operation="delete"),
        room=room_channel(room_id),
    )

    return StatusResponse()
