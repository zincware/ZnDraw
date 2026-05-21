"""Bookmarks REST API endpoints for room frame bookmarks."""

from uuid import UUID

from fastapi import APIRouter
from sqlmodel import select

from zndraw.dependencies import (
    AccessReadDep,
    SessionDep,
    SioDep,
    WritableRoomDep,
    room_channel,
)
from zndraw.exceptions import (
    BookmarkNotFound,
    NotAuthenticated,
    RoomLocked,
    RoomNotFound,
    problem_responses,
)
from zndraw.models import RoomBookmark
from zndraw.schemas import (
    BookmarkCreateRequest,
    BookmarkResponse,
    BookmarksResponse,
    StatusResponse,
)
from zndraw.socket_events import BookmarksInvalidate

router = APIRouter(prefix="/v1/rooms/{owner_id}/{room_name}/bookmarks", tags=["bookmarks"])


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_bookmarks(
    session: SessionDep,
    access: AccessReadDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
) -> BookmarksResponse:
    """Get all bookmarks for a room."""
    room_id = access.room.id
    result = await session.exec(
        select(RoomBookmark).where(RoomBookmark.room_id == room_id)
    )
    rows = result.all()
    bookmarks = {str(row.frame_index): row.label for row in rows}
    return BookmarksResponse(items=bookmarks)


@router.get(
    "/{index}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, BookmarkNotFound),
)
async def get_bookmark(
    session: SessionDep,
    access: AccessReadDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
    index: int,
) -> BookmarkResponse:
    """Get a single bookmark by frame index."""
    room_id = access.room.id
    row = await session.get(RoomBookmark, (room_id, index))
    if row is None:
        raise BookmarkNotFound.exception(f"Bookmark '{index}' not found")
    return BookmarkResponse(index=index, label=row.label)


@router.put(
    "/{index}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, RoomLocked),
)
async def set_bookmark(
    session: SessionDep,
    sio: SioDep,
    _room: WritableRoomDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
    index: int,
    request: BookmarkCreateRequest,
) -> BookmarkResponse:
    """Create or update a bookmark."""
    room_id = _room.id

    row = await session.get(RoomBookmark, (room_id, index))
    if row is None:
        row = RoomBookmark(room_id=room_id, frame_index=index, label=request.label)
        session.add(row)
    else:
        row.label = request.label
    await session.commit()

    await sio.emit(
        BookmarksInvalidate(room_id=room_id, index=index, operation="set"),
        room=room_channel(room_id),
    )

    return BookmarkResponse(index=index, label=request.label)


@router.delete(
    "/{index}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, RoomLocked, BookmarkNotFound
    ),
)
async def delete_bookmark(
    session: SessionDep,
    sio: SioDep,
    _room: WritableRoomDep,
    owner_id: UUID,  # noqa: ARG001
    room_name: str,  # noqa: ARG001
    index: int,
) -> StatusResponse:
    """Delete a bookmark."""
    room_id = _room.id

    row = await session.get(RoomBookmark, (room_id, index))
    if row is None:
        raise BookmarkNotFound.exception(f"Bookmark '{index}' not found")
    await session.delete(row)
    await session.commit()

    await sio.emit(
        BookmarksInvalidate(room_id=room_id, index=index, operation="delete"),
        room=room_channel(room_id),
    )

    return StatusResponse()
