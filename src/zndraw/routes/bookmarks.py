"""Bookmarks REST API endpoints for room frame bookmarks."""

from fastapi import APIRouter, status
from sqlmodel import select

from zndraw.dependencies import (
    CurrentUserDep,
    SessionDep,
    SioDep,
    WritableRoomDep,
    room_channel,
    verify_room,
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

router = APIRouter(prefix="/v1/rooms/{room_id}/bookmarks", tags=["bookmarks"])


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_bookmarks(
    session: SessionDep,
    current_user: CurrentUserDep,
    room_id: str,
) -> BookmarksResponse:
    """Get all bookmarks for a room."""
    await verify_room(session, room_id)
    result = await session.execute(
        select(RoomBookmark).where(RoomBookmark.room_id == room_id)
    )
    rows = result.scalars().all()
    bookmarks = {str(row.frame_index): row.label for row in rows}
    return BookmarksResponse(items=bookmarks)


@router.get(
    "/{index}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, BookmarkNotFound),
)
async def get_bookmark(
    session: SessionDep,
    current_user: CurrentUserDep,
    room_id: str,
    index: int,
) -> BookmarkResponse:
    """Get a single bookmark by frame index."""
    await verify_room(session, room_id)
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
    room: WritableRoomDep,
    room_id: str,
    index: int,
    request: BookmarkCreateRequest,
) -> BookmarkResponse:
    """Create or update a bookmark."""

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
    room: WritableRoomDep,
    room_id: str,
    index: int,
) -> StatusResponse:
    """Delete a bookmark."""

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
