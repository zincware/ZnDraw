"""Chat REST API endpoints for room messages."""

from datetime import UTC, datetime
from typing import Annotated

from fastapi import APIRouter, Query, status
from sqlalchemy import func
from sqlmodel import col, select

from zndraw.broadcast import broadcast_to_room
from zndraw.dependencies import (
    AccessEditDep,
    AccessReadDep,
    CurrentUserDep,
    SessionDep,
    SioDep,
)
from zndraw.exceptions import (
    MessageNotFound,
    NotAuthenticated,
    NotMessageOwner,
    RoomNotFound,
    problem_responses,
)
from zndraw.models import Message
from zndraw.schemas import (
    MessageCreate,
    MessageEditRequest,
    MessageResponse,
    MessagesMetadata,
    MessagesResponse,
)
from zndraw.socket_events import MessageEdited, MessageNew
from zndraw_auth import User

router = APIRouter(
    prefix="/v1/rooms/{owner}/{room_name}/chat/messages", tags=["chat"]
)


def _datetime_to_unix_ms(dt: datetime) -> int:
    """Convert a datetime to unix milliseconds."""
    return int(dt.timestamp() * 1000)


def _message_to_response(
    msg: Message, display_name: str | None = None
) -> MessageResponse:
    """Convert a Message model to MessageResponse."""
    return MessageResponse(
        id=msg.id,  # type: ignore[arg-type]
        room_id=msg.room_id,
        user_id=msg.user_id,
        content=msg.content,
        created_at=msg.created_at,
        updated_at=msg.updated_at,
        display_name=display_name,
    )


@router.get(
    "",
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def list_messages(
    session: SessionDep,
    access: AccessReadDep,
    limit: Annotated[int, Query(ge=1, le=100)] = 30,
    before: Annotated[int | None, Query(description="Unix ms cursor")] = None,
) -> MessagesResponse:
    """List messages with cursor pagination (newest first)."""
    room_id = access.room.id

    stmt = select(Message).where(Message.room_id == room_id)

    if before is not None:
        before_dt = datetime.fromtimestamp(before / 1000, tz=UTC)
        stmt = stmt.where(col(Message.created_at) < before_dt)

    # Fetch one extra to determine has_more
    stmt = stmt.order_by(col(Message.created_at).desc()).limit(limit + 1)
    result = await session.exec(stmt)
    rows = list(result.all())

    has_more = len(rows) > limit
    rows = rows[:limit]

    # Get total count
    count_stmt = (
        select(func.count()).select_from(Message).where(Message.room_id == room_id)
    )
    total_count = (await session.exec(count_stmt)).one()

    # Look up display names for all user_ids
    user_ids = {row.user_id for row in rows}
    display_name_map: dict[str, str | None] = {}
    if user_ids:
        users_result = await session.exec(
            select(User).where(col(User.id).in_(user_ids))
        )
        for user in users_result.all():
            display_name_map[str(user.id)] = user.display_name

    messages = [
        _message_to_response(row, display_name_map.get(str(row.user_id)))
        for row in rows
    ]

    oldest_ts = _datetime_to_unix_ms(rows[-1].created_at) if rows else None
    newest_ts = _datetime_to_unix_ms(rows[0].created_at) if rows else None

    return MessagesResponse(
        items=messages,
        metadata=MessagesMetadata(
            has_more=has_more,
            total_count=total_count,
            oldest_timestamp=oldest_ts,
            newest_timestamp=newest_ts,
        ),
    )


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def create_message(
    session: SessionDep,
    sio: SioDep,
    current_user: CurrentUserDep,
    access: AccessEditDep,
    request: MessageCreate,
) -> MessageResponse:
    """Create a new chat message."""
    room_id = access.room.id

    msg = Message(
        room_id=room_id,
        user_id=current_user.id,  # type: ignore[arg-type]
        content=request.content,
    )
    session.add(msg)
    await session.commit()
    await session.refresh(msg)

    display_name = current_user.display_name

    await broadcast_to_room(
        sio,
        MessageNew.for_room(
            access.room,
            id=msg.id,  # type: ignore[arg-type]
            user_id=current_user.id,  # type: ignore[arg-type]
            content=msg.content,
            created_at=msg.created_at,
            email=display_name,
        ),
        access.room,
    )

    return _message_to_response(msg, display_name)


@router.patch(
    "/{message_id}",
    responses=problem_responses(
        NotAuthenticated, RoomNotFound, MessageNotFound, NotMessageOwner
    ),
)
async def edit_message(
    session: SessionDep,
    sio: SioDep,
    current_user: CurrentUserDep,
    access: AccessReadDep,
    message_id: int,
    request: MessageEditRequest,
) -> MessageResponse:
    """Edit an existing chat message. Only the author can edit."""
    room_id = access.room.id

    msg = await session.get(Message, message_id)
    if msg is None or msg.room_id != room_id:
        raise MessageNotFound.exception(f"Message {message_id} not found")

    if msg.user_id != current_user.id:
        raise NotMessageOwner.exception("Only the message author can edit")

    msg.content = request.content
    msg.updated_at = datetime.now(UTC)
    await session.commit()
    await session.refresh(msg)

    await broadcast_to_room(
        sio,
        MessageEdited.for_room(
            access.room,
            id=msg.id,  # type: ignore[arg-type]
            content=msg.content,
            updated_at=msg.updated_at,
        ),
        access.room,
    )

    return _message_to_response(msg, current_user.display_name)
