"""Shared room-lookup helpers for joblib emit machinery.

Both ``router.py`` and ``sweeper.py`` need to map a joblib ``room_id``
(real UUID, composed ``<uuid>/<name>``, ``@global``, ``@internal``, or
unknown garbage from the wire) to either the persisted ``Room`` row or
the address that should be surfaced on the wire when no row exists.

Joblib-only test environments may not have the ``Room`` table — the
``OperationalError`` / ``ProgrammingError`` catch keeps these helpers
usable in those harnesses.
"""

from __future__ import annotations

from typing import TYPE_CHECKING
from uuid import UUID

from sqlalchemy.exc import OperationalError, ProgrammingError

if TYPE_CHECKING:
    from sqlmodel.ext.asyncio.session import AsyncSession

    from zndraw.models import Room


async def fetch_room(session: AsyncSession, room_id: str) -> Room | None:
    """Return the Room for ``room_id``, or None.

    Accepts a UUID string, a composed ``<owner>/<name>`` address (where
    ``<owner>`` is either a UUID or a display-name / group-name slug), or
    a sigil. Returns None for sigils, unknown ids, and missing rooms.
    """
    if room_id in ("@global", "@internal"):
        return None
    from zndraw.models import Room

    try:
        try:
            UUID(room_id)
        except ValueError:
            if "/" not in room_id:
                return None
            owner_part, _, name_part = room_id.partition("/")
            owner_uuid: UUID | None
            try:
                owner_uuid = UUID(owner_part)
            except ValueError:
                # Display-name / group-name slug — resolve to UUID.
                from zndraw.dependencies import get_owner_uuid_from_segment

                try:
                    owner_uuid = await get_owner_uuid_from_segment(
                        session, owner_part
                    )
                except Exception:
                    return None
            from sqlmodel import col, or_, select as sql_select

            result = await session.exec(
                sql_select(Room).where(
                    or_(
                        Room.owner_user_id == owner_uuid,
                        Room.owner_group_id == owner_uuid,
                    ),
                    col(Room.room_name) == name_part,
                )
            )
            return result.one_or_none()
        return await session.get(Room, room_id)
    except (OperationalError, ProgrammingError):
        return None


async def room_address_for(session: AsyncSession, room_id: str) -> str:
    """Return display-name composed address for a known room, else echo ``room_id``."""
    if room_id in ("@global", "@internal"):
        return room_id
    room = await fetch_room(session, room_id)
    if room is None:
        return room_id
    from zndraw.models import build_public_address

    return await build_public_address(session, room)
