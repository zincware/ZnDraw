"""Unit tests for the Room namespace data model."""

from __future__ import annotations

import pytest
import pytest_asyncio
from sqlalchemy.exc import IntegrityError
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw.access import Visibility
from zndraw.models import Group, Room
from zndraw_auth import User


@pytest_asyncio.fixture
async def alice(session: AsyncSession) -> User:
    user = User(email="alice@example.com", hashed_password="x")
    session.add(user)
    await session.commit()
    await session.refresh(user)
    return user


async def test_room_public_address_user_owned(
    session: AsyncSession, alice: User
) -> None:
    room = Room(
        room_name="my-experiment",
        owner_user_id=alice.id,
        visibility=Visibility.PRIVATE,
        created_by_id=alice.id,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)

    assert room.public_address == f"{alice.id}/my-experiment"


async def test_room_public_address_group_owned(
    session: AsyncSession, alice: User
) -> None:
    group = Group(name="g1", created_by_id=alice.id)
    session.add(group)
    await session.commit()
    await session.refresh(group)

    room = Room(
        room_name="shared",
        owner_group_id=group.id,
        visibility=Visibility.GROUP,
        created_by_id=alice.id,
    )
    session.add(room)
    await session.commit()
    await session.refresh(room)

    assert room.public_address == f"{group.id}/shared"


async def test_unique_per_owner(session: AsyncSession, alice: User) -> None:
    session.add(
        Room(
            room_name="dup",
            owner_user_id=alice.id,
            visibility=Visibility.PRIVATE,
            created_by_id=alice.id,
        )
    )
    await session.commit()

    session.add(
        Room(
            room_name="dup",
            owner_user_id=alice.id,
            visibility=Visibility.PRIVATE,
            created_by_id=alice.id,
        )
    )
    with pytest.raises(IntegrityError):
        await session.commit()


async def test_same_name_different_owners(session: AsyncSession, alice: User) -> None:
    bob = User(email="bob@example.com", hashed_password="x")
    session.add(bob)
    await session.commit()
    await session.refresh(bob)

    session.add(
        Room(
            room_name="r",
            owner_user_id=alice.id,
            visibility=Visibility.PRIVATE,
            created_by_id=alice.id,
        )
    )
    session.add(
        Room(
            room_name="r",
            owner_user_id=bob.id,
            visibility=Visibility.PRIVATE,
            created_by_id=bob.id,
        )
    )
    await session.commit()  # no IntegrityError
