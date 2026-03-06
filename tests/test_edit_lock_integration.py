"""Integration tests for edit lock behavior via the ZnDraw Python client.

Tests that lock acquisition blocks other clients from modifying,
and that the lock holder can proceed.
"""

import uuid

import ase
import pytest

from zndraw import ZnDraw
from zndraw.client import RoomLockedError


def _make_atoms() -> ase.Atoms:
    return ase.Atoms("H", positions=[[0, 0, 0]])


def test_lock_blocks_other_client(server: str):
    """Lock held by Client A blocks modifications from Client B."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    a.append(_make_atoms())

    with a.get_lock(msg="editing"), pytest.raises(RoomLockedError):
        b.append(_make_atoms())

    a.disconnect()
    b.disconnect()


def test_lock_holder_can_modify(server: str):
    """Lock holder can still modify the room."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)

    a.append(_make_atoms())

    with a.get_lock(msg="editing"):
        a.append(_make_atoms())
        assert len(a) == 2

    a.disconnect()


def test_lock_release_unblocks(server: str):
    """After lock release, other clients can modify again."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    a.append(_make_atoms())

    with a.get_lock(msg="editing"):
        pass  # Lock released after context

    # B should now be able to modify
    b.append(_make_atoms())
    assert len(b) == 2

    a.disconnect()
    b.disconnect()


def test_multiple_operations_under_lock(server: str):
    """Multiple operations by the lock holder succeed."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)

    with a.get_lock(msg="batch edit"):
        for i in range(5):
            atoms = ase.Atoms("H", positions=[[float(i), 0, 0]])
            a.append(atoms)

    assert len(a) == 5

    a.disconnect()
