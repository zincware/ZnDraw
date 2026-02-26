"""Tests for multi-client synchronization via the ZnDraw Python client.

Two clients connected to the same room should see each other's changes
via REST API (no Socket.IO event propagation tested here — just data consistency).
"""

import uuid

import ase
import numpy as np

from zndraw import ZnDraw


def _make_atoms(x: float) -> ase.Atoms:
    """Create a simple H atom at position (x, 0, 0) with info tag."""
    atoms = ase.Atoms("H", positions=[[x, 0, 0]])
    atoms.info["x"] = x
    return atoms


def test_two_clients_see_appends(server: str):
    """Client B sees frames appended by Client A."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    a.append(_make_atoms(1.0))
    a.append(_make_atoms(2.0))

    assert len(b) == 2
    np.testing.assert_allclose(b[0].positions[0, 0], 1.0, atol=1e-6)
    np.testing.assert_allclose(b[1].positions[0, 0], 2.0, atol=1e-6)

    a.disconnect()
    b.disconnect()


def test_two_clients_see_extend(server: str):
    """Client B sees frames extended by Client A."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    frames = [_make_atoms(float(i)) for i in range(5)]
    a.extend(frames)

    assert len(b) == 5
    for i in range(5):
        np.testing.assert_allclose(b[i].positions[0, 0], float(i), atol=1e-6)

    a.disconnect()
    b.disconnect()


def test_two_clients_see_setitem(server: str):
    """Client B sees frame updates by Client A."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    a.extend([_make_atoms(0.0), _make_atoms(1.0)])

    replacement = _make_atoms(99.0)
    a[0] = replacement

    np.testing.assert_allclose(b[0].positions[0, 0], 99.0, atol=1e-6)

    a.disconnect()
    b.disconnect()


def test_two_clients_see_delete(server: str):
    """Client B sees frame deletion by Client A."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    a.extend([_make_atoms(float(i)) for i in range(5)])
    del a[2]

    assert len(b) == 4
    # Frame at index 2 should now be what was index 3
    np.testing.assert_allclose(b[2].positions[0, 0], 3.0, atol=1e-6)

    a.disconnect()
    b.disconnect()


def test_step_sync(server: str):
    """Step set by Client A is visible to Client B."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    a.extend([_make_atoms(float(i)) for i in range(5)])
    a.step = 3

    assert b.step == 3

    a.disconnect()
    b.disconnect()


def test_bookmark_sync(server: str):
    """Bookmarks set by Client A are visible to Client B."""
    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    a.append(_make_atoms(0.0))
    a.bookmarks[0] = "origin"

    assert b.bookmarks[0] == "origin"

    a.disconnect()
    b.disconnect()


def test_geometry_sync(server: str):
    """Geometries set by Client A are visible to Client B."""
    from zndraw.geometries import Sphere

    room_id = uuid.uuid4().hex
    a = ZnDraw(url=server, room=room_id)
    b = ZnDraw(url=server, room=room_id)

    a.geometries["sphere1"] = Sphere(radius=[2.0])

    retrieved = b.geometries["sphere1"]
    assert isinstance(retrieved, Sphere)

    a.disconnect()
    b.disconnect()
