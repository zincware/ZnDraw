"""E2E test: Client data persists across disconnect/reconnect cycles.

Uses a real uvicorn server via the server_factory fixture.
"""

import uuid

import ase
import numpy as np

from zndraw import ZnDraw


def _make_atoms(x: float) -> ase.Atoms:
    """Create a simple H atom at position (x, 0, 0)."""
    atoms = ase.Atoms("H", positions=[[x, 0, 0]])
    atoms.info["x"] = x
    return atoms


def test_frames_persist_after_disconnect(server: str):
    """Data written by a client survives after it disconnects.

    A second client connecting to the same room should see the frames
    that were written before the first client disconnected.
    """
    room_id = uuid.uuid4().hex

    # First client: write frames then disconnect
    client_a = ZnDraw(url=server, room=room_id)
    client_a.extend([_make_atoms(float(i)) for i in range(3)])
    client_a.disconnect()

    # Second client: reconnect and read frames
    client_b = ZnDraw(url=server, room=room_id)
    assert len(client_b) == 3
    np.testing.assert_allclose(client_b[0].positions[0, 0], 0.0, atol=1e-6)
    np.testing.assert_allclose(client_b[1].positions[0, 0], 1.0, atol=1e-6)
    np.testing.assert_allclose(client_b[2].positions[0, 0], 2.0, atol=1e-6)
    client_b.disconnect()


def test_step_persists_after_disconnect(server: str):
    """Step value set by a client survives after it disconnects."""
    room_id = uuid.uuid4().hex

    client_a = ZnDraw(url=server, room=room_id)
    client_a.extend([_make_atoms(float(i)) for i in range(5)])
    client_a.step = 4
    client_a.disconnect()

    client_b = ZnDraw(url=server, room=room_id)
    assert client_b.step == 4
    client_b.disconnect()


def test_multiple_reconnect_cycles(server: str):
    """Data remains consistent across multiple disconnect/reconnect cycles."""
    room_id = uuid.uuid4().hex

    # First write
    c1 = ZnDraw(url=server, room=room_id)
    c1.append(_make_atoms(10.0))
    c1.disconnect()

    # Reconnect and append more
    c2 = ZnDraw(url=server, room=room_id)
    assert len(c2) == 1
    c2.append(_make_atoms(20.0))
    c2.disconnect()

    # Final reconnect: should see both frames
    c3 = ZnDraw(url=server, room=room_id)
    assert len(c3) == 2
    np.testing.assert_allclose(c3[0].positions[0, 0], 10.0, atol=1e-6)
    np.testing.assert_allclose(c3[1].positions[0, 0], 20.0, atol=1e-6)
    c3.disconnect()


def test_guest_token_reconnect(server: str):
    """Client reconnects with a fresh guest token and still sees the room data."""
    import httpx

    room_id = uuid.uuid4().hex

    # Write with first guest token
    token_resp = httpx.post(f"{server}/v1/auth/guest", timeout=10.0)
    assert token_resp.status_code == 200
    token_a = token_resp.json()["access_token"]

    c1 = ZnDraw(url=server, room=room_id, token=token_a)
    c1.append(_make_atoms(42.0))
    c1.disconnect()

    # Reconnect with a different guest token
    token_resp2 = httpx.post(f"{server}/v1/auth/guest", timeout=10.0)
    assert token_resp2.status_code == 200
    token_b = token_resp2.json()["access_token"]

    c2 = ZnDraw(url=server, room=room_id, token=token_b)
    assert len(c2) == 1
    np.testing.assert_allclose(c2[0].positions[0, 0], 42.0, atol=1e-6)
    c2.disconnect()
