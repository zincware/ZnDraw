"""Integration test for template room isolation.

Verifies that when a room is set as the default template:
1. New rooms copy frames from the template
2. Modifications to new rooms don't affect the template
3. The template remains unchanged after operations on derived rooms
"""

import pathlib
import uuid
from collections.abc import Callable
from typing import Any

import ase
import httpx
import numpy as np
import pytest

from zndraw import ZnDraw


def _make_atoms(x: float, symbol: str = "H") -> ase.Atoms:
    """Create a simple atom at position (x, 0, 0) with info tag."""
    atoms = ase.Atoms(symbol, positions=[[x, 0, 0]])
    atoms.info["x"] = x
    return atoms


def _set_default_room(server_url: str, token: str, room_id: str) -> None:
    """Set a room as the server's default template."""
    response = httpx.put(
        f"{server_url}/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
        json={"room_id": room_id},
    )
    response.raise_for_status()


def _unset_default_room(server_url: str, token: str) -> None:
    """Unset the server's default template."""
    response = httpx.delete(
        f"{server_url}/v1/server-settings/default-room",
        headers={"Authorization": f"Bearer {token}"},
    )
    response.raise_for_status()


@pytest.mark.parametrize("storage_type", ["memory", "lmdb"])
def test_template_room_isolation(
    server_factory: Callable[[dict[str, str]], Any],
    tmp_path: pathlib.Path,
    storage_type: str,
) -> None:
    """Template room frames are not affected by modifications to derived rooms.

    Scenario:
    1. Admin creates room A with 5 frames, locks it, sets as default template
    2. Guest creates room B (should copy A's 5 frames)
    3. Guest appends 3 frames to room B
    4. Verify:
       - Room A still has exactly 5 frames (unchanged)
       - Room B has exactly 8 frames (5 + 3)
       - Room B's first 5 frames match room A's frames
       - Room A's frames are byte-for-byte identical to before
    """
    # Configure server with auth and storage backend
    env_overrides = {
        "ZNDRAW_AUTH_DEFAULT_ADMIN_EMAIL": "admin@local.test",
        "ZNDRAW_AUTH_DEFAULT_ADMIN_PASSWORD": "adminpassword",
    }
    if storage_type == "lmdb":
        lmdb_path = tmp_path / "test_storage.lmdb"
        env_overrides["ZNDRAW_SERVER_STORAGE"] = str(lmdb_path)

    server_instance = server_factory(env_overrides)
    server_url = server_instance.url

    template_room_id = uuid.uuid4().hex
    template_frame_count = 5
    append_frame_count = 3

    # Step 1: Admin creates template room with 5 frames
    admin = ZnDraw(
        url=server_url,
        room=template_room_id,
        user="admin@local.test",
        password="adminpassword",
    )

    template_frames = [_make_atoms(float(i), "C") for i in range(template_frame_count)]
    admin.extend(template_frames)
    assert len(admin) == template_frame_count

    # Capture template frame data for later comparison
    template_positions = [
        admin[i].positions.copy() for i in range(template_frame_count)
    ]
    template_symbols = [
        admin[i].get_chemical_symbols() for i in range(template_frame_count)
    ]

    # Lock the room (admin lock)
    admin.api.update_room({"locked": True})

    # Set as default template
    token = admin.api.token
    assert token is not None
    _set_default_room(server_url, token, template_room_id)

    # Step 2: Guest creates a new room (should copy from template)
    guest = ZnDraw(url=server_url, copy_from=None)  # None = use server default template
    derived_room_id = guest.room

    assert derived_room_id != template_room_id, "Guest should get a new room"
    assert len(guest) == template_frame_count, (
        f"New room should have {template_frame_count} frames from template, "
        f"got {len(guest)}"
    )

    # Verify derived room has same frames as template
    for i in range(template_frame_count):
        np.testing.assert_allclose(
            guest[i].positions,
            template_positions[i],
            atol=1e-6,
            err_msg=f"Frame {i} positions don't match template",
        )
        assert guest[i].get_chemical_symbols() == template_symbols[i], (
            f"Frame {i} symbols don't match template"
        )

    # Step 3: Guest appends frames to derived room
    new_frames = [_make_atoms(100.0 + i, "O") for i in range(append_frame_count)]
    guest.extend(new_frames)

    expected_derived_length = template_frame_count + append_frame_count
    assert len(guest) == expected_derived_length, (
        f"Derived room should have {expected_derived_length} frames, got {len(guest)}"
    )

    # Step 4: Verify template room is unchanged
    # Re-fetch admin client to ensure fresh data
    admin_check = ZnDraw(
        url=server_url,
        room=template_room_id,
        user="admin@local.test",
        password="adminpassword",
    )

    assert len(admin_check) == template_frame_count, (
        f"Template room should still have {template_frame_count} frames, "
        f"got {len(admin_check)}"
    )

    # Verify template frames are identical
    for i in range(template_frame_count):
        np.testing.assert_allclose(
            admin_check[i].positions,
            template_positions[i],
            atol=1e-6,
            err_msg=f"Template frame {i} positions were modified!",
        )
        assert admin_check[i].get_chemical_symbols() == template_symbols[i], (
            f"Template frame {i} symbols were modified!"
        )

    # Verify derived room has all expected frames
    # First 5 should match template
    for i in range(template_frame_count):
        np.testing.assert_allclose(
            guest[i].positions,
            template_positions[i],
            atol=1e-6,
            err_msg=f"Derived room frame {i} doesn't match template",
        )

    # Last 3 should be the appended oxygen atoms
    for i in range(append_frame_count):
        frame_idx = template_frame_count + i
        expected_x = 100.0 + i
        np.testing.assert_allclose(
            guest[frame_idx].positions[0, 0],
            expected_x,
            atol=1e-6,
            err_msg=f"Appended frame {i} has wrong position",
        )
        assert guest[frame_idx].get_chemical_symbols() == ["O"], (
            f"Appended frame {i} should be oxygen"
        )

    # Cleanup
    _unset_default_room(server_url, token)
    admin.disconnect()
    admin_check.disconnect()
    guest.disconnect()
