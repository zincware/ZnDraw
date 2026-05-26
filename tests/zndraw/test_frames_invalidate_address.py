"""Regression: FramesInvalidate carries composed room_address (finding #2).

Before the wire-convention refactor, the route emitted ``room_id=room.id``
(surrogate) but the frontend keyed React Query by the composed address, so
the predicate never matched. This test pins ``room_address`` to the composed
form.
"""

import ase
import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.client import atoms_to_json_dict


def _make_json_frame(formula: str = "H2") -> dict:
    atoms = ase.Atoms(
        formula,
        positions=[
            [i, 0, 0] for i in range(ase.Atoms(formula).get_global_number_of_atoms())
        ],
    )
    return atoms_to_json_dict(atoms)


@pytest.mark.asyncio
async def test_append_frame_emits_frames_invalidate_with_composed_address(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.post(
        f"/v1/rooms/{room.public_address}/frames",
        json={"frames": [_make_json_frame("H2")]},
        headers=auth_header(token),
    )
    assert response.status_code == 201, response.text

    invalidates = [e for e in mock_sio.emitted if e["event"] == "frames_invalidate"]
    assert len(invalidates) == 1, (
        f"expected 1 frames_invalidate, got {len(invalidates)}"
    )
    captured = invalidates[0]
    data = captured["data"]
    assert "room_id" in data, (
        f"FramesInvalidate payload missing room_id; got keys: {list(data)}"
    )
    assert "room_address" in data, (
        f"FramesInvalidate payload missing room_address; got keys: {list(data)}"
    )
    assert str(data["room_id"]) == room.id
    assert data["room_address"] == room.public_address
    assert captured["room"] == f"room:{room.id}"
