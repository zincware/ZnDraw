"""Layer 2: runtime broadcast contract sweep across emit families.

For each representative action against a real room, the captured event must:
  * be a RoomScopedEvent subclass
  * carry ``room_id == UUID(room.id)`` and ``room_address == room.public_address``
  * be emitted on channel ``f"room:{room.id}"``

Sigil channels (``@global``, ``@internal``) are whitelisted: this test only
drives real-room actions, so no sigil channel should appear.
"""

from __future__ import annotations

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
from uuid import UUID

import zndraw.socket_events  # noqa: F401 — registers all subclasses
import zndraw_joblib.events  # noqa: F401 — registers joblib subclasses

from zndraw.client import atoms_to_json_dict
from zndraw.socket_events import RoomScopedEvent


def _make_json_frame(formula: str = "H2") -> dict:
    atoms = ase.Atoms(
        formula,
        positions=[[i, 0, 0] for i in range(ase.Atoms(formula).get_global_number_of_atoms())],
    )
    return atoms_to_json_dict(atoms)


def _all_subclasses(cls: type) -> set[type]:
    out: set[type] = set()
    stack: list[type] = list(cls.__subclasses__())
    while stack:
        sub = stack.pop()
        if sub in out:
            continue
        out.add(sub)
        stack.extend(sub.__subclasses__())
    return out


def _room_scoped_event_names() -> set[str]:
    names: set[str] = set()
    for cls in _all_subclasses(RoomScopedEvent):
        snake = "".join(
            f"_{c.lower()}" if c.isupper() else c for c in cls.__name__
        ).lstrip("_")
        names.add(snake)
    return names


@pytest.mark.protected
@pytest.mark.asyncio
async def test_every_route_action_emits_consistent_room_scoped_events(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    base = f"/v1/rooms/{room.public_address}"
    headers = auth_header(token)

    # frames family
    r = await client.post(
        f"{base}/frames", json={"frames": [_make_json_frame()]}, headers=headers
    )
    assert r.status_code == 201, r.text

    # bookmarks family
    r = await client.put(
        f"{base}/bookmarks/0", json={"label": "x"}, headers=headers
    )
    assert r.status_code == 200, r.text

    # figures family — data must be a JSON string per FigureData schema
    r = await client.post(
        f"{base}/figures/fig",
        json={"figure": {"type": "plotly", "data": '{"data":[]}'}},
        headers=headers,
    )
    assert r.status_code == 201, r.text

    # chat family
    r = await client.post(
        f"{base}/chat/messages", json={"content": "hi"}, headers=headers
    )
    assert r.status_code == 201, r.text

    # progress family
    r = await client.post(
        f"{base}/progress",
        json={"progress_id": "p", "description": "d", "unit": "it"},
        headers=headers,
    )
    assert r.status_code == 201, r.text

    # joblib family
    r = await client.put(
        f"/v1/joblib/rooms/{room.public_address}/jobs",
        json={
            "category": "analysis",
            "name": "noop",
            "schema_": {"type": "object"},
        },
        headers=headers,
    )
    assert r.status_code in (200, 201), r.text

    rs_event_names = _room_scoped_event_names()
    captured = [e for e in mock_sio.emitted if e["event"] in rs_event_names]
    assert captured, "no RoomScopedEvent emissions observed across the matrix"

    room_uuid = UUID(room.id)
    for emit in captured:
        data = emit["data"]
        observed = data["room_id"]
        observed_uuid = observed if isinstance(observed, UUID) else UUID(str(observed))
        assert observed_uuid == room_uuid, (
            f"{emit['event']}: room_id={observed!r} != {room.id!r}"
        )
        assert data["room_address"] == room.public_address, (
            f"{emit['event']}: room_address={data['room_address']!r} != "
            f"{room.public_address!r}"
        )
        assert emit["room"] == f"room:{room.id}", (
            f"{emit['event']}: channel={emit['room']!r} != room:{room.id}"
        )

    # Sigil whitelist: no sigil-channel emit should occur during a real-room sweep.
    for emit in mock_sio.emitted:
        ch = emit.get("room")
        if isinstance(ch, str) and ch.startswith("room:@"):
            pytest.fail(
                f"unexpected sigil channel emit during real-room sweep: "
                f"event={emit['event']} channel={ch}"
            )
