"""Regression: joblib emits to the surrogate channel (finding #1).

Before the wire-convention refactor, ``register_job`` emitted to
``room:<owner_uuid>/<name>`` (composed) while in-room sockets joined
``room:<surrogate>``. This test pins the surrogate channel.
"""

from __future__ import annotations

import pytest
from helpers import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession
from zndraw_joblib.settings import JobLibSettings


@pytest.mark.asyncio
async def test_register_job_emits_to_surrogate_channel(
    client: AsyncClient,
    session: AsyncSession,
    mock_sio: MockSioServer,
) -> None:
    from zndraw.app import app

    app.state.joblib_settings = JobLibSettings()
    app.state.tsio = mock_sio

    user, token = await create_test_user_in_db(session)
    room = await create_test_room(session, user)

    response = await client.put(
        f"/v1/joblib/rooms/{room.public_address}/jobs",
        json={
            "category": "analysis",
            "name": "noop",
            "schema_": {"type": "object"},
        },
        headers=auth_header(token),
    )
    assert response.status_code in (200, 201), response.text

    invalidate_emits = [
        e for e in mock_sio.emitted if e["event"] == "jobs_invalidate"
    ]
    assert invalidate_emits, "register_job did not emit jobs_invalidate"
    captured = invalidate_emits[-1]
    assert captured["room"] == f"room:{room.id}", (
        f"expected channel room:{room.id}, got {captured['room']}"
    )
    # MockSioServer captures via model_dump(); UUID fields surface as UUID instances.
    assert str(captured["data"]["room_id"]) == room.id
    assert captured["data"]["room_address"] == room.public_address
