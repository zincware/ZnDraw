"""Regression tests: non-superuser stranger must not read content from a private room.

Covers the existence-leak + unauthorized-read vulnerability in frames /
trajectory / screenshots / isosurface / progress / frame-selection routes.
"""

import pytest
from httpx import AsyncClient
from sqlmodel.ext.asyncio.session import AsyncSession

from helpers import _register_and_login, create_test_user_in_db


@pytest.mark.asyncio
async def test_stranger_cannot_read_frames_from_private_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """A stranger must receive 404 (not 200 or 403) for all content GETs on a private room.

    Parameters
    ----------
    client
        Async test client with real Redis and DB.
    session
        Async database session (shared with client fixture).
    """
    owner = await _register_and_login(client, "cg-own@test.com")
    _s, stranger = await create_test_user_in_db(
        session, email="cg-str@test.com", is_superuser=False
    )
    r = await client.post(
        "/v1/rooms",
        json={"room_id": "cg-priv", "visibility": "private"},
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code in (200, 201), f"Room creation failed: {r.status_code} {r.text}"

    # Each of these should 404 (existence hidden) for the stranger.
    # Note: isosurface is at /frames/{index}/isosurface — tested via frames path.
    # Note: progress only has POST/PATCH/DELETE — no GET list route exists.
    for path in [
        "/v1/rooms/cg-priv/frames",
        "/v1/rooms/cg-priv/frames/0",
        "/v1/rooms/cg-priv/frames/0/metadata",
        "/v1/rooms/cg-priv/frames/0/isosurface?cube_key=k&isovalue=0.5",
        "/v1/rooms/cg-priv/trajectory",
        "/v1/rooms/cg-priv/screenshots",
        "/v1/rooms/cg-priv/frame-selection",
        "/v1/rooms/cg-priv/edit-lock",
    ]:
        r = await client.get(
            path, headers={"Authorization": f"Bearer {stranger}"}
        )
        assert r.status_code == 404, f"{path} returned {r.status_code}: {r.text}"


@pytest.mark.asyncio
async def test_owner_can_read_own_private_room_content(
    client: AsyncClient,
) -> None:
    """Owner must be able to read content from their own private room.

    Parameters
    ----------
    client
        Async test client with real Redis and DB.
    """
    owner = await _register_and_login(client, "cg-own-ok@test.com")
    r = await client.post(
        "/v1/rooms",
        json={"room_id": "cg-own-priv", "visibility": "private"},
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code in (200, 201), f"Room creation failed: {r.status_code} {r.text}"

    # Owner can read at least /frames (sanity check that the positive path works)
    r = await client.get(
        "/v1/rooms/cg-own-priv/frames",
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code == 200, f"Owner got {r.status_code}: {r.text}"


@pytest.mark.asyncio
async def test_stranger_cannot_write_to_private_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """A stranger must receive 404 for write operations on a private room.

    Parameters
    ----------
    client
        Async test client with real Redis and DB.
    session
        Async database session.
    """
    owner = await _register_and_login(client, "cg-write-own@test.com")
    _s, stranger = await create_test_user_in_db(
        session, email="cg-write-str@test.com", is_superuser=False
    )
    r = await client.post(
        "/v1/rooms",
        json={"room_id": "cg-write-priv", "visibility": "private"},
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code in (200, 201)

    # POST progress — stranger should see 404
    r = await client.post(
        "/v1/rooms/cg-write-priv/progress",
        json={"progress_id": "p1", "description": "test"},
        headers={"Authorization": f"Bearer {stranger}"},
    )
    assert r.status_code == 404, f"POST /progress returned {r.status_code}: {r.text}"

    # POST /trajectory/download-tokens — stranger should see 404
    r = await client.post(
        "/v1/rooms/cg-write-priv/trajectory/download-tokens",
        headers={"Authorization": f"Bearer {stranger}"},
    )
    assert r.status_code == 404, f"POST /download-tokens returned {r.status_code}: {r.text}"
