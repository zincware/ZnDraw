"""Regression tests: non-superuser stranger must not read content from a private room.

Covers the existence-leak + unauthorized-read vulnerability in frames /
trajectory / screenshots / isosurface / progress / frame-selection routes.
"""

import pytest
from helpers import _register_and_login, create_room_via_api, create_test_user_in_db
from httpx import AsyncClient
from sqlmodel.ext.asyncio.session import AsyncSession


@pytest.mark.asyncio
async def test_stranger_cannot_read_frames_from_private_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Stranger gets 404 (not 200/403) for all content GETs on a private room."""
    owner = await _register_and_login(client, "cg-own@test.com")
    _s, stranger = await create_test_user_in_db(
        session, email="cg-str@test.com", is_superuser=False
    )
    room_id = await create_room_via_api(client, owner, "cg-priv", visibility="private")

    for path in [
        f"/v1/rooms/{room_id}/frames",
        f"/v1/rooms/{room_id}/frames/0",
        f"/v1/rooms/{room_id}/frames/0/metadata",
        f"/v1/rooms/{room_id}/frames/0/isosurface?cube_key=k&isovalue=0.5",
        f"/v1/rooms/{room_id}/trajectory",
        f"/v1/rooms/{room_id}/screenshots",
        f"/v1/rooms/{room_id}/frame-selection",
        f"/v1/rooms/{room_id}/edit-lock",
    ]:
        r = await client.get(path, headers={"Authorization": f"Bearer {stranger}"})
        assert r.status_code == 404, f"{path} returned {r.status_code}: {r.text}"


@pytest.mark.asyncio
async def test_owner_can_read_own_private_room_content(
    client: AsyncClient,
) -> None:
    """Owner must be able to read content from their own private room."""
    owner = await _register_and_login(client, "cg-own-ok@test.com")
    room_id = await create_room_via_api(
        client, owner, "cg-own-priv", visibility="private"
    )

    # Owner can read at least /frames (sanity check that the positive path works)
    r = await client.get(
        f"/v1/rooms/{room_id}/frames",
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code == 200, f"Owner got {r.status_code}: {r.text}"


@pytest.mark.asyncio
async def test_stranger_cannot_write_to_private_room(
    client: AsyncClient, session: AsyncSession
) -> None:
    """A stranger must receive 404 for write operations on a private room."""
    owner = await _register_and_login(client, "cg-write-own@test.com")
    _s, stranger = await create_test_user_in_db(
        session, email="cg-write-str@test.com", is_superuser=False
    )
    room_id = await create_room_via_api(
        client, owner, "cg-write-priv", visibility="private"
    )

    # POST progress — stranger should see 404
    r = await client.post(
        f"/v1/rooms/{room_id}/progress",
        json={"progress_id": "p1", "description": "test"},
        headers={"Authorization": f"Bearer {stranger}"},
    )
    assert r.status_code == 404, f"POST /progress returned {r.status_code}: {r.text}"

    # POST /trajectory/download-tokens — stranger should see 404
    r = await client.post(
        f"/v1/rooms/{room_id}/trajectory/download-tokens",
        headers={"Authorization": f"Bearer {stranger}"},
    )
    assert r.status_code == 404, (
        f"POST /download-tokens returned {r.status_code}: {r.text}"
    )
