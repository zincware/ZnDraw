"""Share-link REST + permission-compose tests."""

import pytest
from helpers import (
    _register_and_login,
    create_room_via_api,
    create_test_user_in_db,
)
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_create_and_use_view_link(client: AsyncClient, session) -> None:
    owner = await _register_and_login(client, "sl1-own@test.com")
    _guest_user, guest = await create_test_user_in_db(
        session, email="sl1-gst@test.com", is_superuser=False
    )
    room_id = await create_room_via_api(
        client, owner, "sl-room-1", visibility="private"
    )
    # 404 without share token
    r = await client.get(
        f"/v1/rooms/{room_id}", headers={"Authorization": f"Bearer {guest}"}
    )
    assert r.status_code == 404

    # Create view link
    r = await client.post(
        f"/v1/rooms/{room_id}/share-links",
        json={"access": "view"},
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code == 201
    token = r.json()["token"]

    # Guest can read with token
    r = await client.get(
        f"/v1/rooms/{room_id}",
        headers={"Authorization": f"Bearer {guest}", "X-Room-Share-Token": token},
    )
    assert r.status_code == 200


@pytest.mark.asyncio
async def test_list_share_links_requires_manage(client: AsyncClient, session) -> None:
    owner = await _register_and_login(client, "sl-list-own@test.com")
    _other_user, other = await create_test_user_in_db(
        session, email="sl-list-other@test.com", is_superuser=False
    )
    room_id = await create_room_via_api(client, owner, "sl-list-1", visibility="public")
    # Owner can list (empty so far)
    r = await client.get(
        f"/v1/rooms/{room_id}/share-links",
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code == 200
    # Non-owner on public room cannot list (403 — can_read passes, can_manage doesn't)
    r = await client.get(
        f"/v1/rooms/{room_id}/share-links",
        headers={"Authorization": f"Bearer {other}"},
    )
    assert r.status_code == 403


@pytest.mark.asyncio
async def test_revoke_link(client: AsyncClient, session) -> None:
    owner = await _register_and_login(client, "sl2-own@test.com")
    _guest_user, guest = await create_test_user_in_db(
        session, email="sl2-gst@test.com", is_superuser=False
    )
    room_id = await create_room_via_api(
        client, owner, "sl-room-2", visibility="private"
    )
    link = (
        await client.post(
            f"/v1/rooms/{room_id}/share-links",
            json={"access": "view"},
            headers={"Authorization": f"Bearer {owner}"},
        )
    ).json()
    r = await client.delete(
        f"/v1/rooms/{room_id}/share-links/{link['id']}",
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code == 204

    # Guest can no longer use the revoked token
    r = await client.get(
        f"/v1/rooms/{room_id}",
        headers={
            "Authorization": f"Bearer {guest}",
            "X-Room-Share-Token": link["token"],
        },
    )
    assert r.status_code == 404


@pytest.mark.asyncio
async def test_wrong_room_token_rejected(client: AsyncClient, session) -> None:
    owner = await _register_and_login(client, "sl4-own@test.com")
    _guest_user, guest = await create_test_user_in_db(
        session, email="sl4-gst@test.com", is_superuser=False
    )
    room_id_a = await create_room_via_api(
        client, owner, "sl-room-4a", visibility="private"
    )
    room_id_b = await create_room_via_api(
        client, owner, "sl-room-4b", visibility="private"
    )
    link = (
        await client.post(
            f"/v1/rooms/{room_id_a}/share-links",
            json={"access": "view"},
            headers={"Authorization": f"Bearer {owner}"},
        )
    ).json()
    r = await client.get(
        f"/v1/rooms/{room_id_b}",
        headers={
            "Authorization": f"Bearer {guest}",
            "X-Room-Share-Token": link["token"],
        },
    )
    assert r.status_code == 404


@pytest.mark.asyncio
async def test_revoke_nonexistent_returns_404(
    client: AsyncClient,
) -> None:
    from uuid import uuid4

    owner = await _register_and_login(client, "sl-rn-own@test.com")
    room_id = await create_room_via_api(
        client, owner, "sl-rn-room", visibility="public"
    )
    r = await client.delete(
        f"/v1/rooms/{room_id}/share-links/{uuid4()}",
        headers={"Authorization": f"Bearer {owner}"},
    )
    assert r.status_code == 404
    assert r.json()["type"].endswith("/share-link-not-found")
