"""Scope-filtered room listing tests."""

import pytest
from helpers import (
    _register_and_login,
    create_room_via_api,
    create_test_user_in_db,
    get_user_id,
)
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_list_rooms_union(client: AsyncClient) -> None:
    token_a = await _register_and_login(client, "list-a@test.com")
    token_b = await _register_and_login(client, "list-b@test.com")
    owner_id_a = await get_user_id(client, token_a)

    # User A creates: public + private rooms
    pub_id = await create_room_via_api(client, token_a, "r-pub", visibility="public")
    priv_id = await create_room_via_api(client, token_a, "r-priv", visibility="private")

    # User B creates a group; A joins as a viewer; B creates a group-room.
    group_resp = await client.post(
        "/v1/groups",
        json={"name": "listing-test"},
        headers={"Authorization": f"Bearer {token_b}"},
    )
    group_body = group_resp.json()
    gid = group_body["id"]
    group_name = group_body["name"]
    add_r = await client.post(
        f"/v1/groups/{gid}/members",
        json={"user_id": owner_id_a, "role": "viewer"},
        headers={"Authorization": f"Bearer {token_b}"},
    )
    assert add_r.status_code == 201
    grp_id = await create_room_via_api(
        client, token_b, "r-grp", owner=group_name, visibility="group"
    )

    # A sees: public (r-pub), own private (r-priv), group member (r-grp)
    r = await client.get("/v1/rooms", headers={"Authorization": f"Bearer {token_a}"})
    ids = {it["room_id"] for it in r.json()["items"]}
    assert {pub_id, priv_id, grp_id} <= ids

    # B sees: public (r-pub), own group (r-grp); not A's private
    r = await client.get("/v1/rooms", headers={"Authorization": f"Bearer {token_b}"})
    ids = {it["room_id"] for it in r.json()["items"]}
    assert {pub_id, grp_id} <= ids
    assert priv_id not in ids


@pytest.mark.asyncio
async def test_list_rooms_requires_auth(client: AsyncClient) -> None:
    r = await client.get("/v1/rooms")
    assert r.status_code == 401


@pytest.mark.asyncio
async def test_create_private_room_sets_owner(client: AsyncClient) -> None:
    token = await _register_and_login(client, "priv-own@test.com")
    room_id = await create_room_via_api(
        client, token, "priv-own-1", visibility="private"
    )
    r = await client.get(
        f"/v1/rooms/{room_id}",
        headers={"Authorization": f"Bearer {token}"},
    )
    body = r.json()
    assert body["visibility"] == "private"
    assert body["owner"] is not None
    assert body["owner_kind"] == "user"


@pytest.mark.asyncio
async def test_create_group_room_requires_membership(
    client: AsyncClient, session
) -> None:
    """Non-superuser outsider cannot create a group-room in a group they don't join."""
    _outsider_user, outsider = await create_test_user_in_db(
        session, email="goro-out@test.com", is_superuser=False
    )
    owner = await _register_and_login(client, "goro@test.com")
    group_body = (
        await client.post(
            "/v1/groups",
            json={"name": "gr-owned"},
            headers={"Authorization": f"Bearer {owner}"},
        )
    ).json()

    r = await client.post(
        "/v1/rooms",
        json={
            "owner": group_body["name"],
            "name": "gr-owned-room",
            "visibility": "group",
        },
        headers={"Authorization": f"Bearer {outsider}"},
    )
    assert r.status_code in (403, 404, 409)


@pytest.mark.asyncio
async def test_patch_room_visibility(client: AsyncClient) -> None:
    token = await _register_and_login(client, "patch-vis@test.com")
    room_id = await create_room_via_api(
        client, token, "patch-vis-r", visibility="public"
    )
    r = await client.patch(
        f"/v1/rooms/{room_id}",
        json={"visibility": "private"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert r.status_code == 200
    r = await client.get(
        f"/v1/rooms/{room_id}",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert r.json()["visibility"] == "private"


@pytest.mark.asyncio
async def test_get_private_room_404_to_non_owner(client: AsyncClient, session) -> None:
    _stranger_user, stranger = await create_test_user_in_db(
        session, email="pgs@test.com", is_superuser=False
    )
    owner = await _register_and_login(client, "pgo@test.com")
    room_id = await create_room_via_api(client, owner, "pg-room", visibility="private")
    r = await client.get(
        f"/v1/rooms/{room_id}",
        headers={"Authorization": f"Bearer {stranger}"},
    )
    assert r.status_code == 404
    assert r.json()["type"].endswith("/room-not-found")
