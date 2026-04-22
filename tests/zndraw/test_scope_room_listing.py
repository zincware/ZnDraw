"""Scope-filtered room listing tests."""

import pytest
from helpers import _register_and_login
from httpx import AsyncClient


@pytest.mark.asyncio
async def test_list_rooms_union(client: AsyncClient) -> None:
    token_a = await _register_and_login(client, "list-a@test.com")
    token_b = await _register_and_login(client, "list-b@test.com")

    # User A creates: public + private rooms
    r = await client.post(
        "/v1/rooms",
        json={"room_id": "r-pub", "visibility": "public"},
        headers={"Authorization": f"Bearer {token_a}"},
    )
    assert r.status_code == 201
    r = await client.post(
        "/v1/rooms",
        json={"room_id": "r-priv", "visibility": "private"},
        headers={"Authorization": f"Bearer {token_a}"},
    )
    assert r.status_code == 201

    # User B creates a group; A joins as a viewer; B creates a group-room.
    gid = (
        await client.post(
            "/v1/groups",
            json={"name": "listing-test"},
            headers={"Authorization": f"Bearer {token_b}"},
        )
    ).json()["id"]
    me_a = (
        await client.get(
            "/v1/auth/users/me", headers={"Authorization": f"Bearer {token_a}"}
        )
    ).json()
    add_r = await client.post(
        f"/v1/groups/{gid}/members",
        json={"user_id": me_a["id"], "role": "viewer"},
        headers={"Authorization": f"Bearer {token_b}"},
    )
    assert add_r.status_code == 201
    r = await client.post(
        "/v1/rooms",
        json={"room_id": "r-grp", "visibility": "group", "owner_group_id": gid},
        headers={"Authorization": f"Bearer {token_b}"},
    )
    assert r.status_code == 201

    # A sees: public (r-pub), own private (r-priv), group member (r-grp)
    r = await client.get("/v1/rooms", headers={"Authorization": f"Bearer {token_a}"})
    ids = {it["id"] for it in r.json()["items"]}
    assert {"r-pub", "r-priv", "r-grp"} <= ids

    # B sees: public (r-pub), own group (r-grp); not A's private
    r = await client.get("/v1/rooms", headers={"Authorization": f"Bearer {token_b}"})
    ids = {it["id"] for it in r.json()["items"]}
    assert {"r-pub", "r-grp"} <= ids
    assert "r-priv" not in ids


@pytest.mark.asyncio
async def test_list_rooms_requires_auth(client: AsyncClient) -> None:
    r = await client.get("/v1/rooms")
    assert r.status_code == 401


@pytest.mark.asyncio
async def test_create_private_room_sets_owner(client: AsyncClient) -> None:
    token = await _register_and_login(client, "priv-own@test.com")
    r = await client.post(
        "/v1/rooms",
        json={"room_id": "priv-own-1", "visibility": "private"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert r.status_code == 201
    r = await client.get(
        "/v1/rooms/priv-own-1",
        headers={"Authorization": f"Bearer {token}"},
    )
    body = r.json()
    assert body["visibility"] == "private"
    assert body["owner_user_id"] is not None
    assert body["owner_group_id"] is None


@pytest.mark.asyncio
async def test_create_group_room_requires_membership(
    client: AsyncClient, session
) -> None:
    """Non-superuser outsider cannot create a group-room in a group
    they don't belong to."""
    from helpers import create_test_user_in_db

    owner = await _register_and_login(client, "goro@test.com")
    # Bypass register flow (dev-mode auto-promotes to superuser)
    _outsider_user, outsider = await create_test_user_in_db(
        session, email="goro-out@test.com", is_superuser=False
    )
    gid = (
        await client.post(
            "/v1/groups",
            json={"name": "gr-owned"},
            headers={"Authorization": f"Bearer {owner}"},
        )
    ).json()["id"]

    r = await client.post(
        "/v1/rooms",
        json={
            "room_id": "gr-owned-room",
            "visibility": "group",
            "owner_group_id": gid,
        },
        headers={"Authorization": f"Bearer {outsider}"},
    )
    assert r.status_code == 409
    assert r.json()["type"].endswith("/transfer-target-invalid")


@pytest.mark.asyncio
async def test_patch_room_visibility(client: AsyncClient) -> None:
    token = await _register_and_login(client, "patch-vis@test.com")
    await client.post(
        "/v1/rooms",
        json={"room_id": "patch-vis-r", "visibility": "public"},
        headers={"Authorization": f"Bearer {token}"},
    )
    r = await client.patch(
        "/v1/rooms/patch-vis-r",
        json={"visibility": "private"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert r.status_code == 200
    r = await client.get(
        "/v1/rooms/patch-vis-r",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert r.json()["visibility"] == "private"


@pytest.mark.asyncio
async def test_get_private_room_404_to_non_owner(client: AsyncClient, session) -> None:
    from helpers import create_test_user_in_db

    owner = await _register_and_login(client, "pgo@test.com")
    # Bypass register flow (dev-mode auto-promotes to superuser)
    _stranger_user, stranger = await create_test_user_in_db(
        session, email="pgs@test.com", is_superuser=False
    )
    await client.post(
        "/v1/rooms",
        json={"room_id": "pg-room", "visibility": "private"},
        headers={"Authorization": f"Bearer {owner}"},
    )
    r = await client.get(
        "/v1/rooms/pg-room",
        headers={"Authorization": f"Bearer {stranger}"},
    )
    assert r.status_code == 404
    assert r.json()["type"].endswith("/room-not-found")
