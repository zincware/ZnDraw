"""Integration tests for /v1/groups endpoints."""

import pytest
from helpers import _register_and_login, auth_header, create_test_user_in_db
from httpx import AsyncClient
from sqlmodel.ext.asyncio.session import AsyncSession


@pytest.mark.asyncio
async def test_create_group_creator_is_admin(client: AsyncClient) -> None:
    token = await _register_and_login(client, "gc1@test.com")
    r = await client.post(
        "/v1/groups",
        json={"name": "alpha"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert r.status_code == 201
    gid = r.json()["id"]

    r = await client.get(
        f"/v1/groups/{gid}", headers={"Authorization": f"Bearer {token}"}
    )
    assert r.status_code == 200
    assert r.json()["my_role"] == "admin"


@pytest.mark.asyncio
async def test_group_name_unique(client: AsyncClient) -> None:
    t1 = await _register_and_login(client, "gc2a@test.com")
    t2 = await _register_and_login(client, "gc2b@test.com")
    r = await client.post(
        "/v1/groups", json={"name": "same"}, headers={"Authorization": f"Bearer {t1}"}
    )
    assert r.status_code == 201
    r = await client.post(
        "/v1/groups", json={"name": "same"}, headers={"Authorization": f"Bearer {t2}"}
    )
    assert r.status_code == 409
    assert r.json()["type"].endswith("/group-name-taken")


@pytest.mark.asyncio
async def test_list_my_groups(client: AsyncClient) -> None:
    t = await _register_and_login(client, "gc3@test.com")
    for n in ["a", "b"]:
        await client.post(
            "/v1/groups",
            json={"name": f"list-{n}"},
            headers={"Authorization": f"Bearer {t}"},
        )
    r = await client.get("/v1/groups", headers={"Authorization": f"Bearer {t}"})
    assert r.status_code == 200
    names = {g["name"] for g in r.json()["items"]}
    assert {"list-a", "list-b"} <= names


@pytest.mark.asyncio
async def test_add_member_admin_only(
    client: AsyncClient, session: AsyncSession
) -> None:
    admin = await _register_and_login(client, "gadm@test.com")
    # Insert joiner directly with is_superuser=False to bypass dev-mode promotion
    joiner_user, joiner_token = await create_test_user_in_db(
        session, "gjoin@test.com", is_superuser=False
    )
    gid = (
        await client.post(
            "/v1/groups",
            json={"name": "add-test"},
            headers={"Authorization": f"Bearer {admin}"},
        )
    ).json()["id"]

    # Non-admin attempt → 403
    r = await client.post(
        f"/v1/groups/{gid}/members",
        json={"user_id": str(joiner_user.id), "role": "member"},
        headers=auth_header(joiner_token),
    )
    assert r.status_code == 403

    # Admin adds → 201
    r = await client.post(
        f"/v1/groups/{gid}/members",
        json={"user_id": str(joiner_user.id), "role": "member"},
        headers={"Authorization": f"Bearer {admin}"},
    )
    assert r.status_code == 201


@pytest.mark.asyncio
async def test_last_admin_cannot_self_demote(client: AsyncClient) -> None:
    admin = await _register_and_login(client, "la1@test.com")
    me = (
        await client.get(
            "/v1/auth/users/me", headers={"Authorization": f"Bearer {admin}"}
        )
    ).json()
    gid = (
        await client.post(
            "/v1/groups",
            json={"name": "last-admin"},
            headers={"Authorization": f"Bearer {admin}"},
        )
    ).json()["id"]

    # Try to demote yourself → 409 LastGroupAdmin
    r = await client.patch(
        f"/v1/groups/{gid}/members/{me['id']}",
        json={"role": "member"},
        headers={"Authorization": f"Bearer {admin}"},
    )
    assert r.status_code == 409
    assert r.json()["type"].endswith("/last-group-admin")


@pytest.mark.asyncio
async def test_delete_group_with_rooms_blocked(
    client: AsyncClient, session: AsyncSession
) -> None:
    """Reassign-before-delete invariant — insert Room directly via session
    (the bridge routes/rooms.py doesn't yet honor owner_group_id requests)."""
    from zndraw.access import Visibility
    from zndraw.models import Room

    admin = await _register_and_login(client, "dg1@test.com")
    gid = (
        await client.post(
            "/v1/groups",
            json={"name": "has-rooms"},
            headers={"Authorization": f"Bearer {admin}"},
        )
    ).json()["id"]

    # Insert a group-owned room directly via the session fixture
    from uuid import UUID

    session.add(
        Room(
            id="g-room",
            owner_group_id=UUID(gid),
            visibility=Visibility.GROUP,
        )
    )
    await session.commit()

    r = await client.delete(
        f"/v1/groups/{gid}", headers={"Authorization": f"Bearer {admin}"}
    )
    assert r.status_code == 409
    assert r.json()["type"].endswith("/group-has-rooms")


@pytest.mark.asyncio
async def test_default_role_on_add_is_viewer(client: AsyncClient) -> None:
    admin = await _register_and_login(client, "dv1@test.com")
    joiner = await _register_and_login(client, "dv2@test.com")
    gid = (
        await client.post(
            "/v1/groups",
            json={"name": "def-role"},
            headers={"Authorization": f"Bearer {admin}"},
        )
    ).json()["id"]
    me = (
        await client.get(
            "/v1/auth/users/me", headers={"Authorization": f"Bearer {joiner}"}
        )
    ).json()
    # Omit role → default VIEWER
    await client.post(
        f"/v1/groups/{gid}/members",
        json={"user_id": me["id"]},
        headers={"Authorization": f"Bearer {admin}"},
    )
    r = await client.get(
        f"/v1/groups/{gid}", headers={"Authorization": f"Bearer {joiner}"}
    )
    assert r.json()["my_role"] == "viewer"
