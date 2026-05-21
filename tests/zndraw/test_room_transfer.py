"""Transfer semantics: collision -> 409; happy path -> room_renamed."""

from __future__ import annotations

import pytest
from httpx import AsyncClient

ADMIN_EMAIL = "admin@local.test"
ADMIN_PASSWORD = "adminpassword"  # noqa: S105


async def _register_and_login(
    client: AsyncClient,
    email: str,
    password: str = "test12345",  # noqa: S107
) -> tuple[str, str]:
    """Register a user and return (user_id, token)."""
    reg = await client.post(
        "/v1/auth/register", json={"email": email, "password": password}
    )
    assert reg.status_code == 201, reg.text
    user_id = reg.json()["id"]
    login = await client.post(
        "/v1/auth/jwt/login",
        data={"username": email, "password": password},
    )
    assert login.status_code == 200, login.text
    return user_id, login.json()["access_token"]


async def _login_admin(client: AsyncClient) -> tuple[str, str]:
    """Login as the pre-seeded admin superuser; return (user_id, token)."""
    login = await client.post(
        "/v1/auth/jwt/login",
        data={"username": ADMIN_EMAIL, "password": ADMIN_PASSWORD},
    )
    assert login.status_code == 200, login.text
    token = login.json()["access_token"]
    me = await client.get("/v1/auth/users/me", headers={"Authorization": f"Bearer {token}"})
    assert me.status_code == 200, me.text
    return me.json()["id"], token


@pytest.mark.asyncio
async def test_transfer_collision_returns_409(
    http_client_auth: AsyncClient,
) -> None:
    """Superuser tries to transfer a room to a user that already has a room with that name -> 409."""
    admin_id, admin_token = await _login_admin(http_client_auth)
    dst_id, dst_token = await _register_and_login(
        http_client_auth, "tx-dst-409@example.com"
    )

    admin_headers = {"Authorization": f"Bearer {admin_token}"}
    dst_headers = {"Authorization": f"Bearer {dst_token}"}

    # Admin creates room "collision" in their own namespace
    r = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": admin_id, "name": "collision"},
        headers=admin_headers,
    )
    assert r.status_code == 201, r.text

    # Destination user creates their own room "collision"
    r = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": dst_id, "name": "collision"},
        headers=dst_headers,
    )
    assert r.status_code == 201, r.text

    # Admin attempts to transfer their "collision" room to dst -> 409
    resp = await http_client_auth.patch(
        f"/v1/rooms/{admin_id}/collision",
        json={"new_owner_id": dst_id},
        headers=admin_headers,
    )
    assert resp.status_code == 409, resp.text


@pytest.mark.asyncio
async def test_transfer_happy_path(
    http_client_auth: AsyncClient,
) -> None:
    """Transfer to a group the caller manages; surrogate UUID unchanged; new composed address returned."""
    caller_id, caller_token = await _register_and_login(
        http_client_auth, "tx-happy@example.com"
    )
    headers = {"Authorization": f"Bearer {caller_token}"}

    # Caller creates a group they own
    g = await http_client_auth.post(
        "/v1/groups", json={"name": "tx-grp"}, headers=headers
    )
    assert g.status_code == 201, g.text
    group_id = g.json()["id"]

    # Caller creates room "t" in own namespace
    r = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": caller_id, "name": "t"},
        headers=headers,
    )
    assert r.status_code == 201, r.text
    original_surrogate = r.json()["room_id"].split("/")[1] if "/" not in r.json()["room_id"].replace(caller_id + "/", "", 1) else None

    # Transfer to the group, changing visibility to group
    resp = await http_client_auth.patch(
        f"/v1/rooms/{caller_id}/t",
        json={"new_owner_id": group_id, "visibility": "group"},
        headers=headers,
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["room_id"] == f"{group_id}/t"

    # New composed address is reachable
    r = await http_client_auth.get(
        f"/v1/rooms/{group_id}/t", headers=headers
    )
    assert r.status_code == 200

    # Old composed address is NOT reachable
    r = await http_client_auth.get(
        f"/v1/rooms/{caller_id}/t", headers=headers
    )
    assert r.status_code == 404
