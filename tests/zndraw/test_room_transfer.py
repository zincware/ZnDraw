"""Transfer semantics: collision -> 409; happy path -> room_renamed."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

if TYPE_CHECKING:
    from httpx import AsyncClient

ADMIN_EMAIL = "admin@local.test"
ADMIN_PASSWORD = "adminpassword"


async def _register_and_login(
    client: AsyncClient,
    email: str,
    password: str = "test12345",
) -> tuple[str, str]:
    """Register a user and return ``(display_name, access_token)``."""
    reg = await client.post(
        "/v1/auth/register", json={"email": email, "password": password}
    )
    assert reg.status_code == 201, reg.text
    display_name = reg.json()["display_name"]
    login = await client.post(
        "/v1/auth/jwt/login",
        data={"username": email, "password": password},
    )
    assert login.status_code == 200, login.text
    return display_name, login.json()["access_token"]


async def _login_admin(client: AsyncClient) -> tuple[str, str]:
    """Login as the pre-seeded admin superuser; return ``(display_name, token)``."""
    login = await client.post(
        "/v1/auth/jwt/login",
        data={"username": ADMIN_EMAIL, "password": ADMIN_PASSWORD},
    )
    assert login.status_code == 200, login.text
    token = login.json()["access_token"]
    me = await client.get(
        "/v1/auth/users/me", headers={"Authorization": f"Bearer {token}"}
    )
    assert me.status_code == 200, me.text
    return me.json()["display_name"], token


@pytest.mark.asyncio
async def test_transfer_collision_returns_409(
    http_client_auth: AsyncClient,
) -> None:
    """Transfer to a user who already has a room with that name returns 409."""
    admin_name, admin_token = await _login_admin(http_client_auth)
    dst_name, dst_token = await _register_and_login(
        http_client_auth, "tx-dst-409@example.com"
    )

    admin_headers = {"Authorization": f"Bearer {admin_token}"}
    dst_headers = {"Authorization": f"Bearer {dst_token}"}

    # Admin creates room "collision" in their own namespace
    r = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": admin_name, "name": "collision"},
        headers=admin_headers,
    )
    assert r.status_code == 201, r.text

    # Destination user creates their own room "collision"
    r = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": dst_name, "name": "collision"},
        headers=dst_headers,
    )
    assert r.status_code == 201, r.text

    # Admin attempts to transfer their "collision" room to dst -> 409
    resp = await http_client_auth.patch(
        f"/v1/rooms/{admin_name}/collision",
        json={"new_owner": dst_name},
        headers=admin_headers,
    )
    assert resp.status_code == 409, resp.text


@pytest.mark.asyncio
async def test_transfer_happy_path(
    http_client_auth: AsyncClient,
) -> None:
    """Transfer to managed group keeps surrogate UUID; new composed address returned."""
    caller_name, caller_token = await _register_and_login(
        http_client_auth, "tx-happy@example.com"
    )
    headers = {"Authorization": f"Bearer {caller_token}"}

    # Caller creates a group they own
    g = await http_client_auth.post(
        "/v1/groups", json={"name": "tx-grp"}, headers=headers
    )
    assert g.status_code == 201, g.text
    group_name = g.json()["name"]

    # Caller creates room "t" in own namespace
    r = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": caller_name, "name": "t"},
        headers=headers,
    )
    assert r.status_code == 201, r.text

    # Transfer to the group, changing visibility to group
    resp = await http_client_auth.patch(
        f"/v1/rooms/{caller_name}/t",
        json={"new_owner": group_name, "visibility": "group"},
        headers=headers,
    )
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["room_id"] == f"{group_name}/t"

    # New composed address is reachable
    r = await http_client_auth.get(f"/v1/rooms/{group_name}/t", headers=headers)
    assert r.status_code == 200

    # Old composed address is NOT reachable
    r = await http_client_auth.get(f"/v1/rooms/{caller_name}/t", headers=headers)
    assert r.status_code == 404
