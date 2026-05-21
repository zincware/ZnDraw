"""GET /v1/rooms emits composed addresses and resolved owner labels."""

from __future__ import annotations

import pytest
from httpx import AsyncClient


async def _register_and_login(client: AsyncClient, email: str) -> tuple[str, str]:
    r = await client.post(
        "/v1/auth/register",
        json={"email": email, "password": "password123"},
    )
    user_id = r.json()["id"]
    r = await client.post(
        "/v1/auth/jwt/login",
        data={"username": email, "password": "password123"},
    )
    return user_id, r.json()["access_token"]


@pytest.mark.asyncio
async def test_list_returns_composed_room_id(
    http_client_auth: AsyncClient,
) -> None:
    user_id, token = await _register_and_login(
        http_client_auth, "list-user@example.com"
    )
    headers = {"Authorization": f"Bearer {token}"}
    create = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": user_id, "name": "r1"},
        headers=headers,
    )
    assert create.status_code == 201, create.text

    resp = await http_client_auth.get("/v1/rooms", headers=headers)
    assert resp.status_code == 200
    items = resp.json()["items"]
    item = next(it for it in items if it["room_id"] == f"{user_id}/r1")
    assert item["owner_kind"] == "user"
    assert item["owner_label"] == "list-user@example.com"
    assert item["owner_id"] == user_id


@pytest.mark.asyncio
async def test_get_room_returns_composed_address(
    http_client_auth: AsyncClient,
) -> None:
    user_id, token = await _register_and_login(http_client_auth, "get-user@example.com")
    headers = {"Authorization": f"Bearer {token}"}
    create = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": user_id, "name": "getme"},
        headers=headers,
    )
    assert create.status_code == 201

    resp = await http_client_auth.get(f"/v1/rooms/{user_id}/getme", headers=headers)
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["room_id"] == f"{user_id}/getme"
    assert body["owner_kind"] == "user"
