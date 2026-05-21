"""Existence-hiding invariants for the namespaced POST /v1/rooms."""

from __future__ import annotations

from uuid import uuid4

import pytest
from httpx import AsyncClient


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


@pytest.mark.asyncio
async def test_create_in_own_namespace(http_client_auth: AsyncClient) -> None:
    user_id, token = await _register_and_login(
        http_client_auth, "ns-create@example.com"
    )
    resp = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": user_id, "name": "my-room"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 201
    body = resp.json()
    assert body["created"] is True
    assert body["room_id"] == f"{user_id}/my-room"


@pytest.mark.asyncio
async def test_idempotent_reuse_in_own_namespace(
    http_client_auth: AsyncClient,
) -> None:
    user_id, token = await _register_and_login(
        http_client_auth, "ns-dup@example.com"
    )
    payload = {"owner_id": user_id, "name": "dup"}
    headers = {"Authorization": f"Bearer {token}"}
    first = await http_client_auth.post("/v1/rooms", json=payload, headers=headers)
    second = await http_client_auth.post("/v1/rooms", json=payload, headers=headers)
    assert first.status_code == 201
    assert second.status_code == 200
    assert second.json()["created"] is False
    assert first.json()["room_id"] == second.json()["room_id"]


@pytest.mark.asyncio
async def test_cross_namespace_post_returns_403(
    http_client_auth: AsyncClient,
) -> None:
    _user_id, token = await _register_and_login(
        http_client_auth, "ns-cross@example.com"
    )
    foreign_owner = str(uuid4())
    resp = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": foreign_owner, "name": "anything"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 403


@pytest.mark.asyncio
async def test_cross_namespace_post_byte_identical_regardless_of_state(
    http_client_auth: AsyncClient,
) -> None:
    """403 body is byte-identical regardless of whether the foreign namespace
    has a matching room or not, closing the existence-leak side channel."""
    caller_id, caller_token = await _register_and_login(
        http_client_auth, "ns-byte-caller@example.com"
    )
    headers = {"Authorization": f"Bearer {caller_token}"}

    # POST to a completely unknown (non-existent) owner_id
    resp_unknown = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": str(uuid4()), "name": "x"},
        headers=headers,
    )

    # Register another user and have them create a room "x"
    other_id, other_token = await _register_and_login(
        http_client_auth, "ns-byte-other@example.com"
    )
    r = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": other_id, "name": "x"},
        headers={"Authorization": f"Bearer {other_token}"},
    )
    assert r.status_code == 201, r.text

    # Caller attempts to POST to other's namespace where "x" now exists
    resp_taken = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": other_id, "name": "x"},
        headers=headers,
    )

    assert resp_unknown.status_code == resp_taken.status_code == 403
    assert resp_unknown.json() == resp_taken.json()


@pytest.mark.asyncio
async def test_group_post_requires_membership(
    http_client_auth: AsyncClient,
) -> None:
    owner_id, owner_token = await _register_and_login(
        http_client_auth, "ns-grp-owner@example.com"
    )
    caller_id, caller_token = await _register_and_login(
        http_client_auth, "ns-grp-caller@example.com"
    )

    # Owner creates a group
    grp_resp = await http_client_auth.post(
        "/v1/groups",
        json={"name": "ns-test-group"},
        headers={"Authorization": f"Bearer {owner_token}"},
    )
    assert grp_resp.status_code == 201, grp_resp.text
    group_id = grp_resp.json()["id"]

    # Caller is NOT a member — should get 403
    resp = await http_client_auth.post(
        "/v1/rooms",
        json={"owner_id": group_id, "name": "shared"},
        headers={"Authorization": f"Bearer {caller_token}"},
    )
    assert resp.status_code == 403
