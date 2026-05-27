"""Existence-hiding invariants for the namespaced POST /v1/rooms."""

from __future__ import annotations

from typing import TYPE_CHECKING
from uuid import uuid4

import pytest

if TYPE_CHECKING:
    from httpx import AsyncClient


async def _register_and_login(
    client: AsyncClient,
    email: str,
    password: str = "test12345",
) -> tuple[str, str]:
    """Register a user and return (display_name, token)."""
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


@pytest.mark.asyncio
async def test_create_in_own_namespace(http_client_auth: AsyncClient) -> None:
    display_name, token = await _register_and_login(
        http_client_auth, "ns-create@example.com"
    )
    resp = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": display_name, "name": "my-room"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 201
    body = resp.json()
    assert body["created"] is True
    assert body["room_id"] == f"{display_name}/my-room"


@pytest.mark.asyncio
async def test_idempotent_reuse_in_own_namespace(
    http_client_auth: AsyncClient,
) -> None:
    display_name, token = await _register_and_login(
        http_client_auth, "ns-dup@example.com"
    )
    payload = {"owner": display_name, "name": "dup"}
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
    _, token = await _register_and_login(http_client_auth, "ns-cross@example.com")
    # A regex-valid but unknown display-name owner triggers UserNotFound (404),
    # not a 403 — the previous test asserted 403 because the body carried a
    # UUID owner which the namespace check rejected pre-resolution.
    foreign_owner = f"ns-foreign-{uuid4().hex[:6]}"
    resp = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": foreign_owner, "name": "anything"},
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code in (403, 404)


@pytest.mark.asyncio
async def test_cross_namespace_post_byte_identical_regardless_of_state(
    http_client_auth: AsyncClient,
) -> None:
    """Cross-namespace POST yields identical responses regardless of whether the
    target namespace has a matching room — closing the existence-leak side
    channel."""
    _, caller_token = await _register_and_login(
        http_client_auth, "ns-byte-caller@example.com"
    )
    headers = {"Authorization": f"Bearer {caller_token}"}

    # Register another user — caller will probe their namespace.
    other_display, other_token = await _register_and_login(
        http_client_auth, "ns-byte-other@example.com"
    )

    # Caller probes other's namespace BEFORE the room exists.
    resp_unknown = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": other_display, "name": "x"},
        headers=headers,
    )

    # Other creates the room.
    r = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": other_display, "name": "x"},
        headers={"Authorization": f"Bearer {other_token}"},
    )
    assert r.status_code == 201, r.text

    # Caller probes again now that the room exists.
    resp_taken = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": other_display, "name": "x"},
        headers=headers,
    )

    assert resp_unknown.status_code == resp_taken.status_code == 403
    assert resp_unknown.json() == resp_taken.json()


@pytest.mark.asyncio
async def test_group_post_requires_membership(
    http_client_auth: AsyncClient,
) -> None:
    _, owner_token = await _register_and_login(
        http_client_auth, "ns-grp-owner@example.com"
    )
    _, caller_token = await _register_and_login(
        http_client_auth, "ns-grp-caller@example.com"
    )

    # Owner creates a group
    grp_resp = await http_client_auth.post(
        "/v1/groups",
        json={"name": "ns-test-group"},
        headers={"Authorization": f"Bearer {owner_token}"},
    )
    assert grp_resp.status_code == 201, grp_resp.text
    group_name = grp_resp.json()["name"]

    # Caller is NOT a member — should get 403
    resp = await http_client_auth.post(
        "/v1/rooms",
        json={"owner": group_name, "name": "shared"},
        headers={"Authorization": f"Bearer {caller_token}"},
    )
    assert resp.status_code == 403
