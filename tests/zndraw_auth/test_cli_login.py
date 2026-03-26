"""Tests for CLI login (device-code) flow."""

from datetime import UTC, datetime, timedelta
from unittest.mock import patch

import pytest
from httpx import AsyncClient

from zndraw_auth import TokenResponse, UserCreate


@pytest.mark.asyncio
async def test_create_cli_login_challenge(client: AsyncClient) -> None:
    """POST /auth/cli-login creates a challenge with code and secret."""
    response = await client.post("/auth/cli-login")
    assert response.status_code == 200
    data = response.json()
    assert "code" in data
    assert "secret" in data
    assert "approve_url" in data
    assert len(data["code"]) == 8


@pytest.mark.asyncio
async def test_poll_pending_challenge(client: AsyncClient) -> None:
    """GET /auth/cli-login/{code} returns pending when not yet approved."""
    create = await client.post("/auth/cli-login")
    data = create.json()

    response = await client.get(
        f"/auth/cli-login/{data['code']}",
        params={"secret": data["secret"]},
    )
    assert response.status_code == 200
    assert response.json()["status"] == "pending"
    assert response.json().get("token") is None


@pytest.mark.asyncio
async def test_poll_wrong_secret(client: AsyncClient) -> None:
    """GET /auth/cli-login/{code} with wrong secret returns 404."""
    create = await client.post("/auth/cli-login")
    data = create.json()

    response = await client.get(
        f"/auth/cli-login/{data['code']}",
        params={"secret": "wrong-secret"},
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_poll_nonexistent_code(client: AsyncClient) -> None:
    """GET /auth/cli-login/{code} with unknown code returns 404."""
    response = await client.get(
        "/auth/cli-login/ZZZZZZZZ",
        params={"secret": "doesnt-matter"},
    )
    assert response.status_code == 404


async def _get_auth_header(
    client: AsyncClient, email: str, password: str
) -> dict[str, str]:
    """Register, login, return auth header."""
    user_data = UserCreate(email=email, password=password)
    await client.post("/auth/register", json=user_data.model_dump())
    response = await client.post(
        "/auth/jwt/login",
        data={
            "username": email,
            "password": password,
            "grant_type": "password",
        },
    )
    token = TokenResponse.model_validate(response.json())
    return {"Authorization": f"Bearer {token.access_token}"}


@pytest.mark.asyncio
async def test_approve_challenge(client: AsyncClient) -> None:
    """PATCH /auth/cli-login/{code} approves and mints a token."""
    create = await client.post("/auth/cli-login")
    data = create.json()

    headers = await _get_auth_header(client, "cli-user@test.com", "password123")
    response = await client.patch(
        f"/auth/cli-login/{data['code']}",
        headers=headers,
    )
    assert response.status_code == 200


@pytest.mark.asyncio
async def test_approve_requires_auth(client: AsyncClient) -> None:
    """PATCH /auth/cli-login/{code} without auth returns 401."""
    create = await client.post("/auth/cli-login")
    data = create.json()

    response = await client.patch(f"/auth/cli-login/{data['code']}")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_full_cli_login_flow(client: AsyncClient) -> None:
    """Full flow: create -> approve -> poll gets token -> token works."""
    create = await client.post("/auth/cli-login")
    challenge = create.json()

    headers = await _get_auth_header(client, "browser@test.com", "password123")
    approve = await client.patch(
        f"/auth/cli-login/{challenge['code']}",
        headers=headers,
    )
    assert approve.status_code == 200

    poll = await client.get(
        f"/auth/cli-login/{challenge['code']}",
        params={"secret": challenge["secret"]},
    )
    assert poll.status_code == 200
    poll_data = poll.json()
    assert poll_data["status"] == "approved"
    assert poll_data["token"] is not None

    cli_headers = {"Authorization": f"Bearer {poll_data['token']}"}
    me = await client.get("/test/protected", headers=cli_headers)
    assert me.status_code == 200
    assert me.json()["email"] == "browser@test.com"


@pytest.mark.asyncio
async def test_poll_after_redeem_returns_404(client: AsyncClient) -> None:
    """Second poll after redeem returns 404 (one-time retrieval)."""
    create = await client.post("/auth/cli-login")
    challenge = create.json()
    headers = await _get_auth_header(client, "redeem@test.com", "password123")
    await client.patch(
        f"/auth/cli-login/{challenge['code']}",
        headers=headers,
    )

    poll1 = await client.get(
        f"/auth/cli-login/{challenge['code']}",
        params={"secret": challenge["secret"]},
    )
    assert poll1.status_code == 200
    assert poll1.json()["status"] == "approved"

    poll2 = await client.get(
        f"/auth/cli-login/{challenge['code']}",
        params={"secret": challenge["secret"]},
    )
    assert poll2.status_code == 404


@pytest.mark.asyncio
async def test_reject_challenge(client: AsyncClient) -> None:
    """DELETE /auth/cli-login/{code} rejects, subsequent poll returns 404."""
    create = await client.post("/auth/cli-login")
    challenge = create.json()

    headers = await _get_auth_header(client, "rejector@test.com", "password123")
    response = await client.delete(
        f"/auth/cli-login/{challenge['code']}",
        headers=headers,
    )
    assert response.status_code == 204

    poll = await client.get(
        f"/auth/cli-login/{challenge['code']}",
        params={"secret": challenge["secret"]},
    )
    assert poll.status_code == 404


@pytest.mark.asyncio
async def test_reject_requires_auth(client: AsyncClient) -> None:
    """DELETE /auth/cli-login/{code} without auth returns 401."""
    create = await client.post("/auth/cli-login")
    challenge = create.json()

    response = await client.delete(f"/auth/cli-login/{challenge['code']}")
    assert response.status_code == 401


@pytest.mark.asyncio
async def test_poll_expired_challenge(client: AsyncClient) -> None:
    """GET /auth/cli-login/{code} returns 410 when challenge is expired."""
    # Create challenge (stores real timestamps)
    create = await client.post("/auth/cli-login")
    challenge = create.json()

    # Mock datetime.now to return future time for the poll
    future = datetime.now(UTC).replace(tzinfo=None) + timedelta(minutes=10)
    with patch("zndraw_auth.cli_login.datetime") as mock_dt:
        mock_dt.now.return_value = future

        response = await client.get(
            f"/auth/cli-login/{challenge['code']}",
            params={"secret": challenge["secret"]},
        )
    assert response.status_code == 410


@pytest.mark.asyncio
async def test_routers_importable_from_package() -> None:
    """cli_login_router and admin_token_router are importable from zndraw_auth."""
    from zndraw_auth import (
        CLILoginChallenge,
        CLILoginCreateResponse,
        CLILoginStatusResponse,
        ImpersonationTokenResponse,
        admin_token_router,
        cli_login_router,
    )

    assert cli_login_router is not None
    assert admin_token_router is not None
    assert CLILoginChallenge is not None
    assert CLILoginCreateResponse is not None
    assert CLILoginStatusResponse is not None
    assert ImpersonationTokenResponse is not None


@pytest.mark.asyncio
async def test_approve_expired_challenge(client: AsyncClient) -> None:
    """PATCH /auth/cli-login/{code} returns 410 when challenge is expired."""
    create = await client.post("/auth/cli-login")
    challenge = create.json()

    headers = await _get_auth_header(client, "expired@test.com", "password123")

    future = datetime.now(UTC).replace(tzinfo=None) + timedelta(minutes=10)
    with patch("zndraw_auth.cli_login.datetime") as mock_dt:
        mock_dt.now.return_value = future

        response = await client.patch(
            f"/auth/cli-login/{challenge['code']}",
            headers=headers,
        )
    assert response.status_code == 410
