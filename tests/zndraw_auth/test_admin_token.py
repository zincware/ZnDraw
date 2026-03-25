"""Tests for admin token minting."""

import pytest
from httpx import AsyncClient

from zndraw_auth import TokenResponse, UserCreate, UserRead


async def _get_auth_header(
    client: AsyncClient, email: str, password: str
) -> dict[str, str]:
    """Register, login, return auth header."""
    user_data = UserCreate(email=email, password=password)
    await client.post("/auth/register", json=user_data.model_dump())
    response = await client.post(
        "/auth/jwt/login",
        data={"username": email, "password": password, "grant_type": "password"},
    )
    token = TokenResponse.model_validate(response.json())
    return {"Authorization": f"Bearer {token.access_token}"}


async def _get_admin_header(client: AsyncClient) -> dict[str, str]:
    """Login as the pre-created admin (admin@test.com / admin-password)."""
    response = await client.post(
        "/auth/jwt/login",
        data={
            "username": "admin@test.com",
            "password": "admin-password",
            "grant_type": "password",
        },
    )
    token = TokenResponse.model_validate(response.json())
    return {"Authorization": f"Bearer {token.access_token}"}


@pytest.mark.asyncio
async def test_admin_mint_token_for_user(client: AsyncClient) -> None:
    """Superuser can mint a token for another user."""
    user_data = UserCreate(email="target@test.com", password="password123")
    reg = await client.post("/auth/register", json=user_data.model_dump())
    target_user = UserRead.model_validate(reg.json())

    admin_headers = await _get_admin_header(client)
    response = await client.post(
        f"/admin/users/{target_user.id}/token",
        headers=admin_headers,
    )
    assert response.status_code == 200
    data = response.json()
    assert "access_token" in data
    assert data["token_type"] == "bearer"

    # The minted token works as the target user
    cli_headers = {"Authorization": f"Bearer {data['access_token']}"}
    me = await client.get("/test/protected", headers=cli_headers)
    assert me.status_code == 200
    assert me.json()["email"] == "target@test.com"


@pytest.mark.asyncio
async def test_admin_mint_token_non_superuser_forbidden(
    client: AsyncClient,
) -> None:
    """Non-superuser cannot mint tokens."""
    headers = await _get_auth_header(client, "regular@test.com", "password123")

    user_data = UserCreate(email="other@test.com", password="password123")
    reg = await client.post("/auth/register", json=user_data.model_dump())
    target = UserRead.model_validate(reg.json())

    response = await client.post(
        f"/admin/users/{target.id}/token",
        headers=headers,
    )
    assert response.status_code == 403


@pytest.mark.asyncio
async def test_admin_mint_token_target_not_found(client: AsyncClient) -> None:
    """Minting token for nonexistent user returns 404."""
    admin_headers = await _get_admin_header(client)
    fake_id = "00000000-0000-0000-0000-000000000000"
    response = await client.post(
        f"/admin/users/{fake_id}/token",
        headers=admin_headers,
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_admin_mint_token_has_impersonated_by_claim(
    client: AsyncClient,
) -> None:
    """Minted token contains impersonated_by claim with admin's UUID."""
    import jwt as pyjwt

    user_data = UserCreate(email="claimed@test.com", password="password123")
    reg = await client.post("/auth/register", json=user_data.model_dump())
    target = UserRead.model_validate(reg.json())

    admin_headers = await _get_admin_header(client)
    response = await client.post(
        f"/admin/users/{target.id}/token",
        headers=admin_headers,
    )
    data = response.json()

    claims = pyjwt.decode(
        data["access_token"],
        options={"verify_signature": False},
    )
    assert claims["sub"] == str(target.id)
    assert "impersonated_by" in claims


@pytest.mark.asyncio
async def test_admin_mint_token_unauthenticated(client: AsyncClient) -> None:
    """Unauthenticated request to admin endpoint returns 401."""
    fake_id = "00000000-0000-0000-0000-000000000000"
    response = await client.post(f"/admin/users/{fake_id}/token")
    assert response.status_code == 401
