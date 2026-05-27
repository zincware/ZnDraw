"""Tests for authentication endpoints (zndraw-auth / fastapi-users)."""

import pytest
from httpx import AsyncClient
from sqlmodel.ext.asyncio.session import AsyncSession

from zndraw_auth import User

# =============================================================================
# Login Tests
# =============================================================================


@pytest.mark.asyncio
async def test_login(client: AsyncClient, test_user: User) -> None:
    """Test login with valid credentials returns JWT."""
    response = await client.post(
        "/v1/auth/jwt/login",
        data={"username": "testuser@local.test", "password": "testpassword"},
    )
    assert response.status_code == 200

    body = response.json()
    assert body["token_type"] == "bearer"
    assert body["access_token"]


@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("username", "password"),
    [
        ("testuser@local.test", "wrongpassword"),
        ("nonexistent@local.test", "password"),
    ],
    ids=["invalid_password", "nonexistent_user"],
)
async def test_login_fails(
    client: AsyncClient, test_user: User, username: str, password: str
) -> None:
    """Login returns 400 for invalid credentials."""
    response = await client.post(
        "/v1/auth/jwt/login",
        data={"username": username, "password": password},
    )
    assert response.status_code == 400


# =============================================================================
# Guest Session Tests
# =============================================================================


@pytest.mark.asyncio
async def test_guest_session(client: AsyncClient) -> None:
    """Test creating a guest session returns JWT."""
    response = await client.post("/v1/auth/guest")
    assert response.status_code == 200

    body = response.json()
    assert body["token_type"] == "bearer"
    assert body["access_token"]


@pytest.mark.asyncio
async def test_guest_sessions_unique_tokens(client: AsyncClient) -> None:
    """Test that multiple guest sessions get unique tokens."""
    response1 = await client.post("/v1/auth/guest")
    response2 = await client.post("/v1/auth/guest")

    assert response1.status_code == 200
    assert response2.status_code == 200

    assert response1.json()["access_token"] != response2.json()["access_token"]


@pytest.mark.asyncio
async def test_guest_user_has_is_guest_true(client: AsyncClient, session) -> None:
    """POST /v1/auth/guest must mark the created user with is_guest=True."""
    from sqlmodel import select

    from zndraw_auth import User

    r = await client.post("/v1/auth/guest")
    assert r.status_code == 200
    email = r.json()["email"]
    result = await session.exec(select(User).where(User.email == email))
    user = result.one()
    assert user.is_guest is True


# =============================================================================
# Registration Tests
# =============================================================================


@pytest.mark.asyncio
async def test_register(client: AsyncClient) -> None:
    """Test registering a new user."""
    response = await client.post(
        "/v1/auth/register",
        json={"email": "newuser@example.com", "password": "newpassword"},
    )
    assert response.status_code == 201

    body = response.json()
    assert body["email"] == "newuser@example.com"
    assert "id" in body


@pytest.mark.asyncio
@pytest.mark.parametrize(
    ("setup_email", "email", "password", "expected_status"),
    [
        ("duplicate@example.com", "duplicate@example.com", "password456", 400),
        (None, "newuser@example.com", None, 422),
    ],
    ids=["duplicate_email", "missing_password"],
)
async def test_register_fails(
    client: AsyncClient,
    setup_email: str | None,
    email: str,
    password: str | None,
    expected_status: int,
) -> None:
    """Registration returns error for invalid input."""
    if setup_email:
        setup_resp = await client.post(
            "/v1/auth/register",
            json={"email": setup_email, "password": "password123"},
        )
        assert setup_resp.status_code == 201
    body = {"email": email}
    if password is not None:
        body["password"] = password
    response = await client.post("/v1/auth/register", json=body)
    assert response.status_code == expected_status


# =============================================================================
# is_verified Permission Tests
# =============================================================================


@pytest.mark.asyncio
async def test_non_superuser_cannot_set_is_verified(
    client: AsyncClient, session
) -> None:
    """Non-superusers must not be able to flip is_verified on themselves.

    fastapi-users' ``safe=True`` default on the users router drops privileged
    fields from non-superuser updates — this test catches regressions if
    that behavior changes (e.g., someone re-registers the router with
    ``safe=False``).
    """
    from helpers import create_test_user_in_db

    # Bypass dev-mode auto-promote; helper creates users with is_verified=True
    user, token = await create_test_user_in_db(
        session, email="notsu@test.com", is_superuser=False
    )
    original_verified = user.is_verified  # True by default in the helper

    # Attempt to flip is_verified to the opposite value
    flipped = not original_verified
    r = await client.patch(
        "/v1/auth/users/me",
        json={"is_verified": flipped},
        headers={"Authorization": f"Bearer {token}"},
    )
    # fastapi-users silently drops unsafe fields for non-superusers, so the
    # request itself succeeds (200), but is_verified is unchanged.
    assert r.status_code == 200

    me = (
        await client.get(
            "/v1/auth/users/me",
            headers={"Authorization": f"Bearer {token}"},
        )
    ).json()
    assert me["is_verified"] is original_verified


@pytest.mark.asyncio
async def test_superuser_can_set_is_verified(client: AsyncClient, session) -> None:
    """Superusers can verify other users via PATCH."""
    from helpers import create_test_user_in_db

    _admin_user, admin_token = await create_test_user_in_db(
        session, email="admin-verify@test.com", is_superuser=True
    )
    # Helper creates users with is_verified=True; flip to False to test the
    # superuser's ability to change the field.
    target_user, _target_token = await create_test_user_in_db(
        session, email="target-verify@test.com", is_superuser=False
    )
    assert target_user.is_verified is True

    r = await client.patch(
        f"/v1/auth/users/{target_user.id}",
        json={"is_verified": False},
        headers={"Authorization": f"Bearer {admin_token}"},
    )
    assert r.status_code == 200
    assert r.json()["is_verified"] is False


# =============================================================================
# Display Name Suggestion Tests
# =============================================================================


@pytest.mark.asyncio
async def test_available_display_name_returns_valid_slug(
    client: AsyncClient,
) -> None:
    resp = await client.get("/v1/users/available-display-name")
    assert resp.status_code == 200, resp.text
    body = resp.json()
    from zndraw_auth.display_names import DISPLAY_NAME_PATTERN
    assert DISPLAY_NAME_PATTERN.fullmatch(body["display_name"])


@pytest.mark.asyncio
async def test_available_display_name_never_collides_with_existing(
    client: AsyncClient, session: AsyncSession
) -> None:
    # Register one user, then ask for a suggestion — must differ.
    reg = await client.post(
        "/v1/auth/register",
        json={
            "email": "eve@example.com",
            "password": "very-strong-passw0rd",
            "display_name": "eve-the-curious",
        },
    )
    assert reg.status_code == 201, reg.text
    resp = await client.get("/v1/users/available-display-name")
    assert resp.json()["display_name"] != "eve-the-curious"
