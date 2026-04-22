"""End-to-end scope refactor smoke test hitting a real uvicorn server.

Uses ``server_auth`` fixture so dev-mode auto-promote is OFF; registered
users are ordinary active users.
"""

import pytest
from httpx import AsyncClient


async def _register_and_login(
    client: AsyncClient, email: str, password: str = "test12345"
) -> str:
    """Register a user and return a JWT access token.

    Parameters
    ----------
    client
        Async HTTP client pointing at the test server.
    email
        Email address to register.
    password
        Password to use for registration and login.

    Returns
    -------
    str
        JWT access token for the registered user.
    """
    await client.post("/v1/auth/register", json={"email": email, "password": password})
    r = await client.post(
        "/v1/auth/jwt/login",
        data={"username": email, "password": password},
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    r.raise_for_status()
    return r.json()["access_token"]


@pytest.mark.asyncio
async def test_full_group_workflow(server_auth: str) -> None:
    """Exercise the full group/room/share-link permission workflow.

    Steps
    -----
    1. Admin registers, creates a group, adds a member (role=MEMBER).
    2. Admin creates a group-owned room (visibility=group).
    3. Member can read the group room.
    4. Member cannot PATCH (manage-gated) → 403.
    5. Outsider cannot see the group room → 404.
    6. Share link (view) lets outsider read the room.
    7. Admin PATCHes visibility to public.
    8. Outsider now sees it without the token.
    9. Admin revokes the share link → 204.

    Parameters
    ----------
    server_auth
        Base URL of the running auth-enabled uvicorn server.
    """
    async with AsyncClient(base_url=server_auth) as client:
        admin = await _register_and_login(client, "e2e-admin@test.com")
        member = await _register_and_login(client, "e2e-mem@test.com")
        outsider = await _register_and_login(client, "e2e-out@test.com")

        # Create group + add member (role=MEMBER to enable edit, not manage)
        gid = (
            await client.post(
                "/v1/groups",
                json={"name": "e2e-team"},
                headers={"Authorization": f"Bearer {admin}"},
            )
        ).json()["id"]
        me = (
            await client.get(
                "/v1/auth/users/me", headers={"Authorization": f"Bearer {member}"}
            )
        ).json()
        add_r = await client.post(
            f"/v1/groups/{gid}/members",
            json={"user_id": me["id"], "role": "member"},
            headers={"Authorization": f"Bearer {admin}"},
        )
        assert add_r.status_code == 201

        # Group-owned room
        create_r = await client.post(
            "/v1/rooms",
            json={"room_id": "e2e-grp", "visibility": "group", "owner_group_id": gid},
            headers={"Authorization": f"Bearer {admin}"},
        )
        assert create_r.status_code == 201

        # Member can read the room
        r = await client.get(
            "/v1/rooms/e2e-grp", headers={"Authorization": f"Bearer {member}"}
        )
        assert r.status_code == 200

        # Member cannot PATCH (manage-gated) — expect 403
        r = await client.patch(
            "/v1/rooms/e2e-grp",
            json={"description": "from member"},
            headers={"Authorization": f"Bearer {member}"},
        )
        assert r.status_code == 403

        # Outsider cannot see the room (404 — existence hidden)
        r = await client.get(
            "/v1/rooms/e2e-grp",
            headers={"Authorization": f"Bearer {outsider}"},
        )
        assert r.status_code == 404

        # Share link (view) lets outsider read
        link = (
            await client.post(
                "/v1/rooms/e2e-grp/share-links",
                json={"access": "view"},
                headers={"Authorization": f"Bearer {admin}"},
            )
        ).json()
        r = await client.get(
            "/v1/rooms/e2e-grp",
            headers={
                "Authorization": f"Bearer {outsider}",
                "X-Room-Share-Token": link["token"],
            },
        )
        assert r.status_code == 200

        # Admin can PATCH visibility (manage role)
        r = await client.patch(
            "/v1/rooms/e2e-grp",
            json={"visibility": "public"},
            headers={"Authorization": f"Bearer {admin}"},
        )
        assert r.status_code == 200

        # Outsider now sees it without the token
        r = await client.get(
            "/v1/rooms/e2e-grp",
            headers={"Authorization": f"Bearer {outsider}"},
        )
        assert r.status_code == 200
        assert r.json()["visibility"] == "public"

        # Admin revokes the link
        r = await client.delete(
            f"/v1/rooms/e2e-grp/share-links/{link['id']}",
            headers={"Authorization": f"Bearer {admin}"},
        )
        assert r.status_code == 204
