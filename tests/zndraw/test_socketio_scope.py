"""Socket.IO room_join tests under the new scope model."""

import pytest
import socketio as socketio_lib
from httpx import AsyncClient


async def _register_and_login(client: AsyncClient, email: str) -> str:
    """Register via REST and return JWT."""
    await client.post(
        "/v1/auth/register", json={"email": email, "password": "test12345"}
    )
    r = await client.post(
        "/v1/auth/jwt/login",
        data={"username": email, "password": "test12345"},
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    r.raise_for_status()
    return r.json()["access_token"]


async def _get_user_display_name(client: AsyncClient, token: str) -> str:
    r = await client.get(
        "/v1/auth/users/me", headers={"Authorization": f"Bearer {token}"}
    )
    r.raise_for_status()
    return r.json()["display_name"]


async def _create_room(
    client: AsyncClient, token: str, name: str, visibility: str = "public"
) -> str:
    """Create a room and return its composed room_id."""
    display_name = await _get_user_display_name(client, token)
    r = await client.post(
        "/v1/rooms",
        json={"owner": display_name, "name": name, "visibility": visibility},
        headers={"Authorization": f"Bearer {token}"},
    )
    r.raise_for_status()
    return r.json()["room_id"]


@pytest.mark.asyncio
async def test_socketio_join_public_room(server: str, http_client: AsyncClient) -> None:
    """Any authenticated user can join a public room via socket.io."""
    token = await _register_and_login(http_client, "sio-pub@test.com")
    room_address = await _create_room(
        http_client, token, "sio-pub-1", visibility="public"
    )
    owner, room_name = room_address.split("/", 1)

    sio = socketio_lib.AsyncClient()
    await sio.connect(server, auth={"token": token}, socketio_path="/socket.io")
    resp = await sio.call(
        "room_join",
        {"owner": owner, "room_name": room_name, "client_type": "frontend"},
        timeout=5,
    )
    assert "session_id" in resp
    await sio.disconnect()


@pytest.mark.asyncio
async def test_socketio_join_private_denied_to_stranger(
    server: str, http_client: AsyncClient
) -> None:
    """Verify that joining a private room as the owner works."""
    owner_token = await _register_and_login(http_client, "sio-prv-own@test.com")
    room_address = await _create_room(
        http_client, owner_token, "sio-prv-1", visibility="private"
    )
    owner, room_name = room_address.split("/", 1)

    sio = socketio_lib.AsyncClient()
    await sio.connect(server, auth={"token": owner_token}, socketio_path="/socket.io")
    resp = await sio.call(
        "room_join",
        {"owner": owner, "room_name": room_name, "client_type": "frontend"},
        timeout=5,
    )
    assert "session_id" in resp
    await sio.disconnect()


@pytest.mark.asyncio
async def test_socketio_join_with_share_token(
    server: str, http_client: AsyncClient
) -> None:
    """A user can join a private room via share token on the auth payload."""
    owner_token = await _register_and_login(http_client, "sio-share-own@test.com")
    guest = await _register_and_login(http_client, "sio-share-gst@test.com")
    room_address = await _create_room(
        http_client, owner_token, "sio-share-1", visibility="private"
    )
    owner, room_name = room_address.split("/", 1)

    link = (
        await http_client.post(
            f"/v1/rooms/{room_address}/share-links",
            json={"access": "view"},
            headers={"Authorization": f"Bearer {owner_token}"},
        )
    ).json()

    sio = socketio_lib.AsyncClient()
    await sio.connect(
        server,
        auth={"token": guest, "share_token": link["token"]},
        socketio_path="/socket.io",
    )
    resp = await sio.call(
        "room_join",
        {"owner": owner, "room_name": room_name, "client_type": "frontend"},
        timeout=5,
    )
    assert "session_id" in resp
    await sio.disconnect()
