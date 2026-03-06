import httpx
import pytest
import socketio
import socketio.exceptions
from zndraw_socketio import wrap

from zndraw.socket_events import UserGet, UserGetResponse


@pytest.mark.asyncio
async def test_socketio_connect_without_token(server: str) -> None:
    """Test that Socket.IO connection is rejected without a token."""
    sio_client = socketio.AsyncClient()

    with pytest.raises(socketio.exceptions.ConnectionError):
        await sio_client.connect(server)


@pytest.mark.asyncio
async def test_socketio_connect_with_invalid_token(server: str) -> None:
    """Test that Socket.IO connection is rejected with invalid token."""
    sio_client = socketio.AsyncClient()

    with pytest.raises(socketio.exceptions.ConnectionError):
        await sio_client.connect(server, auth={"token": "invalid-token"})


@pytest.mark.asyncio
async def test_socketio_get_user_authenticated(server_auth: str) -> None:
    """Test that get_user returns correct data for authenticated user."""
    async with httpx.AsyncClient() as client:
        # Register a user
        register_response = await client.post(
            f"{server_auth}/v1/auth/register",
            json={"email": "testuser@example.com", "password": "testpassword"},
        )
        assert register_response.status_code == 201

        # Login to get a token (OAuth2 form data)
        login_response = await client.post(
            f"{server_auth}/v1/auth/jwt/login",
            data={"username": "testuser@example.com", "password": "testpassword"},
        )
        assert login_response.status_code == 200
        token = login_response.json()["access_token"]

        # Connect to Socket.IO with auth payload
        sio_client = socketio.AsyncClient()
        await sio_client.connect(server_auth, auth={"token": token})
        assert sio_client.connected

        # Call user_get via typed wrapper
        tsio_client = wrap(sio_client)
        response = await tsio_client.call(UserGet(), response_model=UserGetResponse)
        assert response.email == "testuser@example.com"
        assert response.is_superuser is False

        await sio_client.disconnect()


@pytest.mark.asyncio
async def test_socketio_get_user_guest(server_auth: str) -> None:
    """Test that get_user returns correct data for guest user."""
    async with httpx.AsyncClient() as client:
        # Create a guest token
        guest_response = await client.post(f"{server_auth}/v1/auth/guest")
        assert guest_response.status_code == 200
        token = guest_response.json()["access_token"]

        # Connect to Socket.IO with auth payload
        sio_client = socketio.AsyncClient()
        await sio_client.connect(server_auth, auth={"token": token})
        assert sio_client.connected

        # Call user_get via typed wrapper
        tsio_client = wrap(sio_client)
        response = await tsio_client.call(UserGet(), response_model=UserGetResponse)
        assert response.email.endswith("@guest.user")
        assert response.is_superuser is False

        await sio_client.disconnect()
