"""Regression: pyclient ZnDraw refreshes cached_length on FramesInvalidate (finding #3).

Stand up a real server via ``server_factory``, register an admin user, create
a room, and spawn a ``ZnDraw`` client pointed at that room's composed address.
A second authenticated HTTP request appends a frame. The first client's
``cached_length`` must update within a short timeout.
"""

from __future__ import annotations

import time
from collections.abc import Callable
from typing import TYPE_CHECKING

import ase
import httpx

from zndraw.client import atoms_to_json_dict

if TYPE_CHECKING:
    from conftest import ServerInstance

ServerFactory = Callable[[dict[str, str]], "ServerInstance"]


def _make_json_frame() -> dict:
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    return atoms_to_json_dict(atoms)


def test_pyclient_cached_length_updates_on_frame_append(
    server_factory: ServerFactory,
) -> None:
    from zndraw import ZnDraw

    instance = server_factory(
        {
            "ZNDRAW_AUTH_DEFAULT_ADMIN_EMAIL": "admin@local.test",
            "ZNDRAW_AUTH_DEFAULT_ADMIN_PASSWORD": "adminpassword",
        }
    )
    base_url = instance.url

    login = httpx.post(
        f"{base_url}/v1/auth/jwt/login",
        data={"username": "admin@local.test", "password": "adminpassword"},
        timeout=5.0,
    )
    login.raise_for_status()
    token = login.json()["access_token"]
    headers = {"Authorization": f"Bearer {token}"}

    user_resp = httpx.get(f"{base_url}/v1/auth/users/me", headers=headers, timeout=5.0)
    user_resp.raise_for_status()
    user_id = user_resp.json()["id"]

    create = httpx.post(
        f"{base_url}/v1/rooms",
        json={
            "owner_id": user_id,
            "name": "pytest-room",
            "visibility": "public",
            "copy_from": "@none",
        },
        headers=headers,
        timeout=5.0,
    )
    create.raise_for_status()
    composed_room = create.json()["room_id"]

    client = ZnDraw(url=base_url, room=composed_room, token=token)
    try:
        # connect explicitly so the socket is live before we read len
        client.connect()
        assert len(client) == 0

        append = httpx.post(
            f"{base_url}/v1/rooms/{composed_room}/frames",
            json={"frames": [_make_json_frame()]},
            headers=headers,
            timeout=5.0,
        )
        append.raise_for_status()

        deadline = time.monotonic() + 3.0
        while time.monotonic() < deadline:
            if client.cached_length == 1:
                break
            time.sleep(0.05)
        assert client.cached_length == 1, (
            f"cached_length did not refresh; got {client.cached_length}"
        )
    finally:
        client.disconnect()
