"""ZnDraw client accepts and parses composed room addresses."""

from __future__ import annotations

import contextlib

import httpx
import pytest

from zndraw import ZnDraw


def _get_guest_display_name(server: str) -> str:
    """Mint a guest token and return that user's display_name."""
    with httpx.Client(base_url=server) as client:
        token_resp = client.post("/v1/auth/guest")
        token_resp.raise_for_status()
        token = token_resp.json()["access_token"]
        me_resp = client.get(
            "/v1/auth/users/me",
            headers={"Authorization": f"Bearer {token}"},
        )
        me_resp.raise_for_status()
        return me_resp.json()["display_name"]


def test_client_rejects_single_segment(server: str) -> None:
    with pytest.raises(ValueError, match="composed form"):
        ZnDraw(url=server, room="single-segment")


def test_client_accepts_composed(server: str) -> None:
    display_name = _get_guest_display_name(server)
    room_address = f"{display_name}/my-room"
    vis = ZnDraw(url=server, room=room_address)
    assert vis.room == room_address
    with contextlib.suppress(Exception):
        vis.disconnect()


def test_client_rejects_bad_owner_uuid(server: str) -> None:
    # A bare UUID is no longer a valid owner segment — owner must be display_name.
    with pytest.raises(ValueError, match="not a valid display name"):
        ZnDraw(url=server, room="00000000-0000-0000-0000-000000000000/some-name")


def test_client_rejects_bad_room_name(server: str) -> None:
    display_name = _get_guest_display_name(server)
    with pytest.raises(ValueError, match="invalid characters"):
        ZnDraw(url=server, room=f"{display_name}/has spaces")
