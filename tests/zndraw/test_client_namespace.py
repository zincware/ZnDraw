"""ZnDraw client accepts and parses composed room addresses."""

from __future__ import annotations

from uuid import uuid4

import httpx
import pytest

from zndraw import ZnDraw


def _get_guest_user_id(server: str) -> str:
    """Mint a guest token and return that user's ID."""
    with httpx.Client(base_url=server) as client:
        token_resp = client.post("/v1/auth/guest")
        token_resp.raise_for_status()
        token = token_resp.json()["access_token"]
        me_resp = client.get(
            "/v1/auth/users/me",
            headers={"Authorization": f"Bearer {token}"},
        )
        me_resp.raise_for_status()
        return me_resp.json()["id"]


def test_client_rejects_single_segment(server: str) -> None:
    with pytest.raises(ValueError, match="composed form"):
        ZnDraw(url=server, room="single-segment")


def test_client_accepts_composed(server: str) -> None:
    user_id = _get_guest_user_id(server)
    room_address = f"{user_id}/my-room"
    vis = ZnDraw(url=server, room=room_address)
    assert vis.room == room_address
    try:
        vis.disconnect()
    except Exception:
        pass


def test_client_rejects_bad_owner_uuid(server: str) -> None:
    with pytest.raises(ValueError, match="UUID"):
        ZnDraw(url=server, room="not-a-uuid/some-name")


def test_client_rejects_bad_room_name(server: str) -> None:
    user_id = _get_guest_user_id(server)
    with pytest.raises(ValueError, match="invalid characters"):
        ZnDraw(url=server, room=f"{user_id}/has spaces")
