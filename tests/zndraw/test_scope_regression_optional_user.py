"""Routes that previously accepted OptionalUserDep now require auth.

The trajectory GET endpoint is intentionally excluded: it supports
anonymous downloads via temporary download tokens, so it uses optional
auth by design.
"""

import pytest
from httpx import AsyncClient

_FAKE_OWNER = "00000000-0000-0000-0000-000000000001"
_FAKE_ROOM = f"{_FAKE_OWNER}/any-room"

ENDPOINTS = [
    ("GET", "/v1/rooms"),
    ("GET", f"/v1/rooms/{_FAKE_ROOM}/geometries"),
    ("GET", f"/v1/rooms/{_FAKE_ROOM}/geometries/key/selection"),
    ("GET", f"/v1/rooms/{_FAKE_ROOM}/geometries/key"),
    ("GET", f"/v1/rooms/{_FAKE_ROOM}/figures"),
    ("GET", f"/v1/rooms/{_FAKE_ROOM}/figures/k"),
]


@pytest.mark.asyncio
@pytest.mark.parametrize(("method", "path"), ENDPOINTS)
async def test_endpoint_requires_auth(
    client: AsyncClient, method: str, path: str
) -> None:
    r = await client.request(method, path)
    assert r.status_code == 401, f"{method} {path} returned {r.status_code}"
