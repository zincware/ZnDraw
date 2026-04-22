"""Every route that previously accepted OptionalUserDep now requires auth."""
import pytest
from httpx import AsyncClient


ENDPOINTS = [
    ("GET", "/v1/rooms"),
    ("GET", "/v1/rooms/any/geometries"),
    ("GET", "/v1/rooms/any/geometries/key/selection"),
    ("GET", "/v1/rooms/any/geometries/key"),
    ("GET", "/v1/rooms/any/figures"),
    ("GET", "/v1/rooms/any/figures/k"),
    ("GET", "/v1/rooms/any/trajectory"),
]


@pytest.mark.asyncio
@pytest.mark.parametrize("method,path", ENDPOINTS)
async def test_endpoint_requires_auth(
    client: AsyncClient, method: str, path: str
) -> None:
    r = await client.request(method, path)
    assert r.status_code == 401, f"{method} {path} returned {r.status_code}"
