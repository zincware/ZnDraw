"""Every new problem type must appear in at least one route's OpenAPI responses."""

import json

import pytest
from httpx import AsyncClient

NEW_TYPES = {
    "group-not-found",
    "group-name-taken",
    "not-group-member",
    "not-group-admin",
    "last-group-admin",
    "group-has-rooms",
    "transfer-target-invalid",
    "share-link-not-found",
    "share-link-invalid",
}


@pytest.mark.asyncio
async def test_openapi_references_all_new_problem_types(client: AsyncClient) -> None:
    spec = (await client.get("/openapi.json")).json()
    blob = json.dumps(spec)
    missing = {t for t in NEW_TYPES if f"/v1/problems/{t}" not in blob}
    assert not missing, f"Missing from OpenAPI: {missing}"
