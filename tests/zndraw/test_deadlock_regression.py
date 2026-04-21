"""Regression tests for the @internal provider read deadlock (PR #920).

Pre-fix: ``get_worker_token`` used yield-based session DI; combined with
``SessionMakerDep`` in ``read_provider`` it acquired the SQLite serialization
lock twice in one request, deadlocking.  These tests fail if the lock-holding
pattern regresses.
"""

from __future__ import annotations

import asyncio

import httpx
import pytest


@pytest.mark.asyncio
async def test_read_provider_completes_without_deadlock(server_factory):
    server = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": ".",
        }
    )
    async with httpx.AsyncClient(base_url=server.url, timeout=5.0) as client:
        # @internal filesystem provider is seeded by lifespan.
        resp = await client.get(
            "/v1/joblib/rooms/@internal/providers/"
            "@internal:filesystem:FilesystemRead?path=/",
        )
    assert resp.status_code in (200, 401, 403), (
        f"Expected fast response, got {resp.status_code} (likely a deadlock)"
    )


@pytest.mark.asyncio
async def test_concurrent_reads_no_deadlock(server_factory):
    server = server_factory(
        {
            "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
            "ZNDRAW_SERVER_FILEBROWSER_PATH": ".",
        }
    )
    async with httpx.AsyncClient(base_url=server.url, timeout=5.0) as client:
        tasks = [
            client.get(
                "/v1/joblib/rooms/@internal/providers/"
                "@internal:filesystem:FilesystemRead?path=/",
            )
            for _ in range(3)
        ]
        responses = await asyncio.gather(*tasks, return_exceptions=True)
    for r in responses:
        assert not isinstance(r, Exception), r
