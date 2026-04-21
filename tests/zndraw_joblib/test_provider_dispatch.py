"""Regression: a failing taskiq dispatch must release the inflight lock."""

from __future__ import annotations

import asyncio
import contextlib
import uuid


async def _seed_internal_provider(async_session_factory):
    """Create an @internal provider in the DB and return its full_name."""
    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord, Worker

    full_name = "@internal:filesystem:FilesystemRead"
    async with async_session_factory() as session:
        user = User(
            id=uuid.uuid4(),
            email="int-dispatch@test",
            hashed_password="x",
            is_active=True,
            is_superuser=True,
            is_verified=True,
        )
        session.add(user)
        worker = Worker(user_id=user.id)
        session.add(worker)
        await session.flush()
        session.add(
            ProviderRecord(
                room_id="@internal",
                category="filesystem",
                name="FilesystemRead",
                schema_={},
                user_id=user.id,
                worker_id=worker.id,
            )
        )
        await session.commit()
    return full_name


def test_read_provider_releases_inflight_on_dispatch_failure(
    client, app, async_session_factory
):
    """Inflight lock must be released when dispatch raises.

    Pre-fix: if kiq() raises, the lock was never released.  The second
    request finds `acquired=False`, skips dispatch, and waits for a result
    that never comes — wedging for ``provider_inflight_ttl_seconds`` (30s).

    Post-fix: the try/except wrapping ``acquired`` branch releases the lock
    on any exception, so the second request can attempt dispatch again.
    """
    from zndraw_joblib.dependencies import get_result_backend, request_hash
    from zndraw_joblib.registry import InternalProviderRegistry

    asyncio.run(_seed_internal_provider(async_session_factory))

    class _FailingTask:
        """kiq() always raises to simulate broker/transport failure."""

        async def kiq(self, **kwargs):  # noqa: ARG002 — mock signature mirrors taskiq
            raise RuntimeError("simulated dispatch failure")

    registry = InternalProviderRegistry(
        tasks={"@internal:filesystem:FilesystemRead": _FailingTask()},
        providers={},
    )
    app.state.internal_provider_registry = registry

    # Retrieve the InMemoryResultBackend that _build_app installed
    result_backend = app.dependency_overrides[get_result_backend]()

    # First request — dispatch raises; starlette TestClient re-raises server
    # exceptions by default.  We catch it here to verify the inflight lock is
    # still released (the try/except+raise in the handler runs before the error
    # propagates up through the ASGI stack).
    with contextlib.suppress(RuntimeError):
        # Dispatch failure propagates through TestClient — that's expected.
        client.get(
            "/v1/joblib/rooms/@internal/providers/@internal:filesystem:FilesystemRead"
            "?path=/data",
            headers={"Prefer": "wait=0"},
        )

    # Inflight lock MUST be released after the dispatch failure
    params = {"path": "/data"}
    rhash = request_hash(params)
    inflight_key = f"provider-inflight:@internal:filesystem:FilesystemRead:{rhash}"
    assert inflight_key not in result_backend._inflight, (
        "Inflight lock was not released after dispatch failure — "
        "second request would wedge for inflight_ttl seconds"
    )

    # Second request — must also reach dispatch (not skip due to stuck lock).
    # We verify by counting kiq calls: if the lock was NOT released, the second
    # request would find acquired=False, skip dispatch entirely, and wait for a
    # result that never comes (504 timeout or just returning without a kiq call).
    kiq_call_count = 0

    class _CountingTask:
        async def kiq(self, **kwargs):  # noqa: ARG002 — mock signature
            nonlocal kiq_call_count
            kiq_call_count += 1
            raise RuntimeError("simulated dispatch failure")

    app.state.internal_provider_registry = InternalProviderRegistry(
        tasks={"@internal:filesystem:FilesystemRead": _CountingTask()},
        providers={},
    )

    with contextlib.suppress(RuntimeError):
        client.get(
            "/v1/joblib/rooms/@internal/providers/@internal:filesystem:FilesystemRead"
            "?path=/data",
            headers={"Prefer": "wait=0"},
        )

    assert kiq_call_count == 1, (
        f"Expected second request to dispatch (kiq_call_count=1), "
        f"got kiq_call_count={kiq_call_count} — inflight lock was NOT released"
    )
