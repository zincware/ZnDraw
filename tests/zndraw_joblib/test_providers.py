# tests/test_providers.py
"""Tests for provider system endpoints."""

import asyncio
import json
import uuid

from zndraw_auth import User
from zndraw_joblib.dependencies import request_hash
from zndraw_joblib.models import ProviderRecord
from zndraw_joblib.schemas import (
    PaginatedResponse,
    ProviderResponse,
)


def _register_provider(client, room_id="@global", **overrides):
    """Helper to register a provider with sensible defaults."""
    payload = {
        "category": "filesystem",
        "name": "local",
        "schema": {"path": {"type": "string", "default": "/"}},
    }
    payload.update(overrides)
    return client.put(f"/v1/joblib/rooms/{room_id}/providers", json=payload)


# --- Registration ---


def test_register_provider(client):
    resp = _register_provider(client)
    assert resp.status_code == 201
    data = ProviderResponse.model_validate(resp.json())
    assert data.category == "filesystem"
    assert data.name == "local"
    assert data.room_id == "@global"
    assert data.full_name == "@global:filesystem:local"
    assert data.worker_id is not None


def test_register_provider_auto_creates_worker(client):
    resp = _register_provider(client)
    assert resp.status_code == 201
    data = resp.json()
    assert data["worker_id"] is not None
    # Worker should exist
    worker_id = data["worker_id"]
    worker_resp = client.get("/v1/joblib/workers")
    assert worker_resp.status_code == 200
    worker_ids = [w["id"] for w in worker_resp.json()["items"]]
    assert worker_id in worker_ids


def test_register_provider_with_existing_worker(client):
    # Create worker first
    worker_resp = client.post("/v1/joblib/workers")
    worker_id = worker_resp.json()["id"]

    resp = _register_provider(client, worker_id=worker_id)
    assert resp.status_code == 201
    data = resp.json()
    assert data["worker_id"] == worker_id


def test_register_provider_idempotent(client):
    resp1 = _register_provider(client)
    assert resp1.status_code == 201
    id1 = resp1.json()["id"]

    resp2 = _register_provider(client)
    assert resp2.status_code == 200
    id2 = resp2.json()["id"]
    assert id1 == id2


def test_register_provider_different_names(client):
    resp1 = _register_provider(client, name="local")
    assert resp1.status_code == 201
    resp2 = _register_provider(client, name="s3-bucket")
    assert resp2.status_code == 201
    assert resp1.json()["id"] != resp2.json()["id"]


def test_register_provider_different_rooms(client):
    resp1 = _register_provider(client, room_id="@global")
    assert resp1.status_code == 201
    resp2 = _register_provider(client, room_id="room-42")
    assert resp2.status_code == 201
    assert resp1.json()["id"] != resp2.json()["id"]


# --- Listing ---


def test_list_providers_empty(client):
    resp = client.get("/v1/joblib/rooms/@global/providers")
    assert resp.status_code == 200
    data = PaginatedResponse[ProviderResponse].model_validate(resp.json())
    assert data.total == 0
    assert data.items == []


def test_list_providers_global(client):
    _register_provider(client)
    resp = client.get("/v1/joblib/rooms/@global/providers")
    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 1
    assert data["items"][0]["category"] == "filesystem"


def test_list_providers_room_includes_global(client):
    _register_provider(client, room_id="@global")
    resp = client.get("/v1/joblib/rooms/room-42/providers")
    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 1


def test_list_providers_room_scoped(client):
    _register_provider(client, room_id="room-42")
    resp = client.get("/v1/joblib/rooms/room-42/providers")
    assert resp.status_code == 200
    assert resp.json()["total"] == 1

    # Not visible from other rooms
    resp2 = client.get("/v1/joblib/rooms/room-99/providers")
    assert resp2.status_code == 200
    assert resp2.json()["total"] == 0


def test_list_providers_mixed_scopes(client):
    _register_provider(client, room_id="@global", name="global-fs")
    _register_provider(client, room_id="room-42", name="room-fs")
    _register_provider(client, room_id="room-99", name="other-fs")

    # room-42 sees global + room-42
    resp = client.get("/v1/joblib/rooms/room-42/providers")
    assert resp.json()["total"] == 2

    # @global only sees global
    resp = client.get("/v1/joblib/rooms/@global/providers")
    assert resp.json()["total"] == 1


def test_list_providers_includes_internal(client, async_session_factory):
    """Internal providers are visible from every room (and from @global)."""

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.flush()
            session.add(
                ProviderRecord(
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemRead",
                    schema_={},
                    user_id=user.id,
                    worker_id=None,
                )
            )
            await session.commit()

    asyncio.run(seed())

    # A normal room sees the @internal provider
    resp = client.get("/v1/joblib/rooms/room-42/providers")
    assert resp.status_code == 200
    items = resp.json()["items"]
    assert any(p["full_name"] == "@internal:filesystem:FilesystemRead" for p in items)

    # @global does NOT see @internal (explicit scope isolation)
    resp = client.get("/v1/joblib/rooms/@global/providers")
    assert all(
        p["full_name"] != "@internal:filesystem:FilesystemRead"
        for p in resp.json()["items"]
    )


def test_get_provider_info_internal_visible_from_room(client, async_session_factory):
    """A normal room can fetch info on an @internal provider."""

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.flush()
            session.add(
                ProviderRecord(
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemRead",
                    schema_={"path": {"type": "string"}},
                    user_id=user.id,
                    worker_id=None,
                )
            )
            await session.commit()

    asyncio.run(seed())

    resp = client.get(
        "/v1/joblib/rooms/room-42/providers/@internal:filesystem:FilesystemRead/info"
    )
    assert resp.status_code == 200
    assert resp.json()["schema"] == {"path": {"type": "string"}}


# --- Info ---


def test_get_provider_info(client):
    _register_provider(client)
    resp = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local/info"
    )
    assert resp.status_code == 200
    data = ProviderResponse.model_validate(resp.json())
    assert data.category == "filesystem"
    assert data.name == "local"
    assert data.schema_ == {"path": {"type": "string", "default": "/"}}


def test_get_provider_info_not_found(client):
    resp = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:nonexistent/info"
    )
    assert resp.status_code == 404


def test_get_provider_info_room_visibility(client):
    _register_provider(client, room_id="room-42")
    # Visible from room-42
    resp = client.get(
        "/v1/joblib/rooms/room-42/providers/room-42:filesystem:local/info"
    )
    assert resp.status_code == 200
    # NOT visible from room-99
    resp = client.get(
        "/v1/joblib/rooms/room-99/providers/room-42:filesystem:local/info"
    )
    assert resp.status_code == 404


# --- Data Read (long-poll / 200 / 504) ---

IMMEDIATE = {"Prefer": "wait=0"}


def test_read_provider_timeout(client):
    """No cached result + immediate timeout → 504 with Retry-After."""
    _register_provider(client)
    resp = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local?path=/data",
        headers=IMMEDIATE,
    )
    assert resp.status_code == 504
    assert resp.headers["content-type"] == "application/problem+json"
    assert "Retry-After" in resp.headers


def test_read_provider_not_found(client):
    resp = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:nonexistent?path=/data",
        headers=IMMEDIATE,
    )
    assert resp.status_code == 404


def test_read_provider_cached_200(client):
    """Upload result first, then read returns 200 immediately."""
    reg_resp = _register_provider(client)
    provider_id = reg_resp.json()["id"]

    params = {"path": "/data"}
    rhash = request_hash(params)

    # Upload result
    upload_resp = client.post(
        f"/v1/joblib/providers/{provider_id}/results",
        content=json.dumps([{"name": "file.xyz", "size": 42}]).encode(),
        headers={"X-Request-Hash": rhash},
    )
    assert upload_resp.status_code == 204

    # Read → 200 with cached data
    resp = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local", params=params
    )
    assert resp.status_code == 200
    assert resp.json() == [{"name": "file.xyz", "size": 42}]


def test_read_provider_inflight_coalescing(client):
    """Second concurrent read should not re-dispatch."""
    _register_provider(client)
    params = {"path": "/data"}

    # First read acquires inflight, times out
    resp1 = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local",
        params=params,
        headers=IMMEDIATE,
    )
    assert resp1.status_code == 504

    # Second read also times out (inflight already acquired, no re-dispatch)
    resp2 = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local",
        params=params,
        headers=IMMEDIATE,
    )
    assert resp2.status_code == 504


def test_read_provider_cache_hit_skips_long_poll(client):
    """Cache hit returns immediately without Preference-Applied header."""
    reg_resp = _register_provider(client)
    provider_id = reg_resp.json()["id"]

    params = {"path": "/data"}
    rhash = request_hash(params)

    # Upload result first
    client.post(
        f"/v1/joblib/providers/{provider_id}/results",
        content=b'"ok"',
        headers={"X-Request-Hash": rhash},
    )

    # Cache hit returns 200 without Preference-Applied (no wait happened)
    resp = client.get(
        "/v1/joblib/rooms/@global/providers/@global:filesystem:local",
        params=params,
        headers={"Prefer": "wait=5"},
    )
    assert resp.status_code == 200
    assert "Preference-Applied" not in resp.headers


async def test_read_provider_long_poll_wakes_on_upload(async_client):
    """Result uploaded during active long-poll wakes the waiter and returns 200."""
    # Register provider
    reg_resp = await async_client.put(
        "/v1/joblib/rooms/@global/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert reg_resp.status_code == 201
    provider_id = reg_resp.json()["id"]

    params = {"path": "/data"}
    rhash = request_hash(params)
    payload = json.dumps({"files": ["a.xyz"]}).encode()

    async def upload_after_delay():
        await asyncio.sleep(0.1)
        resp = await async_client.post(
            f"/v1/joblib/providers/{provider_id}/results",
            content=payload,
            headers={"X-Request-Hash": rhash},
        )
        assert resp.status_code == 204

    # Long-poll with 5s timeout; upload arrives after 100ms
    read_task = asyncio.create_task(
        async_client.get(
            "/v1/joblib/rooms/@global/providers/@global:filesystem:local",
            params=params,
            headers={"Prefer": "wait=5"},
        )
    )
    upload_task = asyncio.create_task(upload_after_delay())

    read_resp, _ = await asyncio.gather(read_task, upload_task)
    assert read_resp.status_code == 200
    assert read_resp.json() == {"files": ["a.xyz"]}
    assert "Preference-Applied" in read_resp.headers


# --- Delete ---


def test_delete_provider(client):
    reg_resp = _register_provider(client)
    provider_id = reg_resp.json()["id"]

    resp = client.delete(f"/v1/joblib/providers/{provider_id}")
    assert resp.status_code == 204

    # No longer listed
    list_resp = client.get("/v1/joblib/rooms/@global/providers")
    assert list_resp.json()["total"] == 0


def test_delete_provider_not_found(client):
    resp = client.delete(f"/v1/joblib/providers/{uuid.uuid4()}")
    assert resp.status_code == 404


def test_delete_internal_provider_forbidden(client_factory, async_session_factory):
    """@internal providers cannot be deleted by anyone, including superusers."""

    provider_id = uuid.uuid4()

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.flush()
            session.add(
                ProviderRecord(
                    id=provider_id,
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemRead",
                    schema_={},
                    user_id=user.id,
                    worker_id=None,
                )
            )
            await session.commit()

    asyncio.run(seed())

    admin = client_factory("admin", is_superuser=True)
    resp = admin.delete(f"/v1/joblib/providers/{provider_id}")
    assert resp.status_code == 403


# --- Result Upload ---


def test_upload_result_not_found(client):
    resp = client.post(
        f"/v1/joblib/providers/{uuid.uuid4()}/results",
        content=b"{}",
        headers={"X-Request-Hash": "abc"},
    )
    assert resp.status_code == 404


# --- Worker Cascade ---


def test_worker_delete_cascades_providers(client):
    """Deleting a worker should also remove its providers."""
    reg_resp = _register_provider(client)
    worker_id = reg_resp.json()["worker_id"]

    # Provider exists
    list_resp = client.get("/v1/joblib/rooms/@global/providers")
    assert list_resp.json()["total"] == 1

    # Delete worker
    del_resp = client.delete(f"/v1/joblib/workers/{worker_id}")
    assert del_resp.status_code == 204

    # Provider gone
    list_resp = client.get("/v1/joblib/rooms/@global/providers")
    assert list_resp.json()["total"] == 0


# --- Request Hash ---


def test_request_hash_deterministic():
    params = {"path": "/data", "glob": "*.xyz"}
    h1 = request_hash(params)
    h2 = request_hash(params)
    assert h1 == h2


def test_request_hash_order_independent():
    h1 = request_hash({"a": 1, "b": 2})
    h2 = request_hash({"b": 2, "a": 1})
    assert h1 == h2


def test_request_hash_different_params():
    h1 = request_hash({"path": "/data"})
    h2 = request_hash({"path": "/other"})
    assert h1 != h2


# --- Authorization ---


def test_delete_provider_forbidden_other_user(client_factory):
    """A non-superuser cannot delete another user's provider."""
    alice = client_factory("alice", is_superuser=False)
    bob = client_factory("bob", is_superuser=False)

    # Non-superusers can't register @global providers, so use a room.
    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201
    provider_id = resp.json()["id"]

    # Bob cannot delete Alice's provider
    del_resp = bob.delete(f"/v1/joblib/providers/{provider_id}")
    assert del_resp.status_code == 403


def test_upload_result_forbidden_other_user(client_factory):
    """A non-superuser cannot upload results for another user's provider."""
    alice = client_factory("alice", is_superuser=False)
    bob = client_factory("bob", is_superuser=False)

    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201
    provider_id = resp.json()["id"]

    upload_resp = bob.post(
        f"/v1/joblib/providers/{provider_id}/results",
        content=b"{}",
        headers={"X-Request-Hash": "abc"},
    )
    assert upload_resp.status_code == 403


def test_register_provider_forbidden_other_user(client_factory):
    """User B cannot overwrite user A's provider registration."""
    alice = client_factory("alice", is_superuser=False)
    bob = client_factory("bob", is_superuser=False)

    # Alice registers a provider in a room
    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201

    # Bob tries to register the same provider name -> should be rejected
    resp = bob.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 403


def test_register_provider_superuser_can_overwrite(client_factory):
    """Superuser can overwrite another user's provider registration."""
    alice = client_factory("alice", is_superuser=False)
    admin = client_factory("admin", is_superuser=True)

    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201

    # Admin can overwrite
    resp = admin.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 200


def test_register_global_provider_requires_superuser(client_factory):
    """Non-superusers cannot register @global providers."""
    normal = client_factory("normal-user", is_superuser=False)
    resp = normal.put(
        "/v1/joblib/rooms/@global/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 403


def test_register_global_provider_superuser_ok(client_factory):
    """Superusers can register @global providers."""
    admin = client_factory("admin-user", is_superuser=True)
    resp = admin.put(
        "/v1/joblib/rooms/@global/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201


def test_register_internal_provider_forbidden_normal_user(client_factory):
    """@internal is a reserved server-bootstrap namespace — HTTP registration is
    forbidden for normal users. Mirrors the DELETE-side policy at router.py:1337."""
    normal = client_factory("normal-user", is_superuser=False)
    resp = normal.put(
        "/v1/joblib/rooms/@internal/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 403, resp.text


def test_register_internal_provider_forbidden_superuser(client_factory):
    """@internal cannot be registered via HTTP even by a superuser — the namespace
    is seeded only by ``register_internal_providers`` at startup."""
    admin = client_factory("admin-user", is_superuser=True)
    resp = admin.put(
        "/v1/joblib/rooms/@internal/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 403, resp.text


def test_register_provider_category_rejected(app, client):
    """Provider category rejected when allowed_provider_categories is set."""
    app.state.joblib_settings.allowed_provider_categories = ["frames"]
    resp = _register_provider(client, category="filesystem")
    assert resp.status_code == 400
    app.state.joblib_settings.allowed_provider_categories = None  # reset


def test_register_provider_category_allowed(app, client):
    """Provider category accepted when in the allowed list."""
    app.state.joblib_settings.allowed_provider_categories = ["filesystem"]
    resp = _register_provider(client, category="filesystem")
    assert resp.status_code == 201
    app.state.joblib_settings.allowed_provider_categories = None  # reset


def test_read_internal_provider_dispatches_via_taskiq(
    client, app, async_session_factory
):
    """An @internal provider dispatches via the registry's taskiq task, not tsio."""
    import asyncio
    import uuid

    from zndraw_joblib.registry import InternalProviderRegistry

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.flush()
            session.add(
                ProviderRecord(
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemRead",
                    schema_={},
                    user_id=user.id,
                    worker_id=None,
                )
            )
            await session.commit()

    asyncio.run(seed())

    kiq_calls: list[dict] = []

    class _FakeTask:
        async def kiq(self, **kwargs):
            kiq_calls.append(kwargs)

    registry = InternalProviderRegistry(
        tasks={"@internal:filesystem:FilesystemRead": _FakeTask()},
        providers={},
    )
    app.state.internal_provider_registry = registry

    # Immediate timeout so we don't hang — we only care that kiq was called
    resp = client.get(
        "/v1/joblib/rooms/room-42/providers/@internal:filesystem:FilesystemRead"
        "?path=/data",
        headers={"Prefer": "wait=0"},
    )
    # No result uploaded => 504, but dispatch must have happened
    assert resp.status_code == 504
    assert len(kiq_calls) == 1
    call = kiq_calls[0]
    assert json.loads(call["params_json"]) == {"path": "/data"}
    assert "request_id" in call
    assert "provider_id" in call
    # Token is now minted on-demand inside the @internal branch via
    # mint_internal_worker_token — it is a real JWT, not the stub string.
    assert isinstance(call["token"], str)
    assert len(call["token"]) > 0


def test_provider_response_from_record_accepts_null_worker_id():
    """@internal providers have worker_id=None; from_record must not raise,
    and the JSON payload must emit ``"worker_id": null`` (not the string
    "null", not an omitted key) so frontend clients see the field as null.
    """
    import json
    import uuid
    from datetime import UTC, datetime

    from zndraw_joblib.schemas import ProviderResponse

    record = ProviderRecord(
        id=uuid.uuid4(),
        room_id="@internal",
        category="filesystem",
        name="FilesystemRead",
        schema_={},
        content_type="application/json",
        user_id=uuid.uuid4(),
        worker_id=None,
        created_at=datetime.now(UTC),
    )
    response = ProviderResponse.from_record(record)
    assert response.worker_id is None
    assert response.room_id == "@internal"
    # JSON round-trip: field is present with JSON null.
    payload = json.loads(response.model_dump_json())
    assert "worker_id" in payload
    assert payload["worker_id"] is None


def test_global_scope_cannot_resolve_room_provider(client_factory):
    """A @global caller must not reach room-scoped providers (security: LIST
    policy excludes them; resolve must agree)."""
    alice = client_factory("alice-b3", is_superuser=False)
    admin = client_factory("admin-b3", is_superuser=True)

    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201

    # admin, calling with @global scope, must not resolve a room-42 provider.
    resp = admin.get("/v1/joblib/rooms/@global/providers/room-42:filesystem:local")
    assert resp.status_code == 404, resp.text


def test_global_scope_cannot_resolve_internal_provider(
    client_factory, async_session_factory
):
    """A @global caller must not reach @internal providers either — test
    pins the LIST-vs-RESOLVE symmetry documented in
    test_list_providers_includes_internal."""
    import asyncio
    import uuid

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int-b3@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.flush()
            session.add(
                ProviderRecord(
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemReadB3",
                    schema_={},
                    user_id=user.id,
                    worker_id=None,
                )
            )
            await session.commit()

    asyncio.run(seed())

    admin = client_factory("admin-b3-2", is_superuser=True)
    resp = admin.get(
        "/v1/joblib/rooms/@global/providers/@internal:filesystem:FilesystemReadB3"
    )
    assert resp.status_code == 404, resp.text


# --- B5 regression: worker_token must not be resolved for remote providers ---


def test_read_remote_provider_works_with_no_internal_worker_cache(
    unguarded_client_factory,
):
    """A remote provider read must not require app.state.internal_worker_user.

    That cache is only needed for @internal dispatch. On a fresh deploy where
    init_db_on_startup=False and the worker-seed step has not run yet, every
    provider read would 500 if WorkerTokenDep is a route-level dependency —
    including remote reads that never need the token.
    """
    # Build a client where get_worker_token is NOT stubbed and
    # internal_worker_user is None (simulating a cold deploy).
    alice = unguarded_client_factory(
        "alice-b5", is_superuser=False, internal_worker_user=None
    )

    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201, resp.text

    # A remote provider read should reach the dispatch path without touching
    # internal_worker_user. Acceptable outcomes:
    #   504 — no remote worker delivered results within the wait window
    #   409 — no connected worker (NoWorkersAvailable)
    # Unacceptable: 500 (WorkerTokenDep raises RuntimeError before dispatch).
    resp = alice.get(
        "/v1/joblib/rooms/room-42/providers/room-42:filesystem:local?path=/",
        headers={"Prefer": "wait=0"},
    )
    assert resp.status_code in (504, 409), (
        f"Expected 504 or 409 (dispatch reached), got {resp.status_code}: {resp.text}"
    )


def test_legitimate_json_with_error_type_keys_is_not_mis_flagged(
    client_factory,
):
    """A provider returning valid JSON with top-level 'error' and 'type'
    keys (e.g., JSON-Schema) must be returned as-is, not translated to
    an HTTP 400."""
    alice = client_factory("alice-b6a", is_superuser=True)
    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201
    data = resp.json()
    provider_id = data["id"]
    provider_full_name = data["full_name"]

    payload = json.dumps({"type": "object", "error": None, "ok": True}).encode()
    rhash = request_hash({"path": "/"})

    # Seed the result via the upload endpoint (no X-Result-Status = success)
    upload_resp = alice.post(
        f"/v1/joblib/providers/{provider_id}/results",
        content=payload,
        headers={"X-Request-Hash": rhash},
    )
    assert upload_resp.status_code == 204

    resp = alice.get(f"/v1/joblib/rooms/room-42/providers/{provider_full_name}?path=/")
    assert resp.status_code == 200, resp.text
    assert resp.json() == {"type": "object", "error": None, "ok": True}


def test_provider_error_path_returns_problem_detail(
    client_factory,
):
    """When an executor posts an error, read_provider returns RFC 9457
    problem+json with the status from the payload."""
    from zndraw_joblib.exceptions import ProviderExecutionFailed

    alice = client_factory("alice-b6b", is_superuser=True)
    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201
    data = resp.json()
    provider_id = data["id"]
    provider_full_name = data["full_name"]

    problem = ProviderExecutionFailed.create(detail="FileNotFoundError: /nope")
    payload = problem.model_dump_json(exclude_none=True).encode()
    rhash = request_hash({"path": "/nope"})

    # Seed the result via the upload endpoint with X-Result-Status: error
    upload_resp = alice.post(
        f"/v1/joblib/providers/{provider_id}/results",
        content=payload,
        headers={
            "X-Request-Hash": rhash,
            "X-Result-Status": "error",
            "Content-Type": "application/problem+json",
        },
    )
    assert upload_resp.status_code == 204

    resp = alice.get(
        f"/v1/joblib/rooms/room-42/providers/{provider_full_name}?path=/nope"
    )
    assert resp.status_code == 400
    assert resp.headers["content-type"].startswith("application/problem+json")
    body = resp.json()
    assert body["title"] == "Bad Request"
    assert body["status"] == 400
    assert "FileNotFoundError" in body["detail"]


def test_error_status_visible_before_payload_via_notify(client_factory):
    """After upload_result with X-Result-Status: error, the status key
    must exist at every moment cache_key exists — no window where a
    reader can see payload without status."""
    from zndraw_joblib.dependencies import get_result_backend

    alice = client_factory("alice-b6c", is_superuser=True)
    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201
    provider_id = resp.json()["id"]
    provider_full_name = resp.json()["full_name"]

    params = {"path": "/test"}
    rhash = request_hash(params)
    cache_key = f"provider-result:{provider_full_name}:{rhash}"
    status_key = f"{cache_key}:status"

    # Post an error payload through the real upload endpoint.
    error_body = (
        b'{"type":"/v1/problems/provider-execution-failed",'
        b'"title":"Bad Request","status":400,"detail":"X"}'
    )
    upload_resp = alice.post(
        f"/v1/joblib/providers/{provider_id}/results",
        content=error_body,
        headers={
            "X-Request-Hash": rhash,
            "X-Result-Status": "error",
            "Content-Type": "application/problem+json",
        },
    )
    assert upload_resp.status_code == 204

    # Both keys must be present now.
    # Get the result_backend from the app's dependency overrides
    result_backend = alice.app.dependency_overrides[get_result_backend]()

    async def _check() -> None:
        assert await result_backend.get(status_key) == b"error"
        assert await result_backend.get(cache_key) is not None

    asyncio.run(_check())

    # And read_provider must see the error branch, not success.
    resp = alice.get(
        f"/v1/joblib/rooms/room-42/providers/{provider_full_name}",
        params=params,
    )
    assert resp.status_code == 400
    assert resp.headers["content-type"].startswith("application/problem+json")


# --- B7: filebrowser_require_superuser gate ---


def test_internal_filesystem_requires_superuser_by_default(
    client_factory, async_session_factory
):
    """With the default filebrowser_require_superuser=True, a non-superuser
    is 403'd on @internal:filesystem:* read and the provider is absent from LIST.
    """
    import asyncio
    import uuid

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int-b7a@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.flush()
            session.add(
                ProviderRecord(
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemReadB7A",
                    schema_={},
                    user_id=user.id,
                    worker_id=None,
                )
            )
            await session.commit()

    asyncio.run(seed())

    alice = client_factory("alice-b7a", is_superuser=False)

    # Read must return 403
    resp = alice.get(
        "/v1/joblib/rooms/room-42/providers/@internal:filesystem:FilesystemReadB7A?path=/"
    )
    assert resp.status_code == 403, resp.text

    # List must not include the gated provider
    resp = alice.get("/v1/joblib/rooms/room-42/providers")
    assert resp.status_code == 200
    items = resp.json()["items"]
    assert not any(
        p["full_name"].startswith("@internal:filesystem:FilesystemReadB7A")
        for p in items
    )


def test_internal_filesystem_superuser_bypasses_gate(
    client_factory, async_session_factory
):
    """Superusers see and can read @internal:filesystem:* even with the
    default flag on.
    """
    import asyncio
    import uuid

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int-b7b@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.flush()
            session.add(
                ProviderRecord(
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemReadB7B",
                    schema_={},
                    user_id=user.id,
                    worker_id=None,
                )
            )
            await session.commit()

    asyncio.run(seed())

    admin = client_factory("admin-b7b", is_superuser=True)
    resp = admin.get("/v1/joblib/rooms/room-42/providers")
    items = resp.json()["items"]
    assert any(
        p["full_name"] == "@internal:filesystem:FilesystemReadB7B" for p in items
    )


def test_internal_filesystem_gate_disabled_by_flag(
    client_factory, async_session_factory
):
    """With filebrowser_require_superuser=False, non-superusers can access."""
    import asyncio
    import uuid

    alice = client_factory("alice-b7c", is_superuser=False)
    alice.app.state.joblib_settings.filebrowser_require_superuser = False

    try:

        async def seed() -> None:
            async with async_session_factory() as session:
                user = User(
                    id=uuid.uuid4(),
                    email="int-b7c@test",
                    hashed_password="x",
                    is_active=True,
                    is_superuser=True,
                    is_verified=True,
                )
                session.add(user)
                await session.flush()
                session.add(
                    ProviderRecord(
                        room_id="@internal",
                        category="filesystem",
                        name="FilesystemReadB7C",
                        schema_={},
                        user_id=user.id,
                        worker_id=None,
                    )
                )
                await session.commit()

        asyncio.run(seed())

        resp = alice.get("/v1/joblib/rooms/room-42/providers")
        items = resp.json()["items"]
        assert any(
            p["full_name"] == "@internal:filesystem:FilesystemReadB7C" for p in items
        )
    finally:
        alice.app.state.joblib_settings.filebrowser_require_superuser = True


def test_list_providers_pagination_correct_with_gate(
    client_factory, async_session_factory
):
    """When the gate is on and the DB has many @internal:filesystem:*
    rows, pagination.total must reflect the post-filter count, not the
    current page's len(items). Regression for B7's SQL-level filter."""
    import asyncio
    import uuid

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int-b7d@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
            session.add(user)
            await session.flush()
            # Add a visible room-42 provider too — ensure filtering
            # doesn't leak into non-gated rows.
            session.add(
                ProviderRecord(
                    room_id="room-42",
                    category="filesystem",
                    name="VisibleRoomProvider",
                    schema_={},
                    user_id=user.id,
                    worker_id=None,
                )
            )
            # Add 5 gated @internal:filesystem:* rows
            for i in range(5):
                session.add(
                    ProviderRecord(
                        room_id="@internal",
                        category="filesystem",
                        name=f"FilesystemReadB7D{i}",
                        schema_={},
                        user_id=user.id,
                        worker_id=None,
                    )
                )
            await session.commit()

    asyncio.run(seed())

    alice = client_factory("alice-b7d", is_superuser=False)
    resp = alice.get("/v1/joblib/rooms/room-42/providers?limit=100")
    assert resp.status_code == 200
    data = resp.json()
    # None of the 5 gated providers should appear
    names = [p["full_name"] for p in data["items"]]
    assert not any(
        n.startswith("@internal:filesystem:FilesystemReadB7D") for n in names
    )
    # And total must be len(items) — if the old post-filter logic set
    # total correctly we'd see it match. With SQL-level filter they
    # agree by construction.
    assert data["total"] == len(data["items"])
    # The room-42 provider should be visible
    assert any(n == "room-42:filesystem:VisibleRoomProvider" for n in names)
