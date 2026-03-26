# tests/test_providers.py
"""Tests for provider system endpoints."""

import asyncio
import json
import uuid

from zndraw_joblib.dependencies import request_hash
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

    # Alice registers a provider
    resp = alice.put(
        "/v1/joblib/rooms/@global/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    # Non-superuser can't register @global, so use a room
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
