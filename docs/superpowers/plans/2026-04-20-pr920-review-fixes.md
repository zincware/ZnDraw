# PR #920 Code-Review Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Address the 20 review findings documented in `docs/superpowers/specs/2026-04-20-pr920-review-fixes-design.md` across backend (`src/zndraw`, `src/zndraw_joblib`) and frontend (`frontend/src`).

**Architecture:** Each fix is one commit. TDD where a regression test is possible (most backend items); test-after for pure refactors (helper extraction, hook consolidation). Trivial cleanups (file deletion, dep removal) get a single commit without a test. Backend and frontend are independent — no ordering constraint across the split.

**Tech Stack:**
- Backend: Python 3.12 / pytest / FastAPI / SQLModel / Pydantic v2 / httpx / taskiq.
- Frontend: TypeScript / React 19 / vitest / Bun / Zustand / TanStack Query / dockview-react / MUI / Biome.

---

## File Structure

### Backend — created
- *(none — all edits hit existing files)*

### Backend — modified
- `src/zndraw_joblib/exceptions.py` — new `ProviderExecutionFailed` ProblemType
- `src/zndraw_joblib/router.py` — lazy worker_token, visibility fix, superuser gate, status-header preservation
- `src/zndraw_joblib/dependencies.py` — extract token-minting helper callable outside a FastAPI handler
- `src/zndraw/providers/executor.py` — ProblemDetail error payload + `log.exception`
- `src/zndraw/config.py` — `filebrowser_enabled` flag + `filebrowser_require_superuser` flag, drop `env_parse_none_str`
- `src/zndraw/result_backends.py` — `key_prefix` on `StorageResultBackend`
- `src/zndraw/extensions/modifiers.py` — lift `LoadFile` import to top
- `src/zndraw/database.py` — pass `key_prefix` into `StorageResultBackend`, honor `filebrowser_enabled`
- `tests/zndraw_joblib/test_registry.py` — add `token` param to `mock_executor`
- `tests/zndraw_joblib/test_providers.py` — new regression tests
- `tests/zndraw/test_providers_filesystem.py` — delete weak test, rename remaining, use new env vars
- `tests/zndraw/test_result_backends.py` — regression test for storage key_prefix

### Frontend — created
- `frontend/src/stores/dockviewApiStore.ts` — Zustand store holding the DockviewApi
- `frontend/src/panels/dragStyles.ts` — shared `shimmer` keyframes
- `frontend/src/panels/useDragHover.ts` — shared drag-hover hook
- `frontend/src/hooks/useFilesystemProviders.ts` — shared query+socket hook
- `frontend/src/stores/slices/__tests__/activityBarSlice.test.ts` — slice unit tests

### Frontend — modified
- `frontend/src/panels/DockviewLayout.tsx` — remove `sharedApi`, add `resetDockview`/`ensureViewerPanel`
- `frontend/src/panels/PlotsBrowserPanel.tsx` — drop polling, subscribe to store
- `frontend/src/panels/ViewerView.tsx` — subscribe to store
- `frontend/src/panels/ActivityBar.tsx`, `SidebarZone.tsx`, `BottomZone.tsx` — use `useDragHover`, render sliver
- `frontend/src/panels/FilesystemPanel.tsx` — use `useFilesystemProviders`, drop `!` assertions
- `frontend/src/hooks/useHasFilesystemProviders.ts` — wrap shared hook
- `frontend/src/hooks/useLeaveRoom.ts` — read store instead of singleton
- `frontend/src/hooks/socketHandlers/figureHandlers.ts`, `chatHandlers.ts` — read store
- `frontend/src/pages/landingPage.tsx` — use helpers + selector refactor
- `frontend/src/panels/ChatPanel.tsx` — fix exhaustive-deps
- `frontend/package.json`, `frontend/bun.lock` — drop `@vitejs/plugin-react`

### Frontend — deleted
- `frontend/e2e/dockview-layout.spec.ts`
- `frontend/e2e/pr920-fixes.spec.ts`

---

## Backend

### Task B1: `mock_executor` stub adds `token` param (Spec 2.6)

Trivial type fix. Pyright reports 5 call sites passing a 4-param coroutine where a 5-param `InternalExecutor` Protocol is expected. The stub exists deliberately for package-isolation tests — just add the missing param.

**Files:**
- Modify: `tests/zndraw_joblib/test_registry.py:57-63`

- [ ] **Step 1: Add `token: str` param to `mock_executor`**

Change lines 57-63 from:
```python
async def mock_executor(
    extension_cls: type[Extension],
    payload: dict[str, Any],
    room_id: str,
    task_id: str,
) -> None:
    pass
```
to:
```python
async def mock_executor(
    extension_cls: type[Extension],
    payload: dict[str, Any],
    room_id: str,
    task_id: str,
    token: str,
) -> None:
    pass
```

- [ ] **Step 2: Verify pyright clean on the file**

Run: `uv run pyright tests/zndraw_joblib/test_registry.py 2>&1 | grep -E "reportIncompatibleMethodOverride|reportArgumentType"`
Expected: no output (or no diagnostics on `mock_executor`-related lines).

- [ ] **Step 3: Run the affected tests**

Run: `uv run pytest tests/zndraw_joblib/test_registry.py -q`
Expected: all tests pass (behavior unchanged; only signature widened).

- [ ] **Step 4: Commit**

```bash
git add tests/zndraw_joblib/test_registry.py
git commit -m "test(providers): add token param to mock_executor stub

Aligns with InternalExecutor Protocol (registry.py:31-38) which gained
a token: str param when providers started requiring JWTs.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B2: Remove vestigial lazy `LoadFile` import (Spec 2.8)

No cycle exists (verified during brainstorm). The lazy import + `# noqa: E402` at the bottom of `modifiers.py` is pure comment rot.

**Files:**
- Modify: `src/zndraw/extensions/modifiers.py`

- [ ] **Step 1: Move the import to the top**

Current bottom of `src/zndraw/extensions/modifiers.py`:
```python
modifiers.update(molify_modifiers)

# LoadFile reads bytes from a registered filesystem provider and extends the
# room — must be registered as an @internal modifier so the server-side taskiq
# worker can execute it. Imported here (not at the top) to avoid a cycle via
# zndraw_joblib.client.
from zndraw.extensions.filesystem import LoadFile  # noqa: E402

modifiers[LoadFile.__name__] = LoadFile
```

Replace with: delete the comment and lazy import at the bottom, add `from zndraw.extensions.filesystem import LoadFile` near the top of the file after existing imports (alphabetized with other `zndraw.extensions.*` imports). The final lines become:
```python
modifiers.update(molify_modifiers)
modifiers[LoadFile.__name__] = LoadFile
```

- [ ] **Step 2: Verify imports work both ways**

Run: `uv run python -c "from zndraw.extensions.modifiers import modifiers; assert 'LoadFile' in modifiers, 'LoadFile missing'; print('OK')"`
Expected: `OK`

Run: `uv run python -c "from zndraw.extensions.filesystem import LoadFile; from zndraw.extensions.modifiers import modifiers; assert 'LoadFile' in modifiers; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Run extension tests**

Run: `uv run pytest tests/zndraw/test_extensions.py tests/zndraw -q -k "extension or modifier" --no-header 2>&1 | tail -20`
Expected: all collected tests pass.

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/extensions/modifiers.py
git commit -m "chore(extensions): remove vestigial lazy LoadFile import

filesystem.py imports only from zndraw.extensions.abc (which is stdlib+pydantic
only) — no cycle via zndraw_joblib.client exists. The noqa+comment are rot.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B3: `_resolve_provider` mirrors LIST semantics — security fix (Spec 2.3)

`@global` caller can currently resolve providers from any room due to an inverted short-circuit. Mirror the LIST allowlist logic.

**Files:**
- Modify: `src/zndraw_joblib/router.py:984-1005`
- Test: `tests/zndraw_joblib/test_providers.py`

- [ ] **Step 1: Write failing test — `@global` scope cannot resolve room-42 provider**

Add to `tests/zndraw_joblib/test_providers.py`:
```python
def test_global_scope_cannot_resolve_room_provider(client_factory):
    """A @global caller must not reach room-scoped providers (security: LIST
    policy excludes them; resolve must agree)."""
    alice = client_factory("alice", is_superuser=False)
    admin = client_factory("admin", is_superuser=True)

    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201

    # admin, calling with @global scope, must not resolve a room-42 provider.
    resp = admin.get(
        "/v1/joblib/rooms/@global/providers/room-42:filesystem:local"
    )
    assert resp.status_code == 404, resp.text


def test_global_scope_cannot_resolve_internal_provider(
    client_factory, async_session_factory
):
    """A @global caller must not reach @internal providers either — test
    pins the LIST-vs-RESOLVE symmetry documented in
    test_list_providers_includes_internal."""
    import asyncio
    import uuid

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord

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

    admin = client_factory("admin2", is_superuser=True)
    resp = admin.get(
        "/v1/joblib/rooms/@global/providers/@internal:filesystem:FilesystemRead"
    )
    assert resp.status_code == 404, resp.text
```

- [ ] **Step 2: Run the new tests — confirm they fail**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py::test_global_scope_cannot_resolve_room_provider tests/zndraw_joblib/test_providers.py::test_global_scope_cannot_resolve_internal_provider -q`
Expected: FAIL (responses currently resolve successfully).

- [ ] **Step 3: Fix `_resolve_provider`**

Replace `src/zndraw_joblib/router.py:984-1005` (the `_resolve_provider` function). Current:
```python
async def _resolve_provider(
    session: AsyncSession, provider_name: str, room_id: str
) -> ProviderRecord:
    """Resolve a provider full_name into a ProviderRecord, checking room visibility."""
    parts = provider_name.split(":", 2)
    if len(parts) != 3:
        raise ProviderNotFound.exception(
            detail=f"Invalid provider name format: {provider_name}"
        )
    provider_room_id, category, name = parts

    # Visibility check: provider must be in the room or @global
    if room_id not in ("@global",) and provider_room_id not in (
        "@global",
        "@internal",
        room_id,
    ):
        raise ProviderNotFound.exception(
            detail=f"Provider '{provider_name}' not accessible from room '{room_id}'"
        )
    ...
```

New:
```python
async def _resolve_provider(
    session: AsyncSession, provider_name: str, room_id: str
) -> ProviderRecord:
    """Resolve a provider full_name into a ProviderRecord, checking room visibility."""
    parts = provider_name.split(":", 2)
    if len(parts) != 3:
        raise ProviderNotFound.exception(
            detail=f"Invalid provider name format: {provider_name}"
        )
    provider_room_id, category, name = parts

    # Visibility check — mirror _room_provider_filter.
    allowed = (
        {"@global"}
        if room_id == "@global"
        else {"@global", "@internal", room_id}
    )
    if provider_room_id not in allowed:
        raise ProviderNotFound.exception(
            detail=f"Provider '{provider_name}' not accessible from room '{room_id}'"
        )
    ...
```

Leave the rest of the function (the actual DB query after the check) unchanged.

- [ ] **Step 4: Run the new tests — confirm they pass**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py::test_global_scope_cannot_resolve_room_provider tests/zndraw_joblib/test_providers.py::test_global_scope_cannot_resolve_internal_provider -q`
Expected: both PASS.

- [ ] **Step 5: Run the full provider test suite to check for regression**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py tests/zndraw_joblib/test_provider_dispatch.py -q`
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/zndraw_joblib/router.py tests/zndraw_joblib/test_providers.py
git commit -m "fix(providers): @global scope cannot resolve private providers

_resolve_provider's visibility check had an inverted short-circuit:
when caller room_id == '@global', the whole check was skipped, letting
a @global caller resolve providers from any room (Alice's room-42
filesystem provider was reachable via /rooms/@global/providers/...).

Mirror _room_provider_filter: @global callers see only @global.
Regression tests pin both the room-scoped and @internal cases.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B4: `StorageResultBackend` applies `key_prefix` (Spec 2.7)

Parity with `RedisResultBackend`. Two servers sharing a storage backend currently collide on provider cache keys.

**Files:**
- Modify: `src/zndraw/result_backends.py:151-185`
- Modify: `src/zndraw/database.py` (pass prefix through at composition time)
- Test: `tests/zndraw/test_result_backends.py`

- [ ] **Step 1: Write failing test — two backends with different prefixes do not collide**

Add to `tests/zndraw/test_result_backends.py`:
```python
import pytest

from zndraw.result_backends import StorageResultBackend
from zndraw.storage import get_frame_storage


@pytest.mark.asyncio
async def test_storage_result_backend_key_prefix_isolates_instances():
    """Two StorageResultBackends with different key_prefixes must not
    collide when sharing the same underlying FrameStorage."""
    shared_storage = get_frame_storage("memory://")
    backend_a = StorageResultBackend(shared_storage, key_prefix="server-a")
    backend_b = StorageResultBackend(shared_storage, key_prefix="server-b")

    await backend_a.store("result-key", b"payload-a", ttl=0)
    await backend_b.store("result-key", b"payload-b", ttl=0)

    assert await backend_a.get("result-key") == b"payload-a"
    assert await backend_b.get("result-key") == b"payload-b"

    await backend_a.delete("result-key")
    assert await backend_a.get("result-key") is None
    assert await backend_b.get("result-key") == b"payload-b"
```

- [ ] **Step 2: Run the new test — confirm it fails**

Run: `uv run pytest tests/zndraw/test_result_backends.py::test_storage_result_backend_key_prefix_isolates_instances -q`
Expected: FAIL — `TypeError: __init__() got an unexpected keyword argument 'key_prefix'` (or equivalent).

- [ ] **Step 3: Add `key_prefix` to `StorageResultBackend`**

Replace the `StorageResultBackend` class at `src/zndraw/result_backends.py:133-184`. Update `__init__` to accept `key_prefix`, add a `_k` helper, and apply it in `store`, `get`, `delete`:
```python
class StorageResultBackend:
    """Adapt ``FrameStorage`` to the ``ResultBackend`` protocol.

    Uses cache keys as room_ids and stores raw bytes as a single-entry
    frame (``{b"_": data}``).

    Parameters
    ----------
    storage
        The FrameStorage registry to delegate to.
    key_prefix
        Optional namespace prefix. When non-empty every storage key is
        stored as ``{key_prefix}:{key}``. Required for multi-server
        deployments sharing a single FrameStorage instance.
    """

    def __init__(self, storage: FrameStorage, key_prefix: str = "") -> None:
        self._storage = storage
        self._key_prefix = key_prefix

    def _k(self, key: str) -> str:
        """Return the namespaced key, prepending the prefix when set."""
        if self._key_prefix:
            return f"{self._key_prefix}:{key}"
        return key

    async def store(self, key: str, data: bytes, ttl: int) -> None:  # noqa: ARG002
        io = self._storage[self._k(key)]
        await io.clear()
        packed = msgpack.packb(data)
        assert packed is not None
        await io.extend([{b"_": packed}])

    async def get(self, key: str) -> bytes | None:
        io = self._storage[self._k(key)]
        if await io.len() == 0:
            return None
        packed = await io[b"_"][0]
        if packed is None:
            return None
        return msgpack.unpackb(packed)

    async def delete(self, key: str) -> None:
        await self._storage[self._k(key)].clear()

    async def acquire_inflight(self, key: str, ttl: int) -> bool:
        raise NotImplementedError("Use Redis for inflight locks")

    async def release_inflight(self, key: str) -> None:
        raise NotImplementedError("Use Redis for inflight locks")

    async def wait_for_key(self, key: str, timeout: float) -> bytes | None:  # noqa: ASYNC109
        raise NotImplementedError("Use CompositeResultBackend for long-polling")

    async def notify_key(self, key: str) -> None:
        raise NotImplementedError("Use CompositeResultBackend for long-polling")
```

- [ ] **Step 4: Wire `key_prefix` through at composition time**

Find the `StorageResultBackend(...)` construction site in `src/zndraw/database.py` (search `StorageResultBackend`). Change the constructor call to pass `key_prefix=settings.result_backend_key_prefix`:
```python
# before
frames_backend = StorageResultBackend(storage)
# after
frames_backend = StorageResultBackend(
    storage, key_prefix=settings.result_backend_key_prefix
)
```

- [ ] **Step 5: Run the new test — confirm it passes**

Run: `uv run pytest tests/zndraw/test_result_backends.py::test_storage_result_backend_key_prefix_isolates_instances -q`
Expected: PASS.

- [ ] **Step 6: Run the full result-backend suite**

Run: `uv run pytest tests/zndraw/test_result_backends.py -q`
Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/zndraw/result_backends.py src/zndraw/database.py tests/zndraw/test_result_backends.py
git commit -m "fix(backends): StorageResultBackend applies key_prefix

Parity with RedisResultBackend. Previously, two zndraw servers sharing
a FrameStorage instance (Mongo/LMDB) collided on provider cache keys
despite distinct result_backend_key_prefix settings — the prefix was
half-applied.

Wire settings.result_backend_key_prefix through at composition time.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B5: Lazy-resolve `worker_token` on `@internal` branch only (Spec 2.1)

Drop the route-level `WorkerTokenDep` from `read_provider`. Resolve manually inside the `@internal` branch. Refactor `get_worker_token` so it can be called from a plain helper (not only as a FastAPI dep).

**Files:**
- Modify: `src/zndraw_joblib/dependencies.py`
- Modify: `src/zndraw_joblib/router.py:1154-1220`
- Test: `tests/zndraw_joblib/test_providers.py`

- [ ] **Step 1: Inspect `get_worker_token`**

Run: `grep -n "get_worker_token\|internal_worker_user\|WorkerTokenDep" src/zndraw_joblib/dependencies.py`

Expected: shows `get_worker_token(request: Request) -> str` reading `request.app.state.internal_worker_user`. Note that `WorkerTokenDep = Annotated[str, Depends(get_worker_token)]`.

- [ ] **Step 2: Write failing test — remote provider read must work when `internal_worker_user` is `None`**

Add to `tests/zndraw_joblib/test_providers.py`:
```python
def test_read_remote_provider_works_with_no_internal_worker_cache(
    client_factory, monkeypatch_app_state
):
    """A remote provider read must not require app.state.internal_worker_user
    — that cache is only needed for @internal dispatch."""
    alice = client_factory("alice", is_superuser=False)

    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201

    # Simulate a deployment where init_db_on_startup=False and the worker
    # row seed hasn't run yet.
    monkeypatch_app_state("internal_worker_user", None)

    # Short-wait so the test doesn't hang — we only care the request does
    # not 500 before dispatch is attempted.
    resp = alice.get(
        "/v1/joblib/rooms/room-42/providers/room-42:filesystem:local?path=/",
        headers={"Prefer": "wait=0"},
    )
    # Acceptable outcomes:
    #  - 504 (no remote worker uploaded results in time — the desired path)
    #  - 409 NoWorkersAvailable (no connected worker — also acceptable)
    # Unacceptable: 500 (which is what happens today if WorkerTokenDep runs).
    assert resp.status_code in (504, 409), resp.text
```

You will likely need to add a `monkeypatch_app_state` fixture to the relevant `conftest.py` (check `tests/zndraw_joblib/conftest.py` first; add only if not present):
```python
@pytest.fixture
def monkeypatch_app_state(app, monkeypatch):
    """Temporarily override a key on app.state."""
    def _set(key: str, value):
        monkeypatch.setattr(app.state, key, value, raising=False)
    return _set
```

- [ ] **Step 3: Run the new test — confirm it fails**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py::test_read_remote_provider_works_with_no_internal_worker_cache -q`
Expected: FAIL with 500 status (WorkerTokenDep raises RuntimeError).

- [ ] **Step 4: Extract a plain helper in `dependencies.py`**

In `src/zndraw_joblib/dependencies.py`, keep `get_worker_token(request: Request) -> str` as the FastAPI dep but extract a plain `mint_internal_worker_token(app: FastAPI) -> str` helper:
```python
def mint_internal_worker_token(app) -> str:
    """Mint a JWT for the cached internal-worker user.

    Raises
    ------
    RuntimeError
        If ``app.state.internal_worker_user`` is not populated. The cache
        is primed during lifespan; if it is empty the DB row is missing
        (init_db_on_startup=False + the db-init step did not run).
    """
    user = getattr(app.state, "internal_worker_user", None)
    if user is None:
        raise RuntimeError(
            "internal_worker_user is not cached; ensure the db-init step ran"
        )
    # inline what get_worker_token used to do
    auth_settings = app.state.auth_settings
    return _build_jwt(user, auth_settings)


async def get_worker_token(request: Request) -> str:
    return mint_internal_worker_token(request.app)
```

Keep `WorkerTokenDep` pointing at `get_worker_token` — unchanged for any caller that explicitly wanted it as a dep.

- [ ] **Step 5: Remove `WorkerTokenDep` from `read_provider`, mint inside `@internal` branch**

Edit `src/zndraw_joblib/router.py:1154-1220`. Remove `worker_token: WorkerTokenDep` from the signature:
```python
@router.get("/rooms/{room_id}/providers/{provider_name:path}")
async def read_provider(
    room_id: str,
    provider_name: str,
    request: Request,
    session_maker: SessionMakerDep,
    _current_user: CurrentUserFactoryDep,
    result_backend: ResultBackendDep,
    settings: SettingsDep,
    tsio: TsioDep,
    internal_provider_registry: InternalProviderRegistryDep,
    prefer: Annotated[str | None, Header()] = None,
):
```

Inside the `@internal` branch (around line 1203-1220), mint the token just before the `.kiq(...)` call:
```python
if provider.room_id == "@internal":
    if (
        internal_provider_registry is None
        or provider.full_name not in internal_provider_registry.tasks
    ):
        raise InternalJobNotConfigured.exception(  # noqa: TRY301
            detail=(
                f"Internal provider '{provider.full_name}' is registered"
                " in the DB but no executor task is available"
            )
        )
    params_json = json.dumps(params, sort_keys=True, separators=(",", ":"))
    worker_token = mint_internal_worker_token(request.app)
    await internal_provider_registry.tasks[provider.full_name].kiq(
        request_id=rhash,
        provider_id=str(provider.id),
        params_json=params_json,
        token=worker_token,
    )
```

Import `mint_internal_worker_token` at the top of `router.py` alongside the other `dependencies` imports.

- [ ] **Step 6: Run the new test — confirm it passes**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py::test_read_remote_provider_works_with_no_internal_worker_cache -q`
Expected: PASS.

- [ ] **Step 7: Run the full provider suite**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py tests/zndraw_joblib/test_provider_dispatch.py tests/zndraw/test_providers_filesystem.py -q`
Expected: all tests pass.

- [ ] **Step 8: Commit**

```bash
git add src/zndraw_joblib/dependencies.py src/zndraw_joblib/router.py tests/zndraw_joblib/test_providers.py tests/zndraw_joblib/conftest.py
git commit -m "fix(providers): lazy-resolve worker_token for @internal branch only

read_provider previously declared worker_token: WorkerTokenDep as a
route-level dep. FastAPI resolved it on every call, but it was consumed
only in the @internal branch. On a fresh deploy with
init_db_on_startup=False and the worker seed not yet run, every provider
read 500'd for the process lifetime.

Extract mint_internal_worker_token() as a plain helper; call it inside
the @internal branch. Remote reads no longer depend on the cache.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B6: `ProviderExecutionFailed` ProblemType + preserve `X-Result-Status` (Spec 2.2+2.5, merged)

Replace the ad-hoc `{"error", "type"}` error body with RFC 9457 `ProblemDetail`. Preserve the `X-Result-Status` header through upload so the router routes on header, not content sniffing.

**Files:**
- Modify: `src/zndraw_joblib/exceptions.py`
- Modify: `src/zndraw/providers/executor.py`
- Modify: `src/zndraw_joblib/router.py:1180-1346`
- Test: `tests/zndraw/test_providers_executor.py` and/or `tests/zndraw_joblib/test_providers.py`

- [ ] **Step 1: Add `ProviderExecutionFailed` ProblemType**

Add to `src/zndraw_joblib/exceptions.py` (near the other provider problems, after `ProviderTimeout`):
```python
class ProviderExecutionFailed(ProblemType):
    """The provider failed to execute the requested read."""

    title: ClassVar[str] = "Bad Request"
    status: ClassVar[int] = 400
```

- [ ] **Step 2: Write failing test — legitimate JSON with `{"error", "type"}` keys must pass through unaltered**

Add to `tests/zndraw_joblib/test_providers.py`:
```python
def test_legitimate_json_with_error_type_keys_is_not_mis_flagged(
    client_factory, result_backend
):
    """A provider returning valid JSON with top-level 'error' and 'type'
    keys (e.g., JSON-Schema) must be returned as-is, not translated to
    an HTTP 400."""
    import asyncio

    alice = client_factory("alice", is_superuser=True)
    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201
    provider_id = resp.json()["id"]
    provider_full_name = resp.json()["full_name"]

    # Pre-populate the cache with legitimate JSON containing these keys.
    payload = json.dumps({"type": "object", "error": None, "ok": True}).encode()
    rhash = request_hash({"path": "/"})
    cache_key = f"provider-result:{provider_full_name}:{rhash}"

    async def _seed() -> None:
        await result_backend.store(cache_key, payload, 60)
    asyncio.run(_seed())

    resp = alice.get(
        f"/v1/joblib/rooms/room-42/providers/{provider_full_name}?path=/"
    )
    assert resp.status_code == 200, resp.text
    assert resp.json() == {"type": "object", "error": None, "ok": True}
```

- [ ] **Step 3: Run the new test — confirm it fails**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py::test_legitimate_json_with_error_type_keys_is_not_mis_flagged -q`
Expected: FAIL — response is 400 due to content-sniffing.

- [ ] **Step 4: Rewrite executor to emit `ProblemDetail` + set status header**

Replace the `_run()` body in `src/zndraw/providers/executor.py:58-93`. Key discipline: POST failures must propagate (they did in the original) — the outer try/except only wraps the *read work*, not the POST. Use a sentinel (`json_body is not None`) to route the POST:

```python
def _run() -> None:
    json_body: dict[str, Any] | None = None
    content: bytes | None = None

    try:
        handler = self._resolve_handler(provider_cls)
        params = json.loads(params_json) if params_json else {}
        instance = provider_cls(**params)
        result = instance.read(handler)
        if provider_cls.content_type == "application/json":
            content = json.dumps(result).encode()
        else:
            content = result  # type: ignore[assignment]
        headers = {
            "Authorization": f"Bearer {token}",
            "X-Request-Hash": request_id,
        }
    except Exception as err:
        log.exception(
            "InternalProviderExecutor failed for %s",
            provider_cls.__name__,
        )
        from zndraw_joblib.exceptions import ProviderExecutionFailed
        problem = ProviderExecutionFailed.create(
            detail=f"{type(err).__name__}: {err}",
        )
        json_body = problem.model_dump(exclude_none=True)
        headers = {
            "Authorization": f"Bearer {token}",
            "X-Request-Hash": request_id,
            "X-Result-Status": "error",
            "Content-Type": "application/problem+json",
        }

    client_kwargs: dict[str, Any] = {"timeout": self.timeout_seconds}
    if transport is not None:
        client_kwargs["transport"] = transport

    with httpx.Client(**client_kwargs) as client:
        if json_body is not None:
            resp = client.post(
                f"{base_url}/v1/joblib/providers/{provider_id}/results",
                json=json_body,
                headers=headers,
            )
        else:
            resp = client.post(
                f"{base_url}/v1/joblib/providers/{provider_id}/results",
                content=content,
                headers=headers,
            )
        resp.raise_for_status()
```

Note: the outer try/except only wraps `instance.read(handler)` and its serialization. If the POST fails (network, 5xx from server), the error propagates to taskiq — matching the original. A POST failure is not re-POSTed as an executor-failed payload.

- [ ] **Step 5: Update `upload_result` to preserve `X-Result-Status` via paired Redis key**

Edit `src/zndraw_joblib/router.py:upload_result` (around line 1300-1346). Add a `x_result_status` header param and store a sibling key:
```python
@router.post("/providers/{provider_id}/results", status_code=status.HTTP_204_NO_CONTENT)
async def upload_result(
    provider_id: UUID,
    request: Request,
    session_maker: SessionMakerDep,
    user: CurrentUserFactoryDep,
    result_backend: ResultBackendDep,
    settings: SettingsDep,
    tsio: TsioDep,
    x_request_hash: Annotated[str, Header()],
    x_result_status: Annotated[str | None, Header()] = None,
):
    """Provider worker uploads a read result."""
    async with session_maker() as session:
        result = await session.exec(
            select(ProviderRecord).where(ProviderRecord.id == provider_id)
        )
        provider = result.one_or_none()
        if not provider:
            raise ProviderNotFound.exception(
                detail=f"Provider '{provider_id}' not found"
            )

        if provider.user_id != user.id and not user.is_superuser:
            raise Forbidden.exception(
                detail="Not authorized to upload results for this provider"
            )

    cache_key = f"provider-result:{provider.full_name}:{x_request_hash}"
    status_key = f"{cache_key}:status"
    inflight_key = f"provider-inflight:{provider.full_name}:{x_request_hash}"

    data = await request.body()
    await result_backend.store(
        cache_key,
        data,
        settings.provider_result_ttl_seconds,
    )
    if x_result_status == "error":
        await result_backend.store(
            status_key,
            b"error",
            settings.provider_result_ttl_seconds,
        )

    await result_backend.release_inflight(inflight_key)
    await result_backend.notify_key(cache_key)

    await emit(
        tsio,
        {
            Emission(
                ProviderResultReady(
                    provider_name=provider.full_name,
                    request_hash=x_request_hash,
                ),
                f"room:{provider.room_id}",
            )
        },
    )
```

- [ ] **Step 6: Update `read_provider` to route on status header, not content sniff**

Replace both sniff branches in `src/zndraw_joblib/router.py` (around 1180-1194 and 1248-1268). The fast-path cache read becomes:
```python
status_key = f"{cache_key}:status"

cached = await result_backend.get(cache_key)
if cached is not None:
    cached_status = await result_backend.get(status_key)
    if cached_status == b"error":
        # Parse the ProblemDetail to recover the original status code.
        try:
            problem = json.loads(cached)
            status_code = int(problem.get("status", 500))
        except (ValueError, TypeError, KeyError):
            status_code = 500
        return Response(
            content=cached,
            media_type="application/problem+json",
            status_code=status_code,
        )
    return Response(content=cached, media_type=provider.content_type)
```

The long-poll path (around line 1248-1268):
```python
result = await result_backend.wait_for_key(cache_key, timeout)
if result is not None:
    result_status = await result_backend.get(status_key)
    headers = {}
    if requested_wait is not None:
        headers["Preference-Applied"] = f"wait={int(timeout)}"
    if result_status == b"error":
        try:
            problem = json.loads(result)
            status_code = int(problem.get("status", 500))
        except (ValueError, TypeError, KeyError):
            status_code = 500
        return Response(
            content=result,
            media_type="application/problem+json",
            status_code=status_code,
            headers=headers,
        )
    return Response(
        content=result, media_type=provider.content_type, headers=headers
    )
```

Remove the existing content-sniffing blocks entirely.

- [ ] **Step 7: Add a test for the error path (ProblemDetail round-trip)**

Add to `tests/zndraw_joblib/test_providers.py`:
```python
def test_provider_error_path_returns_problem_detail(
    client_factory, result_backend
):
    """When an executor posts an error, read_provider returns RFC 9457
    problem+json with the status from the payload (not a hard-coded 400)."""
    import asyncio

    from zndraw_joblib.exceptions import ProviderExecutionFailed

    alice = client_factory("alice", is_superuser=True)
    resp = alice.put(
        "/v1/joblib/rooms/room-42/providers",
        json={"category": "filesystem", "name": "local", "schema": {}},
    )
    assert resp.status_code == 201
    provider_full_name = resp.json()["full_name"]

    problem = ProviderExecutionFailed.create(detail="FileNotFoundError: /nope")
    payload = problem.model_dump_json(exclude_none=True).encode()
    rhash = request_hash({"path": "/nope"})
    cache_key = f"provider-result:{provider_full_name}:{rhash}"

    async def _seed() -> None:
        await result_backend.store(cache_key, payload, 60)
        await result_backend.store(f"{cache_key}:status", b"error", 60)
    asyncio.run(_seed())

    resp = alice.get(
        f"/v1/joblib/rooms/room-42/providers/{provider_full_name}?path=/nope"
    )
    assert resp.status_code == 400
    assert resp.headers["content-type"].startswith("application/problem+json")
    body = resp.json()
    assert body["title"] == "Bad Request"
    assert body["status"] == 400
    assert "FileNotFoundError" in body["detail"]
```

- [ ] **Step 8: Run all three tests — confirm they pass**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py -q -k "legitimate_json or error_path_returns_problem_detail"`
Expected: all PASS.

- [ ] **Step 9: Run the executor integration test**

Run: `uv run pytest tests/zndraw/test_providers_executor.py tests/zndraw/worker/test_internal_loadfile_e2e.py -q`
Expected: all tests pass. If an assertion depends on the old `{"error", "type"}` shape, update it to match the new `ProblemDetail` shape.

- [ ] **Step 10: Run the full suite**

Run: `uv run pytest tests/zndraw tests/zndraw_joblib -q`
Expected: all tests pass.

- [ ] **Step 11: Commit**

```bash
git add src/zndraw_joblib/exceptions.py src/zndraw/providers/executor.py src/zndraw_joblib/router.py tests/zndraw_joblib/test_providers.py tests/zndraw/test_providers_executor.py
git commit -m "fix(providers): ProviderExecutionFailed + preserve X-Result-Status

The executor emitted {\"error\", \"type\"} and the router content-sniffed
that shape to translate to HTTP 400. Two bugs:
  1. Legitimate payloads with those two top-level keys (JSON-Schema,
     CloudEvents) got mis-flagged as errors.
  2. The error shape diverged from the codebase's RFC 9457 ProblemDetail.

Replace with ProviderExecutionFailed(ProblemType). Executor POSTs via
httpx json= + X-Result-Status header + application/problem+json.
upload_result preserves the header to a sibling :status Redis key.
read_provider routes on the status key, not content sniffing.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B7: `filebrowser_require_superuser` gate (Spec 2.4)

New Pydantic setting. Non-superusers blocked from `@internal:filesystem:*` in three endpoints: list, info, read.

**Files:**
- Modify: `src/zndraw/config.py`
- Modify: `src/zndraw_joblib/router.py`
- Test: `tests/zndraw_joblib/test_providers.py`

- [ ] **Step 1: Add the setting**

Insert into `src/zndraw/config.py` after `filebrowser_path`:
```python
filebrowser_require_superuser: bool = True
"""When True, @internal filesystem providers are accessible only to
superusers. Flip to False to allow all authenticated users.

Secure default — prevents any authenticated user in any room from
reading the directory at ``filebrowser_path``."""
```

- [ ] **Step 2: Write failing test — non-superuser cannot read @internal filesystem**

Add to `tests/zndraw_joblib/test_providers.py`:
```python
def test_internal_filesystem_requires_superuser_by_default(
    client_factory, async_session_factory
):
    """With the default filebrowser_require_superuser=True, a non-superuser
    is 403'd on @internal:filesystem:*."""
    import asyncio
    import uuid

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int2@test",
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

    alice = client_factory("alice", is_superuser=False)

    # Read
    resp = alice.get(
        "/v1/joblib/rooms/room-42/providers/@internal:filesystem:FilesystemRead?path=/"
    )
    assert resp.status_code == 403, resp.text

    # List — no @internal:filesystem:* entries
    resp = alice.get("/v1/joblib/rooms/room-42/providers")
    assert resp.status_code == 200
    items = resp.json()["items"]
    assert not any(
        p["full_name"].startswith("@internal:filesystem:") for p in items
    )

    # Info endpoint gate — add an assertion here only after Step 6 below
    # lands and you have confirmed the info endpoint's path via grep.


def test_internal_filesystem_superuser_can_read(
    client_factory, async_session_factory
):
    """Superusers bypass the gate."""
    import asyncio
    import uuid

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int3@test",
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

    admin = client_factory("admin-su", is_superuser=True)
    resp = admin.get("/v1/joblib/rooms/room-42/providers")
    items = resp.json()["items"]
    assert any(
        p["full_name"] == "@internal:filesystem:FilesystemRead" for p in items
    )
```

- [ ] **Step 3: Run the failing tests**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py -q -k "requires_superuser or superuser_can_read"`
Expected: `test_internal_filesystem_requires_superuser_by_default` FAILs on the read assertion (403 expected, 200/504 returned).

- [ ] **Step 4: Add `_require_internal_filesystem_access` helper**

Add to `src/zndraw_joblib/router.py` near `_resolve_provider`:
```python
def _require_internal_filesystem_access(
    provider: ProviderRecord, user, settings
) -> None:
    """Gate @internal:filesystem:* access on superuser status.

    Raises
    ------
    ProblemError
        403 Forbidden when the caller is non-superuser and the flag is on.
    """
    if (
        provider.room_id == "@internal"
        and provider.category == "filesystem"
        and settings.filebrowser_require_superuser
        and not user.is_superuser
    ):
        raise Forbidden.exception(
            detail="@internal filesystem access requires superuser"
        )
```

- [ ] **Step 5: Enforce at `read_provider`**

In `src/zndraw_joblib/router.py:read_provider`, after the `_resolve_provider` call (line ~1173) add:
```python
_require_internal_filesystem_access(provider, _current_user, settings)
```

Note: `_current_user` is already in scope as a dep. `settings` is already in scope.

- [ ] **Step 6: Enforce at `get_provider_info`**

Locate the provider-info endpoint around `src/zndraw_joblib/router.py:1149` (uses `_resolve_provider`). Ensure the handler has `_current_user: CurrentUserFactoryDep` and `settings: SettingsDep` in its signature; add them if missing. After the `_resolve_provider` call, add:
```python
_require_internal_filesystem_access(provider, _current_user, settings)
```

- [ ] **Step 7: Post-filter list results**

In the `_list_providers` endpoint (grep for `@router.get("/rooms/{room_id}/providers")` in router.py — note the path without the provider_name segment), after fetching the results but before returning, filter:
```python
if settings.filebrowser_require_superuser and not _current_user.is_superuser:
    items = [
        p for p in items
        if not (p.room_id == "@internal" and p.category == "filesystem")
    ]
```

Handle the `total` count accordingly (re-compute from filtered items).

- [ ] **Step 8: Run the new tests — confirm they pass**

Run: `uv run pytest tests/zndraw_joblib/test_providers.py -q -k "requires_superuser or superuser_can_read"`
Expected: both PASS.

- [ ] **Step 9: Add a test for the flag=False escape hatch**

Add to `tests/zndraw_joblib/test_providers.py`:
```python
def test_internal_filesystem_gate_disabled(
    client_factory, async_session_factory, app
):
    """With filebrowser_require_superuser=False, non-superusers can access."""
    import asyncio
    import uuid

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord

    app.state.settings.filebrowser_require_superuser = False

    async def seed() -> None:
        async with async_session_factory() as session:
            user = User(
                id=uuid.uuid4(),
                email="int4@test",
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

    alice = client_factory("alice-gated-off", is_superuser=False)
    resp = alice.get("/v1/joblib/rooms/room-42/providers")
    items = resp.json()["items"]
    assert any(
        p["full_name"] == "@internal:filesystem:FilesystemRead" for p in items
    )

    # Reset for other tests
    app.state.settings.filebrowser_require_superuser = True
```

Run: `uv run pytest tests/zndraw_joblib/test_providers.py::test_internal_filesystem_gate_disabled -q`
Expected: PASS.

- [ ] **Step 10: Run the full provider suite**

Run: `uv run pytest tests/zndraw tests/zndraw_joblib -q`
Expected: all tests pass.

- [ ] **Step 11: Commit**

```bash
git add src/zndraw/config.py src/zndraw_joblib/router.py tests/zndraw_joblib/test_providers.py
git commit -m "feat(providers): filebrowser_require_superuser gate

Secure default (True) — @internal:filesystem:* is restricted to
superusers. Operators who want broader access flip the flag.

Enforced in three endpoints: _list_providers (post-filter),
get_provider_info, read_provider. Helper _require_internal_filesystem_access
centralizes the policy.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task B8: Split `filebrowser_enabled` / `filebrowser_path` + clean up tests (Spec 2.9 + 2.10)

Eliminate the `env_parse_none_str` sentinel. Two explicit fields: one to enable, one for the path.

**Files:**
- Modify: `src/zndraw/config.py`
- Modify: `src/zndraw/database.py` (or wherever `filebrowser_path is None` is checked)
- Modify: `src/zndraw/providers/bootstrap.py`
- Modify: `tests/zndraw/test_providers_filesystem.py`
- Modify: any frontend-relevant check (the icon hiding hook reads from API, not config, so usually none)

- [ ] **Step 1: Audit current callsites**

Run: `grep -rn "filebrowser_path" src/ tests/ | grep -v docs`
Note every site that reads or checks this setting. Typical places: `src/zndraw/providers/bootstrap.py`, `src/zndraw/database.py`, `src/zndraw/providers/executor.py` (the `filebrowser_path` constructor arg — keep as-is; it's a value, not a None-check).

- [ ] **Step 2: Update `config.py`**

Replace in `src/zndraw/config.py`:
```python
# before
model_config = SettingsConfigDict(
    env_prefix="ZNDRAW_SERVER_",
    pyproject_toml_table_header=("tool", "zndraw", "server"),
    env_parse_none_str="none",
)
...
filebrowser_path: str | None = "."
"""Path for the default @internal filesystem provider.
... (long docstring about 'none' sentinel) ...
"""
```
with:
```python
model_config = SettingsConfigDict(
    env_prefix="ZNDRAW_SERVER_",
    pyproject_toml_table_header=("tool", "zndraw", "server"),
)
...
filebrowser_enabled: bool = True
"""Whether the default @internal filesystem provider is registered.

Set ``ZNDRAW_SERVER_FILEBROWSER_ENABLED=false`` to disable — no DB row,
no task registration, and the frontend filesystem activity-bar icon
is hidden (via the providers list being empty for that category)."""

filebrowser_path: str = "."
"""Path rooting the default @internal filesystem provider. Ignored when
``filebrowser_enabled`` is False."""
```

Keep `filebrowser_require_superuser` from Task B7.

- [ ] **Step 3: Update callsites**

For every `settings.filebrowser_path is None` or `if settings.filebrowser_path:` check, replace with `if not settings.filebrowser_enabled:` / `if settings.filebrowser_enabled:`. Typical:
```python
# before
if settings.filebrowser_path is not None:
    ensure_internal_providers(...)
# after
if settings.filebrowser_enabled:
    ensure_internal_providers(...)
```

For places that pass the path to the executor, keep using `settings.filebrowser_path` (it's a path, not a None-check).

- [ ] **Step 4: Update tests — delete weak test, rename remaining two**

In `tests/zndraw/test_providers_filesystem.py`:

Delete `test_filebrowser_path_none_disables_default_provider` (lines 360-368).

Rename and rewrite `test_filebrowser_path_none_disables_default` → `test_filebrowser_disabled_hides_default_provider`:
```python
def test_filebrowser_disabled_hides_default_provider(server_factory):
    """ZNDRAW_SERVER_FILEBROWSER_ENABLED=false drops the @internal provider."""
    instance = server_factory({"ZNDRAW_SERVER_FILEBROWSER_ENABLED": "false"})
    vis = ZnDraw(url=instance.url)
    try:
        providers = _list_providers(vis)
        assert not any(p.room_id == "@internal" for p in providers)
    finally:
        vis.disconnect()
```

Rename `test_filebrowser_path_none_removes_stale_rows` → `test_filebrowser_disabled_removes_stale_rows`. Update the env-var payload in each of its two boots:
```python
# first boot — enabled
server_factory({
    "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "true",
    "ZNDRAW_SERVER_FILEBROWSER_PATH": ".",
    "ZNDRAW_SERVER_DATABASE_URL": db_url,
})
# second boot — disabled
server_factory({
    "ZNDRAW_SERVER_FILEBROWSER_ENABLED": "false",
    "ZNDRAW_SERVER_DATABASE_URL": db_url,
})
```

- [ ] **Step 5: Write a regression test — `guest_password="none"` stays literal**

Add to `tests/zndraw/test_config.py`:
```python
def test_guest_password_literal_none_not_coerced(monkeypatch):
    """Dropping env_parse_none_str means 'none' is a literal string
    everywhere, not an implicit None sentinel."""
    from pydantic import SecretStr

    from zndraw.config import Settings

    monkeypatch.setenv("ZNDRAW_SERVER_GUEST_PASSWORD", "none")
    s = Settings()
    assert isinstance(s.guest_password, SecretStr)
    assert s.guest_password.get_secret_value() == "none"
```

- [ ] **Step 6: Run filesystem + config suites**

Run: `uv run pytest tests/zndraw/test_providers_filesystem.py tests/zndraw/test_config.py -q`
Expected: all tests pass, including the rename.

- [ ] **Step 7: Run the full suite**

Run: `uv run pytest tests/zndraw tests/zndraw_joblib -q`
Expected: all tests pass. Fix any test that still sets `ZNDRAW_SERVER_FILEBROWSER_PATH=none`.

- [ ] **Step 8: Commit**

```bash
git add src/zndraw/config.py src/zndraw/database.py src/zndraw/providers/bootstrap.py tests/zndraw/test_providers_filesystem.py tests/zndraw/test_config.py
git commit -m "refactor(config): split filebrowser_enabled / filebrowser_path

Eliminate env_parse_none_str='none' sentinel. Two explicit fields:
  filebrowser_enabled: bool = True
  filebrowser_path: str = '.'

Removes the model-wide footgun where any string field accepting
'none' (e.g., guest_password) was silently coerced to None.

Delete the weak tautological test; rename remaining tests to match
the new env-var names.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Frontend

### Task F1: Delete temporary e2e spec files (Spec 3.1)

**Files:**
- Delete: `frontend/e2e/dockview-layout.spec.ts`
- Delete: `frontend/e2e/pr920-fixes.spec.ts`

- [ ] **Step 1: Delete the files**

```bash
rm frontend/e2e/dockview-layout.spec.ts frontend/e2e/pr920-fixes.spec.ts
```

- [ ] **Step 2: Verify the remaining e2e suite is unaffected**

Run: `ls frontend/e2e/*.spec.ts`
Expected: 11 remaining spec files (no dockview/pr920 entries).

- [ ] **Step 3: Commit**

```bash
git add -u frontend/e2e/
git commit -m "fix(frontend): remove temporary e2e spec files

dockview-layout.spec.ts and pr920-fixes.spec.ts referenced selectors
(data-sliver-state, text 'Drop to dock right') that do not exist in
the implementation. Frontend e2e is not in CI; the files were
misleading-docs masquerading as coverage.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task F2: Drop unused `@vitejs/plugin-react` dep (Spec 3.3)

**Files:**
- Modify: `frontend/package.json`
- Modify: `frontend/bun.lock`

- [ ] **Step 1: Remove the dep**

Edit `frontend/package.json` — delete the line `"@vitejs/plugin-react": "^4",` from `devDependencies` (keep `"@vitejs/plugin-react-swc": "^4.2.1"`).

- [ ] **Step 2: Regenerate lockfile**

```bash
cd frontend && bun install
```
Expected: `bun.lock` updates, ~7 transitive Babel packages removed.

- [ ] **Step 3: Verify build and dev still work**

Run: `cd frontend && bun run build 2>&1 | tail -5`
Expected: succeeds with `dist/` output.

Run: `cd frontend && bunx vite --version` (dev smoke)
Expected: version prints, no error.

- [ ] **Step 4: Commit**

```bash
git add frontend/package.json frontend/bun.lock
git commit -m "chore(frontend): drop unused @vitejs/plugin-react dep

vite.config.ts uses @vitejs/plugin-react-swc; the non-SWC variant was
installed alongside but never imported. Dropping removes ~7 transitive
Babel packages.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task F3: Zustand `useDockviewApi` store replaces `sharedApi` (Spec 3.2 + 3.6)

Creates the store. All consumers migrate from `getDockviewApi()` (synchronous) to `useDockviewApi((s) => s.api)` (reactive) or `useDockviewApi.getState().api` (imperative).

**Files:**
- Create: `frontend/src/stores/dockviewApiStore.ts`
- Modify: `frontend/src/panels/DockviewLayout.tsx`
- Modify: `frontend/src/panels/PlotsBrowserPanel.tsx`
- Modify: `frontend/src/panels/ViewerView.tsx`
- Modify: `frontend/src/hooks/useLeaveRoom.ts`
- Modify: `frontend/src/hooks/socketHandlers/figureHandlers.ts`
- Modify: `frontend/src/hooks/socketHandlers/chatHandlers.ts`

- [ ] **Step 1: Create the store**

Write `frontend/src/stores/dockviewApiStore.ts`:
```ts
import type { DockviewApi } from "dockview-react";
import { create } from "zustand";

interface DockviewApiStore {
	api: DockviewApi | null;
	setApi: (api: DockviewApi | null) => void;
}

export const useDockviewApi = create<DockviewApiStore>((set) => ({
	api: null,
	setApi: (api) => set({ api }),
}));
```

- [ ] **Step 2: Write failing test on the store**

Write `frontend/src/stores/__tests__/dockviewApiStore.test.ts`:
```ts
import { beforeEach, describe, expect, it } from "vitest";
import { useDockviewApi } from "../dockviewApiStore";

describe("dockviewApiStore", () => {
	beforeEach(() => {
		useDockviewApi.getState().setApi(null);
	});

	it("defaults to null", () => {
		expect(useDockviewApi.getState().api).toBeNull();
	});

	it("setApi updates state", () => {
		const fakeApi = { panels: [] } as unknown as import("dockview-react").DockviewApi;
		useDockviewApi.getState().setApi(fakeApi);
		expect(useDockviewApi.getState().api).toBe(fakeApi);
	});

	it("setApi(null) clears", () => {
		const fakeApi = {} as unknown as import("dockview-react").DockviewApi;
		useDockviewApi.getState().setApi(fakeApi);
		useDockviewApi.getState().setApi(null);
		expect(useDockviewApi.getState().api).toBeNull();
	});
});
```

- [ ] **Step 3: Run the test — confirm it passes**

Run: `cd frontend && bun run test src/stores/__tests__/dockviewApiStore.test.ts`
Expected: 3 passed.

- [ ] **Step 4: Migrate `DockviewLayout.tsx`**

Edit `frontend/src/panels/DockviewLayout.tsx`:

Remove the module-level singleton (lines 25-29):
```ts
// DELETE:
let sharedApi: DockviewApi | null = null;
export function getDockviewApi(): DockviewApi | null {
	return sharedApi;
}
```

Add the import:
```ts
import { useDockviewApi } from "../stores/dockviewApiStore";
```

In `onReady` (line ~57-76), replace `sharedApi = event.api;` with:
```ts
useDockviewApi.getState().setApi(event.api);
```

In the cleanup effect (line 82), replace `sharedApi = null;` with:
```ts
useDockviewApi.getState().setApi(null);
```

- [ ] **Step 5: Migrate `PlotsBrowserPanel.tsx`**

Edit `frontend/src/panels/PlotsBrowserPanel.tsx`. Replace the polling `useOpenPlotKeys` hook (lines 18-64) with a store-subscribing version:
```ts
import { useDockviewApi } from "../stores/dockviewApiStore";

function useOpenPlotKeys(): Set<string> {
	const api = useDockviewApi((s) => s.api);
	const [version, setVersion] = useState(0);

	useEffect(() => {
		if (!api) return;
		const disposables = [
			api.onDidAddPanel(() => setVersion((v) => v + 1)),
			api.onDidRemovePanel(() => setVersion((v) => v + 1)),
		];
		return () => {
			for (const d of disposables) d.dispose();
		};
	}, [api]);

	return useMemo(() => {
		if (!api) return new Set<string>();
		return new Set(
			api.panels
				.filter((p) => p.id.startsWith("plot-"))
				.map((p) => p.id.slice("plot-".length)),
		);
		// biome-ignore lint/correctness/useExhaustiveDependencies: version is a tick trigger
	}, [api, version]);
}
```

In the component body, replace `getDockviewApi()` with the hook value in `onRowClick` — but the hook returns at the component level, so refactor `PlotsBrowserPanel`:
```ts
export function PlotsBrowserPanel() {
	const api = useDockviewApi((s) => s.api);
	const { data, isLoading } = useFigureList();
	const openKeys = useOpenPlotKeys();
	const allKeys = useMemo(() => data?.items ?? [], [data]);

	const onRowClick = (key: string) => {
		if (!api) return;
		openPlotTab(api, key);
	};
	// ... rest unchanged, but any other getDockviewApi() usage uses `api`
}
```

Remove `import { getDockviewApi } from "./DockviewLayout";`.

- [ ] **Step 6: Migrate `ViewerView.tsx`**

Edit `frontend/src/panels/ViewerView.tsx`. Replace `getDockviewApi()` with the reactive hook:
```ts
import { useDockviewApi } from "../stores/dockviewApiStore";

export function ViewerView() {
	const api = useDockviewApi((s) => s.api);
	// ... existing hook references and refs ...
	useEffect(() => {
		if (!api) return;
		const dispose = api.onDidRemovePanel(/* ... existing handler ... */);
		return () => dispose.dispose();
	}, [api]);
	// ...
}
```

- [ ] **Step 7: Migrate `useLeaveRoom.ts`**

Edit `frontend/src/hooks/useLeaveRoom.ts`. Replace any `getDockviewApi()` call with `useDockviewApi.getState().api` (imperative, called from a callback — not a subscription).

- [ ] **Step 8: Migrate socket handlers**

In `frontend/src/hooks/socketHandlers/figureHandlers.ts` and `chatHandlers.ts`, replace `getDockviewApi()` with `useDockviewApi.getState().api`. Handlers are not React components; they cannot use hooks, and they are called at event time — `getState()` is correct.

- [ ] **Step 9: Ensure no `getDockviewApi` references remain**

Run: `grep -rn "getDockviewApi\|sharedApi" frontend/src/`
Expected: no output.

- [ ] **Step 10: Type-check and test**

Run: `cd frontend && bunx tsc --noEmit 2>&1 | tail -20`
Expected: no errors.

Run: `cd frontend && bun run test`
Expected: all tests pass (including the new store test).

- [ ] **Step 11: Commit**

```bash
git add frontend/src/stores/dockviewApiStore.ts frontend/src/stores/__tests__/dockviewApiStore.test.ts frontend/src/panels/DockviewLayout.tsx frontend/src/panels/PlotsBrowserPanel.tsx frontend/src/panels/ViewerView.tsx frontend/src/hooks/useLeaveRoom.ts frontend/src/hooks/socketHandlers/figureHandlers.ts frontend/src/hooks/socketHandlers/chatHandlers.ts
git commit -m "refactor(frontend): useDockviewApi zustand store

Replace the 'let sharedApi' module singleton + setInterval polling
with a Zustand store. Consumers:
  - Components subscribe via useDockviewApi((s) => s.api) and re-render
    when the api binds — no more polling.
  - Non-component callers (socket handlers) use .getState().api.

Also closes the ViewerView.onLeaveRoom stale-reference race: the
subscription semantics replace getting a null singleton mid-cleanup.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task F4: Extract `resetDockview` + `ensureViewerPanel` helpers (Spec 3.7)

**Files:**
- Modify: `frontend/src/panels/DockviewLayout.tsx`
- Modify: `frontend/src/pages/landingPage.tsx:219-230,459-473`

- [ ] **Step 1: Add the helpers to `DockviewLayout.tsx`**

Add below the private `addViewerPanel` helper in `frontend/src/panels/DockviewLayout.tsx`:
```ts
export function resetDockview(api: DockviewApi): void {
	for (const p of api.panels) p.api.close();
	addViewerPanel(api);
}

export function ensureViewerPanel(api: DockviewApi): void {
	if (!api.getPanel("viewer")) addViewerPanel(api);
}
```

- [ ] **Step 2: Migrate `landingPage.tsx:219-230` useEffect**

Edit the roomId-change useEffect to use the helper. Current:
```ts
useEffect(() => {
	const api = getDockviewApi();
	if (!api) return;
	if (!roomId) return;
	if (!api.getPanel("viewer")) {
		api.addPanel({
			id: "viewer",
			component: "viewer",
			title: "3D Viewer",
		});
	}
}, [roomId]);
```

(Note: after Task F3, `getDockviewApi` is gone; replace with store read.)

New:
```ts
const dockApi = useDockviewApi((s) => s.api);
useEffect(() => {
	if (!dockApi) return;
	if (!roomId) return;
	ensureViewerPanel(dockApi);
}, [roomId, dockApi]);
```

- [ ] **Step 3: Migrate `landingPage.tsx:459-473` Reset-layout MenuItem**

Current:
```ts
<MenuItem
	onClick={() => {
		handleProfileClose();
		useAppStore.getState().resetLayout();
		const api = getDockviewApi();
		if (api) {
			for (const p of api.panels) {
				p.api.close();
			}
			api.addPanel({
				id: "viewer",
				component: "viewer",
				title: "3D Viewer",
			});
		}
	}}
>
	<ListItemText>Reset layout</ListItemText>
</MenuItem>
```

New:
```ts
<MenuItem
	onClick={() => {
		handleProfileClose();
		useAppStore.getState().resetLayout();
		const api = useDockviewApi.getState().api;
		if (api) resetDockview(api);
	}}
>
	<ListItemText>Reset layout</ListItemText>
</MenuItem>
```

Add imports at the top of `landingPage.tsx`:
```ts
import { ensureViewerPanel, resetDockview } from "../panels/DockviewLayout";
import { useDockviewApi } from "../stores/dockviewApiStore";
```

- [ ] **Step 4: Type-check and smoke test**

Run: `cd frontend && bunx tsc --noEmit 2>&1 | tail -10`
Expected: no errors.

Run: `cd frontend && bun run test`
Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add frontend/src/panels/DockviewLayout.tsx frontend/src/pages/landingPage.tsx
git commit -m "refactor(panels): resetDockview + ensureViewerPanel helpers

Three paths previously redeclared the viewer panel config (id,
component, title): DockviewLayout.onReady, roomId useEffect, and the
Reset-layout MenuItem. Extract two helpers colocated with addViewerPanel;
consumers import and call.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task F5: Shared `useFilesystemProviders` hook (Spec 3.8)

**Files:**
- Create: `frontend/src/hooks/useFilesystemProviders.ts`
- Modify: `frontend/src/hooks/useHasFilesystemProviders.ts`
- Modify: `frontend/src/panels/FilesystemPanel.tsx`

- [ ] **Step 1: Create the shared hook**

Write `frontend/src/hooks/useFilesystemProviders.ts`:
```ts
import { useQuery, useQueryClient } from "@tanstack/react-query";
import { useEffect } from "react";
import { listProviders } from "../myapi/client";
import { socket } from "../socket";
import { useAppStore } from "../store";

/**
 * Query hook for filesystem providers in the current room.
 *
 * Sets ``staleTime: 5_000`` and subscribes to the
 * ``providers_invalidate`` Socket.IO event so a newly-registered
 * provider triggers a refetch without a page reload. Shared by
 * ``useHasFilesystemProviders`` and ``FilesystemPanel`` so both
 * consumers agree on cache policy and invalidation.
 */
export function useFilesystemProviders() {
	const roomId = useAppStore((s) => s.roomId);
	const queryClient = useQueryClient();

	const result = useQuery({
		queryKey: ["filesystemProviders", roomId],
		queryFn: () => listProviders(roomId!, "filesystem"),
		enabled: !!roomId,
		retry: false,
		staleTime: 5_000,
	});

	useEffect(() => {
		if (!roomId) return;
		const handle = () => {
			queryClient.invalidateQueries({
				queryKey: ["filesystemProviders", roomId],
			});
		};
		socket.on("providers_invalidate", handle);
		return () => {
			socket.off("providers_invalidate", handle);
		};
	}, [roomId, queryClient]);

	return result;
}
```

- [ ] **Step 2: Rewrite `useHasFilesystemProviders`**

Replace `frontend/src/hooks/useHasFilesystemProviders.ts` with:
```ts
import { useFilesystemProviders } from "./useFilesystemProviders";

/**
 * Returns ``true`` when the current room has at least one provider of
 * category ``"filesystem"``.
 */
export function useHasFilesystemProviders(): boolean {
	const { data } = useFilesystemProviders();
	return (data?.length ?? 0) > 0;
}
```

- [ ] **Step 3: Migrate `FilesystemPanel.tsx`**

In `frontend/src/panels/FilesystemPanel.tsx:72-82`, replace the inline `useQuery` with the shared hook:
```ts
import { useFilesystemProviders } from "../hooks/useFilesystemProviders";
// ...
const {
	data: providers,
	isLoading: isLoadingProviders,
	error: providersError,
} = useFilesystemProviders();
```

Remove the `useQuery` import if it is no longer used elsewhere in the file.

- [ ] **Step 4: Type-check**

Run: `cd frontend && bunx tsc --noEmit 2>&1 | tail -10`
Expected: no errors.

- [ ] **Step 5: Run tests**

Run: `cd frontend && bun run test`
Expected: all tests pass.

- [ ] **Step 6: Commit**

```bash
git add frontend/src/hooks/useFilesystemProviders.ts frontend/src/hooks/useHasFilesystemProviders.ts frontend/src/panels/FilesystemPanel.tsx
git commit -m "refactor(hooks): shared useFilesystemProviders hook

Two consumers (useHasFilesystemProviders, FilesystemPanel) used the
same queryKey with different staleTime. Per-observer staleTime means
the Panel's default-0 defeated the Hook's 5_000 dedup whenever the
Panel was open.

Extract one hook owning both query + socket invalidation.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task F6: `useDragHover` hook + shared `shimmer` constant + always-render sliver (Spec 3.4 + 3.5 + 3.12)

Merges three related changes: the shared hook, the shared shimmer constant, and the sliver rendering (which depends on the hook).

**Files:**
- Create: `frontend/src/panels/dragStyles.ts`
- Create: `frontend/src/panels/useDragHover.ts`
- Modify: `frontend/src/panels/ActivityBar.tsx`
- Modify: `frontend/src/panels/SidebarZone.tsx`
- Modify: `frontend/src/panels/BottomZone.tsx`

- [ ] **Step 1: Create the shared shimmer constant**

Write `frontend/src/panels/dragStyles.ts`:
```ts
import { keyframes } from "@mui/system";

export const shimmer = keyframes`
	0%   { background-color: rgba(25, 118, 210, 0.14); }
	50%  { background-color: rgba(25, 118, 210, 0.24); }
	100% { background-color: rgba(25, 118, 210, 0.14); }
`;
```

(Check one of the three existing `keyframes shimmer` declarations first — copy the exact color/stop values so there is zero visual change.)

- [ ] **Step 2: Create the drag-hover hook**

Write `frontend/src/panels/useDragHover.ts`:
```ts
import { useCallback, useRef } from "react";
import { useAppStore } from "../store";

export type BarPosition = "left" | "right" | "bottom";

/**
 * Drag-hover state machine for an ActivityBar / SidebarZone / BottomZone.
 *
 * Tracks `dragDepth` to reliably detect leaving the subtree (nested
 * children generate sibling enter/leave pairs); uses setTimeout(0) on
 * the zero crossing so an immediate re-enter on an adjacent child
 * re-arms the counter before we clear ``hoverBar``.
 */
export function useDragHover(position: BarPosition) {
	const dragDepth = useRef(0);
	const hoverBar = useAppStore((s) => s.hoverBar);
	const setHoverBar = useAppStore((s) => s.setHoverBar);
	const isHovered = hoverBar === position;

	const onDragEnter = useCallback(
		(e: React.DragEvent) => {
			e.preventDefault();
			dragDepth.current++;
			if (dragDepth.current === 1) setHoverBar(position);
		},
		[position, setHoverBar],
	);

	const onDragLeave = useCallback(() => {
		dragDepth.current = Math.max(0, dragDepth.current - 1);
		if (dragDepth.current === 0 && hoverBar === position) {
			setTimeout(() => {
				if (dragDepth.current === 0) setHoverBar(null);
			}, 0);
		}
	}, [position, hoverBar, setHoverBar]);

	const onDragOver = useCallback((e: React.DragEvent) => {
		e.preventDefault();
	}, []);

	return { isHovered, dragHandlers: { onDragEnter, onDragLeave, onDragOver } };
}
```

- [ ] **Step 3: Write a hook unit test**

Write `frontend/src/panels/__tests__/useDragHover.test.tsx`:
```tsx
import { act, renderHook } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { useAppStore } from "../../store";
import { useDragHover } from "../useDragHover";

describe("useDragHover", () => {
	it("sets hoverBar on enter, clears on leave", () => {
		const { result } = renderHook(() => useDragHover("left"));
		expect(result.current.isHovered).toBe(false);

		act(() => {
			result.current.dragHandlers.onDragEnter({ preventDefault() {} } as React.DragEvent);
		});
		expect(useAppStore.getState().hoverBar).toBe("left");
	});

	it("depth counter handles nested enter/leave", async () => {
		const { result } = renderHook(() => useDragHover("right"));
		act(() => {
			result.current.dragHandlers.onDragEnter({ preventDefault() {} } as React.DragEvent);
			result.current.dragHandlers.onDragEnter({ preventDefault() {} } as React.DragEvent);
			result.current.dragHandlers.onDragLeave();
		});
		expect(useAppStore.getState().hoverBar).toBe("right");

		act(() => {
			result.current.dragHandlers.onDragLeave();
		});
		// setTimeout(0) — flush microtasks
		await new Promise((r) => setTimeout(r, 5));
		expect(useAppStore.getState().hoverBar).toBeNull();
	});
});
```

- [ ] **Step 4: Run the hook test — confirm it passes**

Run: `cd frontend && bun run test src/panels/__tests__/useDragHover.test.tsx`
Expected: 2 passed.

- [ ] **Step 5: Migrate `ActivityBar.tsx`**

Edit `frontend/src/panels/ActivityBar.tsx`:

- Remove the module-top `keyframes shimmer` declaration.
- Remove the local `dragDepth`, `onDragEnter`, `onDragOver`, `onDragLeave`, and hover-bar reads (lines 26-134).
- Replace with `const { isHovered, dragHandlers } = useDragHover(position);`.
- Remove the `if (visibleIcons.length === 0 && !isDragActive) return null;` (line 137).
- In the `Box` JSX, spread `{...dragHandlers}` (which takes `onDragEnter`, `onDragLeave`, `onDragOver`).
- Add width-transition sx: when empty+idle, render at 4px; otherwise at the normal bar width. Example:
```ts
const barWidth = visibleIcons.length === 0 && !isDragActive
	? 4  // sliver
	: 48; // normal

// In sx:
{
	width: `${barWidth}px`,
	transition: "width 120ms ease, background-color 120ms ease",
	...dragBgSx,
	...BAR_SX[position],
}
```

Imports:
```ts
import { shimmer } from "./dragStyles";
import { useDragHover } from "./useDragHover";
```

- [ ] **Step 6: Migrate `SidebarZone.tsx`**

Same treatment: remove local drag handlers, import hook+shimmer, call `useDragHover(position)`, spread `dragHandlers`. Apply sliver width logic when appropriate.

- [ ] **Step 7: Migrate `BottomZone.tsx`**

Same treatment (position is `"bottom"`). For vertical sliver, use `height: 4px` when empty+idle.

- [ ] **Step 8: Visually verify there are no remaining `keyframes shimmer` declarations**

Run: `grep -n "keyframes\`" frontend/src/panels/ActivityBar.tsx frontend/src/panels/SidebarZone.tsx frontend/src/panels/BottomZone.tsx`
Expected: no output (all three declarations removed).

Run: `grep -rn "from \"./dragStyles\"" frontend/src/panels/`
Expected: three imports (one per bar file).

- [ ] **Step 9: Type-check and test**

Run: `cd frontend && bunx tsc --noEmit 2>&1 | tail -10`
Expected: no errors.

Run: `cd frontend && bun run test`
Expected: all tests pass.

- [ ] **Step 10: Playwright smoke on dev server**

Per CLAUDE.md, start `uv run zndraw` and `cd frontend && bun run dev`, navigate to a room, open a browser to the dev URL. Manually verify:
- Empty activity bar (drag all icons off one side) stays visible as a narrow strip.
- Dragging an icon highlights the correct bar.
- Dropping adds the icon.

Take screenshots and inspect them.

- [ ] **Step 11: Commit**

```bash
git add frontend/src/panels/dragStyles.ts frontend/src/panels/useDragHover.ts frontend/src/panels/__tests__/useDragHover.test.tsx frontend/src/panels/ActivityBar.tsx frontend/src/panels/SidebarZone.tsx frontend/src/panels/BottomZone.tsx
git commit -m "refactor(panels): useDragHover hook + sliver rendering

Extract ~90 LOC of duplicated drag-depth/setTimeout(0)/onDragEnter
logic into useDragHover(position). Shared keyframes shimmer in
dragStyles.ts (three copies removed).

Empty bars now render as a 4px sliver instead of returning null, so
dragging an icon off a bar still leaves a drop target to drag one
back.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task F7: `activityBarSlice` unit tests (Spec 3.10)

**Files:**
- Create: `frontend/src/stores/slices/__tests__/activityBarSlice.test.ts`

- [ ] **Step 1: Write the test file**

Create `frontend/src/stores/slices/__tests__/activityBarSlice.test.ts`:
```ts
import { beforeEach, describe, expect, it } from "vitest";
import { useAppStore } from "../../../store";

describe("activityBarSlice", () => {
	beforeEach(() => {
		useAppStore.getState().resetLayout();
	});

	it("moveIconToBar moves an icon from one bar to another", () => {
		const s = useAppStore.getState();
		// Assume 'filesystem' defaults to 'left'. Adjust if wrong after
		// reading activityBarSlice.ts initialState().
		expect(s.leftIcons).toContain("filesystem");
		s.moveIconToBar("filesystem", "right");
		expect(useAppStore.getState().leftIcons).not.toContain("filesystem");
		expect(useAppStore.getState().rightIcons).toContain("filesystem");
	});

	it("moveIconToBar with overIdx inserts at position", () => {
		const s = useAppStore.getState();
		s.moveIconToBar("filesystem", "right", 0);
		expect(useAppStore.getState().rightIcons[0]).toBe("filesystem");
	});

	it("dropIconOnPanel removes from bar", () => {
		const s = useAppStore.getState();
		const initial = [...s.leftIcons];
		s.dropIconOnPanel("filesystem");
		expect(useAppStore.getState().leftIcons).not.toContain("filesystem");
		expect(useAppStore.getState().leftIcons.length).toBe(initial.length - 1);
	});

	it("resetLayout restores initial state", () => {
		const s = useAppStore.getState();
		const initialLeft = [...s.leftIcons];
		s.dropIconOnPanel("filesystem");
		s.resetLayout();
		expect(useAppStore.getState().leftIcons).toEqual(initialLeft);
	});

	it("toggleActive opens panel if closed", () => {
		const s = useAppStore.getState();
		expect(s.activeLeft).toBeNull();
		s.toggleActive("left", "filesystem");
		expect(useAppStore.getState().activeLeft).toBe("filesystem");
	});

	it("toggleActive closes panel if already open", () => {
		const s = useAppStore.getState();
		s.toggleActive("left", "filesystem");
		s.toggleActive("left", "filesystem");
		expect(useAppStore.getState().activeLeft).toBeNull();
	});

	it("toggleActive switches active panel", () => {
		const s = useAppStore.getState();
		s.toggleActive("left", "filesystem");
		s.toggleActive("left", "chat");
		expect(useAppStore.getState().activeLeft).toBe("chat");
	});
});
```

Before committing, **read `frontend/src/stores/slices/activityBarSlice.ts`** to verify the exact field names (`leftIcons` vs `icons.left`, `activeLeft` vs `active.left`, etc.) and the default membership of `filesystem`. Adjust the test assertions to match actual shapes — the test file above uses reasonable guesses.

- [ ] **Step 2: Run the tests**

Run: `cd frontend && bun run test src/stores/slices/__tests__/activityBarSlice.test.ts`
Expected: 7 passed.

- [ ] **Step 3: Commit**

```bash
git add frontend/src/stores/slices/__tests__/activityBarSlice.test.ts
git commit -m "test(store): activityBarSlice unit tests

Covers moveIconToBar (with/without overIdx), dropIconOnPanel,
resetLayout, and toggleActive (open/close/switch).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task F8: Selector refactor in `landingPage.tsx:240-249` (Spec 3.11)

**Files:**
- Modify: `frontend/src/pages/landingPage.tsx:240-249`

- [ ] **Step 1: Replace `getState()` block with selectors**

Edit `frontend/src/pages/landingPage.tsx`. The current effect (lines 235-250):
```ts
useEffect(() => {
	const panel = searchParams.get("panel") as PanelId | null;
	if (!panel || !(panel in PANELS)) return;
	const def = PANELS[panel];
	if (def.kind !== "tool" || def.default.bar === "editor") return;
	const state = useAppStore.getState();
	const activeKey =
		def.default.bar === "left"
			? "activeLeft"
			: def.default.bar === "right"
				? "activeRight"
				: "activeBottom";
	if (state[activeKey] !== panel) {
		state.toggleActive(def.default.bar, panel);
	}
}, [searchParams]);
```

New — add the selectors at the top of the component (near other `useAppStore` calls), then rewrite the effect:
```ts
const activeLeft = useAppStore((s) => s.activeLeft);
const activeRight = useAppStore((s) => s.activeRight);
const activeBottom = useAppStore((s) => s.activeBottom);
const toggleActive = useAppStore((s) => s.toggleActive);

useEffect(() => {
	const panel = searchParams.get("panel") as PanelId | null;
	if (!panel || !(panel in PANELS)) return;
	const def = PANELS[panel];
	if (def.kind !== "tool" || def.default.bar === "editor") return;
	const current =
		def.default.bar === "left"
			? activeLeft
			: def.default.bar === "right"
				? activeRight
				: activeBottom;
	if (current !== panel) {
		toggleActive(def.default.bar, panel);
	}
}, [searchParams, activeLeft, activeRight, activeBottom, toggleActive]);
```

- [ ] **Step 2: Type-check and test**

Run: `cd frontend && bunx tsc --noEmit 2>&1 | tail -5`
Expected: no errors.

Run: `cd frontend && bun run test`
Expected: all tests pass.

- [ ] **Step 3: Commit**

```bash
git add frontend/src/pages/landingPage.tsx
git commit -m "refactor(landing): selector-based panel deep-link effect

Replace useAppStore.getState() inside the effect with explicit
selectors in the dep array. Effect now re-runs correctly when the
active bar state changes out-of-band.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task F9: Biome-clean `src/panels` (Spec 3.9)

Must run last so the biome check covers all the other Task F* changes.

**Files:**
- All files under `frontend/src/panels/`

- [ ] **Step 1: Safe autofix pass**

Run: `cd frontend && bunx @biomejs/biome check --write src/panels`
Expected: output reports fixes applied to imports, useOptionalChain, useTemplate.

- [ ] **Step 2: List remaining diagnostics**

Run: `cd frontend && bunx @biomejs/biome check src/panels 2>&1 | tail -40`
Note each remaining error and warning — they need hand fixes.

- [ ] **Step 3: Fix `ChatPanel.tsx` exhaustive-deps (lines 207, 221)**

Read `frontend/src/panels/ChatPanel.tsx:200-230`. For each `useEffect` flagged by `useExhaustiveDependencies`, add the missing deps to the dep array. After adding, verify no infinite re-render by running the frontend dev server and watching the panel.

- [ ] **Step 4: Fix `FilesystemPanel.tsx:79 and 116` noNonNullAssertion**

Read the lines. Replace `roomId!` with an early return:
```ts
// Near the top of the component, after hooks but before return:
if (!roomId) return null;

// Then inside, use `roomId` directly (no `!`).
```

The query already gates on `enabled: !!roomId`, so runtime is safe. The `!` was only there for type coercion — the early return replaces it cleanly.

- [ ] **Step 5: Re-run biome**

Run: `cd frontend && bunx @biomejs/biome check src/panels 2>&1 | tail -20`
Expected: zero errors. (Warnings on `PlotView.tsx` `any` types remain — out of scope per spec §3.9.)

- [ ] **Step 6: Run all tests**

Run: `cd frontend && bun run test`
Expected: all tests pass.

- [ ] **Step 7: Commit**

```bash
git add frontend/src/panels/
git commit -m "chore(panels): biome-clean

Safe autofix pass + hand fixes for noNonNullAssertion (FilesystemPanel
early-returns on null roomId) and useExhaustiveDependencies (ChatPanel).

PlotView.tsx 'any' types are pre-existing legacy debt and explicitly
out of scope per spec §3.9.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Final Verification

- [ ] **Step 1: Backend — full suite**

Run: `uv run pytest tests/zndraw tests/zndraw_joblib -q`
Expected: all tests pass.

- [ ] **Step 2: Frontend — all tests, all type-checks**

Run: `cd frontend && bun run test && bunx tsc --noEmit && bunx @biomejs/biome check src/panels src/hooks src/stores`
Expected: all green.

- [ ] **Step 3: prek pre-commit hooks**

Run: `uvx prek --all-files 2>&1 | tail -20`
Expected: all green.

- [ ] **Step 4: Manual smoke (per CLAUDE.md, dev server)**

Start backend: `uv run zndraw`
Start frontend: `cd frontend && bun run dev`

In a browser (and/or via playwright-cli):
- Navigate to a room. Empty activity bar renders as a sliver.
- Drag a sidebar icon — drop target is highlighted on hover.
- Use the Reset-layout MenuItem — all panels close, viewer reappears.
- Filesystem icon: present as superuser, gated for non-superusers.
- Direct-URL `?panel=filesystem` — opens the panel once the dockview API binds (no polling).

Take screenshots, compare to anticipated behavior.

- [ ] **Step 5: PR update**

Push the branch. On the PR, reply to the code-review comments marking each resolved item, linking to the commit that addressed it.

---

## Out of Scope (Follow-ups)

1. Wire frontend e2e into CI. Infrastructure change — separate PR.
2. `PlotView.tsx` 14 `any` types. Legacy from `FigureWindow`, pre-existing before PR #920.
3. Review finding #11 (`_DummyProvider.category`): confirmed misread; no action.
