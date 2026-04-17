# Internal Providers & Default Filesystem Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an `@internal` Provider path — symmetric to `@internal` extensions — that registers `Provider` subclasses on the taskiq broker, seeds `ProviderRecord` DB rows, and dispatches reads via `taskiq.kiq` instead of Socket.IO. First consumer: a default `FilesystemRead` provider rooted at `ZNDRAW_SERVER_FILEBROWSER_PATH` (default `"."`, disable with `"none"`).

**Architecture:** Mirror the `@internal` extension flow. `zndraw_joblib.registry` gains `register_internal_providers()` + `ensure_internal_providers()`. `read_provider()` in `zndraw_joblib.router` forks on `provider.room_id == "@internal"` — taskiq dispatch instead of `tsio.emit(ProviderRequest, ...)`. A new `InternalProviderExecutor` in `zndraw.providers.executor` resolves the fsspec handler from `settings.filebrowser_path`, runs `provider_cls(**params).read(handler)` in a thread, and POSTs the bytes back to `/v1/joblib/providers/{id}/results` — the same endpoint remote providers use.

**Tech Stack:** Python 3.12, FastAPI, SQLModel, taskiq / taskiq-redis, fsspec, httpx (sync, in `asyncio.to_thread`), Pydantic v2 + pydantic-settings. Frontend: TypeScript/React 19, MUI, TanStack Query, zustand. `uv` for Python, `bun` for frontend.

**Reference spec:** [`docs/superpowers/specs/2026-04-17-internal-providers-design.md`](../specs/2026-04-17-internal-providers-design.md).

---

## File Structure

| File | Action | Responsibility |
|---|---|---|
| `src/zndraw/config.py` | Modify | Add `filebrowser_path: str = "."` (case-insensitive `"none"` disables). |
| `src/zndraw/providers/__init__.py` | Modify | Export a `BUNDLED_PROVIDERS` list of `Provider` subclasses shipped with the server. |
| `src/zndraw/providers/executor.py` | Create | `InternalProviderExecutor` — resolves filesystem handler, runs `Provider.read()`, POSTs result. |
| `src/zndraw_joblib/registry.py` | Modify | Add `register_internal_providers()` and `ensure_internal_providers()`. |
| `src/zndraw_joblib/__init__.py` | Modify | Re-export the two new helpers + `InternalProviderExecutor` protocol. |
| `src/zndraw_joblib/router.py` | Modify | Widen `_room_provider_filter`, `_resolve_provider` visibility; fork `read_provider` on `@internal`; refuse `delete_provider` for `@internal`. |
| `src/zndraw/database.py` | Modify | Seed internal worker `Worker` row; wire `register_internal_providers` + `ensure_internal_providers` in lifespan; add `_collect_providers()`. |
| `src/zndraw/broker.py` | Modify | External taskiq worker also registers internal providers. |
| `frontend/src/panels/ActivityBar.tsx` | Modify | Hide the `filesystem` icon when the provider list for the current room has no filesystem provider. |
| `frontend/src/hooks/useHasFilesystemProviders.ts` | Create | TanStack-Query hook — returns `boolean` for filesystem-category providers in the active room. |
| `tests/zndraw_joblib/test_providers.py` | Modify (append) | Unit: `@internal` widening in `_room_provider_filter`; cache-hit path through internal provider row. |
| `tests/zndraw_joblib/test_registry.py` | Modify (append) | Unit: `register_internal_providers` registers tasks; `ensure_internal_providers` seeds idempotent rows. |
| `tests/zndraw/test_providers_filesystem.py` | Modify (append) | Integration: default server exposes `@internal:filesystem:FilesystemRead`; `filebrowser_path="none"` yields empty provider list. |
| `tests/zndraw/test_config.py` | Modify (append) | Unit: `filebrowser_path` default + env override + `"none"` sentinel. |

No file deletions. No new Python dependencies (fsspec is already a dep via `FilesystemRead`). No new npm dependencies.

---

## Preflight

### Task 0: Verify baseline

Establish a clean starting point so later failures can be attributed to this plan.

**Files:** none.

- [ ] **Step 0.1: Confirm branch and working tree.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git status --short && git rev-parse --abbrev-ref HEAD
```

Expected: HEAD is `spec/dockview-ui-redesign` (or a new branch spun off it). Only expected dirty paths: `CLAUDE.md`, `ISSUES.md`, `docs/...`, `frontend/.playwright-cli/`, `frontend/playwright-report/`, `frontend/test-results/`, `room-default-snapshot.yaml`, modified screenshots, `skills-lock.json`.

- [ ] **Step 0.2: Backend typecheck + lint.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uvx prek run --all-files
```

Expected: zero errors. (If it reports pre-existing failures, record them — anything you *add* must not grow that list.)

- [ ] **Step 0.3: Run the current @internal extension suite to confirm the baseline @internal wiring works.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw/worker/test_internal.py tests/zndraw_joblib/test_providers.py tests/zndraw_joblib/test_registry.py tests/zndraw/test_providers_filesystem.py -q
```

Expected: all tests pass (modulo any pre-existing skips/xfails you noticed in Step 0.2).

- [ ] **Step 0.4: Frontend build sanity.**

```bash
cd /Users/fzills/tools/zndraw-fastapi/frontend && bun install && bun --bun tsc --noEmit
```

Expected: zero errors.

---

## Phase 1: Settings

### Task 1: Add `filebrowser_path` setting

**Files:**
- Modify: `src/zndraw/config.py` (add field after `internal_url`)
- Test: `tests/zndraw/test_config.py` (append a `TestFilebrowserPath` class)

- [ ] **Step 1.1: Write the failing test.**

Append to `/Users/fzills/tools/zndraw-fastapi/tests/zndraw/test_config.py`:

```python
class TestFilebrowserPath:
    """Test filebrowser_path configuration."""

    def test_default_filebrowser_path_is_cwd(self) -> None:
        """Default filebrowser_path should be '.'."""
        settings = Settings()
        assert settings.filebrowser_path == "."

    def test_filebrowser_path_from_env(self) -> None:
        """filebrowser_path should be configurable via ZNDRAW_SERVER_FILEBROWSER_PATH."""
        os.environ["ZNDRAW_SERVER_FILEBROWSER_PATH"] = "/data"
        try:
            settings = Settings()
            assert settings.filebrowser_path == "/data"
        finally:
            os.environ.pop("ZNDRAW_SERVER_FILEBROWSER_PATH", None)

    def test_filebrowser_path_none_sentinel(self) -> None:
        """Sentinel 'none' (case-insensitive) disables the default provider."""
        os.environ["ZNDRAW_SERVER_FILEBROWSER_PATH"] = "NONE"
        try:
            settings = Settings()
            assert settings.filebrowser_path == "NONE"
            assert settings.filebrowser_path.lower() == "none"
        finally:
            os.environ.pop("ZNDRAW_SERVER_FILEBROWSER_PATH", None)
```

- [ ] **Step 1.2: Run the test — it must fail.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw/test_config.py::TestFilebrowserPath -v
```

Expected: all three tests FAIL with `AttributeError: 'Settings' object has no attribute 'filebrowser_path'`.

- [ ] **Step 1.3: Add the field.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw/config.py`, after the `internal_url` line (currently last field in class), add:

```python
    # Filesystem provider
    filebrowser_path: str = "."
    """Path for the default @internal filesystem provider.

    ``"."`` (default) roots at the process cwd (taskiq-worker's cwd in Docker).
    Any other non-``"none"`` string roots the provider at that absolute or
    relative path. The sentinel ``"none"`` (case-insensitive) disables the
    default filesystem provider entirely — no DB row, no task registration,
    and the frontend filesystem activity-bar icon is hidden.
    """
```

- [ ] **Step 1.4: Run the test — must pass.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw/test_config.py::TestFilebrowserPath -v
```

Expected: 3 passed.

- [ ] **Step 1.5: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add src/zndraw/config.py tests/zndraw/test_config.py && git commit -m "$(cat <<'EOF'
feat(config): add ZNDRAW_SERVER_FILEBROWSER_PATH setting

Introduces the filebrowser_path setting that controls the default
@internal filesystem provider root. The sentinel "none" disables it.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 2: Discovery helpers in `zndraw_joblib.registry`

### Task 2: Introduce `InternalProviderExecutor` protocol and `register_internal_providers`

**Files:**
- Modify: `src/zndraw_joblib/registry.py`
- Test: `tests/zndraw_joblib/test_registry.py`

- [ ] **Step 2.1: Write the failing test.**

Read `/Users/fzills/tools/zndraw-fastapi/tests/zndraw_joblib/test_registry.py` first to match the existing style, then append:

```python
# --- Provider registry ---------------------------------------------------

from typing import ClassVar

from zndraw_joblib.provider import Provider


class _DummyProvider(Provider):
    category: ClassVar[str] = "filesystem"
    path: str = "/"

    def read(self, handler):
        return handler.ls(self.path, detail=True)


async def test_register_internal_providers_registers_taskiq_task():
    """Each provider class becomes a task at @internal:<cat>:<ClassName>."""
    from taskiq import InMemoryBroker

    from zndraw_joblib.registry import register_internal_providers

    broker = InMemoryBroker()

    calls = []

    async def executor(cls, params_json, provider_id, request_id, token):
        calls.append((cls, params_json, provider_id, request_id, token))

    reg = register_internal_providers(broker, [_DummyProvider], executor)

    assert "@internal:filesystem:_DummyProvider" in reg.tasks
    assert reg.providers["@internal:filesystem:_DummyProvider"] is _DummyProvider


async def test_register_internal_providers_task_invokes_executor():
    """Calling the registered task fan-outs to the executor with forwarded args."""
    from taskiq import InMemoryBroker

    from zndraw_joblib.registry import register_internal_providers

    broker = InMemoryBroker()
    calls = []

    async def executor(cls, params_json, provider_id, request_id, token):
        calls.append((cls.__name__, params_json, provider_id, request_id, token))

    reg = register_internal_providers(broker, [_DummyProvider], executor)
    task = reg.tasks["@internal:filesystem:_DummyProvider"]

    # Directly invoke the underlying coroutine (InMemoryBroker doesn't run it)
    await task.original_func(
        request_id="abc",
        provider_id="11111111-1111-1111-1111-111111111111",
        params_json='{"path": "/data"}',
        token="tok",
    )

    assert calls == [
        (
            "_DummyProvider",
            '{"path": "/data"}',
            "11111111-1111-1111-1111-111111111111",
            "abc",
            "tok",
        )
    ]
```

Add the `pytest.mark.anyio` or `pytest.mark.asyncio` marker as appropriate — check the top of `tests/zndraw_joblib/test_registry.py` for the existing convention and mirror it. If no async test runner is already configured in the file, import and decorate with `@pytest.mark.asyncio` (that's what the rest of the suite uses — see `test_providers.py::test_read_provider_long_poll_wakes_on_upload`).

- [ ] **Step 2.2: Run the test — must fail.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_registry.py -v -k provider
```

Expected: FAIL — `ImportError: cannot import name 'register_internal_providers' from 'zndraw_joblib.registry'`.

- [ ] **Step 2.3: Implement `register_internal_providers`.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw_joblib/registry.py`, add imports at the top (alongside existing `from zndraw_joblib.client import Extension`):

```python
if TYPE_CHECKING:
    from zndraw_joblib.provider import Provider
```

After `register_internal_tasks` (before `ensure_internal_jobs`), add:

```python
class InternalProviderExecutor(Protocol):
    """Protocol for the host-provided provider executor.

    The server base URL and any handler configuration are captured at
    creation time. Per-request data (params, ids, token) is passed at
    call time.
    """

    async def __call__(
        self,
        provider_cls: type[Provider],
        params_json: str,
        provider_id: str,
        request_id: str,
        token: str,
    ) -> None: ...


@dataclass
class InternalProviderRegistry:
    """Holds taskiq task handles and provider class mappings."""

    tasks: dict[str, Any] = field(default_factory=dict)
    providers: dict[str, type[Provider]] = field(default_factory=dict)
    executor: InternalProviderExecutor | None = None


def register_internal_providers(
    broker: AsyncBroker,
    providers: list[type[Provider]],
    executor: InternalProviderExecutor,
) -> InternalProviderRegistry:
    """Register Provider classes as taskiq tasks on the broker.

    Each class becomes a task named ``@internal:<category>:<ClassName>``.
    The task handler forwards ``(cls, params_json, provider_id, request_id,
    token)`` to *executor*.
    """
    registry = InternalProviderRegistry(executor=executor)

    for prov_cls in providers:
        category = prov_cls.category
        name = prov_cls.__name__
        full_name = f"@internal:{category}:{name}"

        def _make_task_fn(
            cls: type[Provider] = prov_cls,
            ex: InternalProviderExecutor = executor,
        ):
            async def _execute(
                request_id: str,
                provider_id: str,
                params_json: str,
                token: str,
            ) -> None:
                await ex(cls, params_json, provider_id, request_id, token)

            return _execute

        task_handle = broker.register_task(
            _make_task_fn(),
            task_name=full_name,
        )

        registry.tasks[full_name] = task_handle
        registry.providers[full_name] = prov_cls
        logger.debug("Registered internal provider task: %s", full_name)

    logger.info("Registered %d internal provider task(s)", len(providers))
    return registry
```

- [ ] **Step 2.4: Run the test — must pass.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_registry.py -v -k provider
```

Expected: 2 passed.

- [ ] **Step 2.5: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add src/zndraw_joblib/registry.py tests/zndraw_joblib/test_registry.py && git commit -m "$(cat <<'EOF'
feat(joblib): register_internal_providers helper

Registers Provider subclasses as taskiq tasks named
@internal:<category>:<ClassName>, symmetric to
register_internal_tasks for extensions.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: `ensure_internal_providers` seeds DB rows

**Files:**
- Modify: `src/zndraw_joblib/registry.py`
- Test: `tests/zndraw_joblib/test_registry.py`

Internal provider rows need a `user_id` (internal worker user) and `worker_id` (a stable Worker row for that user). Callers must create those *before* calling this helper — the helper itself only upserts `ProviderRecord` rows.

- [ ] **Step 3.1: Write the failing test.**

Append to `/Users/fzills/tools/zndraw-fastapi/tests/zndraw_joblib/test_registry.py`:

```python
async def test_ensure_internal_providers_creates_rows(async_session_factory):
    """Creates a ProviderRecord per provider at room_id='@internal'."""
    import uuid

    from sqlmodel import select

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord, Worker
    from zndraw_joblib.registry import ensure_internal_providers

    # Seed internal user + worker
    user_id = uuid.uuid4()
    async with async_session_factory() as session:
        user = User(
            id=user_id,
            email="internal@test",
            hashed_password="x",
            is_active=True,
            is_superuser=True,
            is_verified=True,
        )
        session.add(user)
        worker = Worker(user_id=user_id)
        session.add(worker)
        await session.commit()
        worker_id = worker.id

    await ensure_internal_providers(
        [_DummyProvider], async_session_factory, user_id=user_id, worker_id=worker_id
    )

    async with async_session_factory() as session:
        result = await session.exec(
            select(ProviderRecord).where(ProviderRecord.room_id == "@internal")
        )
        rows = result.all()
    assert len(rows) == 1
    assert rows[0].category == "filesystem"
    assert rows[0].name == "_DummyProvider"
    assert rows[0].user_id == user_id
    assert rows[0].worker_id == worker_id


async def test_ensure_internal_providers_idempotent(async_session_factory):
    """Running twice leaves exactly one row (upsert on room+category+name)."""
    import uuid

    from sqlmodel import select

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord, Worker
    from zndraw_joblib.registry import ensure_internal_providers

    user_id = uuid.uuid4()
    async with async_session_factory() as session:
        session.add(
            User(
                id=user_id,
                email="internal@test",
                hashed_password="x",
                is_active=True,
                is_superuser=True,
                is_verified=True,
            )
        )
        worker = Worker(user_id=user_id)
        session.add(worker)
        await session.commit()
        worker_id = worker.id

    await ensure_internal_providers(
        [_DummyProvider], async_session_factory, user_id=user_id, worker_id=worker_id
    )
    await ensure_internal_providers(
        [_DummyProvider], async_session_factory, user_id=user_id, worker_id=worker_id
    )

    async with async_session_factory() as session:
        result = await session.exec(
            select(ProviderRecord).where(ProviderRecord.room_id == "@internal")
        )
        rows = result.all()
    assert len(rows) == 1
```

- [ ] **Step 3.2: Run the test — must fail.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_registry.py -v -k ensure_internal_providers
```

Expected: FAIL — `ImportError: cannot import name 'ensure_internal_providers'`.

- [ ] **Step 3.3: Implement `ensure_internal_providers`.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw_joblib/registry.py`, add below `register_internal_providers` (and above `ensure_internal_jobs`):

```python
async def ensure_internal_providers(
    providers: list[type[Provider]],
    session_factory: Callable[[], AbstractAsyncContextManager[AsyncSession]],
    *,
    user_id: UUID,
    worker_id: UUID,
) -> None:
    """Create or update @internal ProviderRecord rows.

    Idempotent — safe to call on every startup. Callers must have already
    seeded a ``User`` (``user_id``) and ``Worker`` (``worker_id``). In
    multi-replica production, call once from ``init_database()`` to avoid
    races on the ``unique_provider`` constraint.
    """
    from sqlmodel import select

    from zndraw_joblib.models import ProviderRecord

    async with session_factory() as session:
        for prov_cls in providers:
            category = prov_cls.category
            name = prov_cls.__name__
            schema = prov_cls.model_json_schema()
            content_type = prov_cls.content_type

            result = await session.exec(
                select(ProviderRecord).where(
                    ProviderRecord.room_id == "@internal",
                    ProviderRecord.category == category,
                    ProviderRecord.name == name,
                )
            )
            existing = result.one_or_none()

            if existing:
                existing.schema_ = schema
                existing.content_type = content_type
                existing.user_id = user_id
                existing.worker_id = worker_id
            else:
                session.add(
                    ProviderRecord(
                        room_id="@internal",
                        category=category,
                        name=name,
                        schema_=schema,
                        content_type=content_type,
                        user_id=user_id,
                        worker_id=worker_id,
                    )
                )

        await session.commit()

    logger.info("Ensured %d @internal provider row(s) in DB", len(providers))
```

Also add the required import at the top of the file (inside the `if TYPE_CHECKING` block, already importing `Callable` / `AbstractAsyncContextManager` / `AsyncSession`):

```python
    from uuid import UUID
```

(Or just import it at module level alongside `from uuid import UUID` — the file already imports from `typing` & `dataclasses`, so one extra top-level import is fine.)

- [ ] **Step 3.4: Run the test — must pass.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_registry.py -v -k ensure_internal_providers
```

Expected: 2 passed.

- [ ] **Step 3.5: Re-export from package.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw_joblib/__init__.py`:

1. Update the `from zndraw_joblib.registry import ...` block to also import:

```python
from zndraw_joblib.registry import (
    InternalExecutor,
    InternalProviderExecutor,
    InternalProviderRegistry,
    InternalRegistry,
    ensure_internal_providers,
    register_internal_jobs,
    register_internal_providers,
    register_internal_tasks,
)
```

2. Add to `__all__` (alphabetical): `"InternalProviderExecutor"`, `"InternalProviderRegistry"`, `"ensure_internal_providers"`, `"register_internal_providers"`.

- [ ] **Step 3.6: Run a broad unit sweep.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib -q
```

Expected: all pass.

- [ ] **Step 3.7: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add src/zndraw_joblib/registry.py src/zndraw_joblib/__init__.py tests/zndraw_joblib/test_registry.py && git commit -m "$(cat <<'EOF'
feat(joblib): ensure_internal_providers seeds ProviderRecord rows

Idempotent upsert of @internal ProviderRecord rows (keyed by
room+category+name). Callers pass the internal worker User and
Worker ids, parallel to ensure_internal_jobs for extensions.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 3: Router — dispatch fork + visibility widening

### Task 4: Widen provider visibility to include `@internal`

**Files:**
- Modify: `src/zndraw_joblib/router.py`
- Test: `tests/zndraw_joblib/test_providers.py`

- [ ] **Step 4.1: Write the failing test.**

Append to `/Users/fzills/tools/zndraw-fastapi/tests/zndraw_joblib/test_providers.py` (after `test_list_providers_mixed_scopes`):

```python
def test_list_providers_includes_internal(client, async_session_factory):
    """Internal providers are visible from every room (and from @global)."""
    import asyncio
    import uuid

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord, Worker

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

    asyncio.run(seed())

    # A normal room sees the @internal provider
    resp = client.get("/v1/joblib/rooms/room-42/providers")
    assert resp.status_code == 200
    items = resp.json()["items"]
    assert any(
        p["full_name"] == "@internal:filesystem:FilesystemRead" for p in items
    )

    # @global also sees only its own scope (not @internal)
    resp = client.get("/v1/joblib/rooms/@global/providers")
    assert all(
        p["full_name"] != "@internal:filesystem:FilesystemRead"
        for p in resp.json()["items"]
    )


def test_get_provider_info_internal_visible_from_room(client, async_session_factory):
    """A normal room can fetch info on an @internal provider."""
    import asyncio
    import uuid

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord, Worker

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
            worker = Worker(user_id=user.id)
            session.add(worker)
            await session.flush()
            session.add(
                ProviderRecord(
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemRead",
                    schema_={"path": {"type": "string"}},
                    user_id=user.id,
                    worker_id=worker.id,
                )
            )
            await session.commit()

    asyncio.run(seed())

    resp = client.get(
        "/v1/joblib/rooms/room-42/providers/"
        "@internal:filesystem:FilesystemRead/info"
    )
    assert resp.status_code == 200
    assert resp.json()["schema"] == {"path": {"type": "string"}}
```

> **Note:** the `async_session_factory` fixture is already defined in `tests/zndraw_joblib/conftest.py`. No new fixtures are needed.

- [ ] **Step 4.2: Run the test — must fail.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_providers.py::test_list_providers_includes_internal tests/zndraw_joblib/test_providers.py::test_get_provider_info_internal_visible_from_room -v
```

Expected: FAIL — list query returns 0 `@internal` rows; info query returns 404.

- [ ] **Step 4.3: Widen `_room_provider_filter`.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw_joblib/router.py`, replace the existing function body (lines 970-974):

```python
def _room_provider_filter(room_id: str):
    """Build a SQLAlchemy filter for providers visible from a given room."""
    if room_id == "@global":
        return ProviderRecord.room_id == "@global"
    if room_id == "@internal":
        return ProviderRecord.room_id == "@internal"
    return ProviderRecord.room_id.in_(["@global", "@internal", room_id])
```

- [ ] **Step 4.4: Widen `_resolve_provider` visibility check.**

In the same file, `_resolve_provider` currently (around line 989) does:

```python
    if room_id not in ("@global",) and provider_room_id not in (
        "@global",
        room_id,
    ):
```

Replace with:

```python
    if room_id not in ("@global",) and provider_room_id not in (
        "@global",
        "@internal",
        room_id,
    ):
```

- [ ] **Step 4.5: Run the tests — must pass.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_providers.py -v
```

Expected: all pass (including the two new tests).

- [ ] **Step 4.6: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add src/zndraw_joblib/router.py tests/zndraw_joblib/test_providers.py && git commit -m "$(cat <<'EOF'
feat(joblib): widen provider visibility to include @internal

Normal rooms now see providers scoped to @global, @internal, and
their own room_id. Mirrors the existing job visibility rules.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Refuse deletion of `@internal` providers

**Files:**
- Modify: `src/zndraw_joblib/router.py`
- Test: `tests/zndraw_joblib/test_providers.py`

- [ ] **Step 5.1: Write the failing test.**

Append to `/Users/fzills/tools/zndraw-fastapi/tests/zndraw_joblib/test_providers.py`:

```python
def test_delete_internal_provider_forbidden(client_factory, async_session_factory):
    """@internal providers cannot be deleted by anyone, including superusers."""
    import asyncio
    import uuid

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord, Worker

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
            worker = Worker(user_id=user.id)
            session.add(worker)
            await session.flush()
            session.add(
                ProviderRecord(
                    id=provider_id,
                    room_id="@internal",
                    category="filesystem",
                    name="FilesystemRead",
                    schema_={},
                    user_id=user.id,
                    worker_id=worker.id,
                )
            )
            await session.commit()

    asyncio.run(seed())

    admin = client_factory("admin", is_superuser=True)
    resp = admin.delete(f"/v1/joblib/providers/{provider_id}")
    assert resp.status_code == 403
```

- [ ] **Step 5.2: Run the test — must fail.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_providers.py::test_delete_internal_provider_forbidden -v
```

Expected: FAIL — status is 204 (superuser can currently delete anything).

- [ ] **Step 5.3: Refuse `@internal` deletion.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw_joblib/router.py`, modify `delete_provider` (around line 1218). After the `ProviderNotFound` check and before the user-ownership check, add:

```python
    if provider.room_id == "@internal":
        raise Forbidden.exception(
            detail="@internal providers cannot be deleted"
        )
```

- [ ] **Step 5.4: Run the test — must pass.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_providers.py -v
```

Expected: all pass.

- [ ] **Step 5.5: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add src/zndraw_joblib/router.py tests/zndraw_joblib/test_providers.py && git commit -m "$(cat <<'EOF'
feat(joblib): refuse deletion of @internal providers

@internal providers are owned by the server, not users. Parallel
to the extension side where @internal Job rows are immutable.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Dispatch fork — `read_provider` → taskiq for `@internal`

**Files:**
- Modify: `src/zndraw_joblib/router.py`
- Modify: `src/zndraw_joblib/dependencies.py` (expose a provider registry dependency — *only if it doesn't already exist*; it doesn't, so add it)
- Test: `tests/zndraw_joblib/test_providers.py`

- [ ] **Step 6.1: Add `get_internal_provider_registry` dependency.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw_joblib/dependencies.py`, next to `get_internal_registry`, add:

```python
from zndraw_joblib.registry import InternalProviderRegistry, InternalRegistry
```

(replace the existing single-class import) and:

```python
async def get_internal_provider_registry(
    request: Request,
) -> InternalProviderRegistry | None:
    """Return the internal provider registry from app.state, or None."""
    return getattr(request.app.state, "internal_provider_registry", None)
```

- [ ] **Step 6.2: Write the failing test.**

Append to `/Users/fzills/tools/zndraw-fastapi/tests/zndraw_joblib/test_providers.py`:

```python
def test_read_internal_provider_dispatches_via_taskiq(
    client, app, async_session_factory, monkeypatch
):
    """An @internal provider dispatches via the registry's taskiq task, not tsio."""
    import asyncio
    import uuid

    from zndraw_auth import User
    from zndraw_joblib.models import ProviderRecord, Worker
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
    assert call["params_json"] == '{"path":"/data"}'
    assert "request_id" in call
    assert "provider_id" in call
    assert call["token"] == "test-worker-token"
```

- [ ] **Step 6.3: Run the test — must fail.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_providers.py::test_read_internal_provider_dispatches_via_taskiq -v
```

Expected: FAIL — current code emits via tsio; no kiq call is made.

- [ ] **Step 6.4: Add the dispatch fork.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw_joblib/router.py`:

1. Add imports (top of file, with the other imports):

```python
import json as _json

from zndraw_joblib.dependencies import (
    ...,
    get_internal_provider_registry,
)
from zndraw_joblib.registry import InternalProviderRegistry, InternalRegistry
```

(Merge these into the existing import blocks; `InternalRegistry` is already imported.)

2. Add a type alias near the other DI aliases (~line 98):

```python
InternalProviderRegistryDep = Annotated[
    InternalProviderRegistry | None, Depends(get_internal_provider_registry)
]
```

3. Modify `read_provider`'s signature to add two dependencies:

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
    worker_token: WorkerTokenDep,
    prefer: Annotated[str | None, Header()] = None,
):
```

4. Replace the dispatch branch (currently `if acquired: await emit(tsio, ...)`) with:

```python
    if acquired:
        if provider.room_id == "@internal":
            if (
                internal_provider_registry is None
                or provider.full_name not in internal_provider_registry.tasks
            ):
                raise InternalJobNotConfigured.exception(
                    detail=(
                        f"Internal provider '{provider.full_name}' is registered"
                        " in the DB but no executor task is available"
                    )
                )
            params_json = _json.dumps(params, sort_keys=True, separators=(",", ":"))
            await internal_provider_registry.tasks[provider.full_name].kiq(
                request_id=rhash,
                provider_id=str(provider.id),
                params_json=params_json,
                token=worker_token,
            )
        else:
            provider_room = f"providers:{provider.full_name}"
            await emit(
                tsio,
                {
                    Emission(
                        ProviderRequest.from_dict_params(
                            request_id=rhash,
                            provider_name=provider.full_name,
                            params=params,
                        ),
                        provider_room,
                    )
                },
            )
```

- [ ] **Step 6.5: Run the tests — must pass.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib/test_providers.py -v
```

Expected: all pass (existing remote-dispatch tests still pass since the `else` branch is unchanged).

- [ ] **Step 6.6: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add src/zndraw_joblib/router.py src/zndraw_joblib/dependencies.py tests/zndraw_joblib/test_providers.py && git commit -m "$(cat <<'EOF'
feat(joblib): dispatch @internal providers via taskiq

read_provider forks on provider.room_id == "@internal" and enqueues
the registered taskiq task instead of emitting ProviderRequest over
Socket.IO. Existing remote-provider path is unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 4: Host app — executor + discovery + lifespan wiring

### Task 7: `InternalProviderExecutor` implementation

**Files:**
- Create: `src/zndraw/providers/executor.py`
- Modify: `src/zndraw/providers/__init__.py`
- Test: new unit test `tests/zndraw/test_providers_executor.py`

- [ ] **Step 7.1: Write the failing test.**

Create `/Users/fzills/tools/zndraw-fastapi/tests/zndraw/test_providers_executor.py`:

```python
"""Unit tests for InternalProviderExecutor."""

import json
from pathlib import Path

import pytest
from httpx import MockTransport, Response

from zndraw.providers.executor import InternalProviderExecutor
from zndraw.providers.filesystem import FilesystemRead


@pytest.fixture
def seeded_dir(tmp_path: Path) -> Path:
    (tmp_path / "a.xyz").touch()
    (tmp_path / "b.xyz").touch()
    return tmp_path


async def test_executor_posts_result_for_filesystem(seeded_dir, monkeypatch):
    """Executor resolves fsspec handler, runs read(), POSTs bytes."""
    captured: dict = {}

    def handler(request):
        captured["url"] = str(request.url)
        captured["headers"] = dict(request.headers)
        captured["body"] = request.content
        return Response(204)

    transport = MockTransport(handler)

    executor = InternalProviderExecutor(
        base_url="http://test",
        filebrowser_path=str(seeded_dir),
        _transport=transport,  # test seam
    )

    await executor(
        FilesystemRead,
        params_json=json.dumps({"path": "/"}),
        provider_id="11111111-1111-1111-1111-111111111111",
        request_id="abc",
        token="tok",
    )

    assert captured["url"].endswith(
        "/v1/joblib/providers/11111111-1111-1111-1111-111111111111/results"
    )
    assert captured["headers"]["authorization"] == "Bearer tok"
    assert captured["headers"]["x-request-hash"] == "abc"
    body = json.loads(captured["body"])
    names = {item["name"] for item in body}
    assert names == {"a.xyz", "b.xyz"}
```

- [ ] **Step 7.2: Run the test — must fail.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw/test_providers_executor.py -v
```

Expected: FAIL — module doesn't exist yet.

- [ ] **Step 7.3: Implement `InternalProviderExecutor`.**

Create `/Users/fzills/tools/zndraw-fastapi/src/zndraw/providers/executor.py`:

```python
"""Taskiq-side executor for @internal providers.

Resolves the filesystem handler from configured settings, invokes
``Provider.read(handler)``, and POSTs the result to the server via
the same ``/v1/joblib/providers/{id}/results`` endpoint that remote
providers use.
"""

from __future__ import annotations

import asyncio
import json
import logging
from dataclasses import dataclass, field
from typing import Any

import httpx

from zndraw_joblib.provider import Provider

log = logging.getLogger(__name__)


@dataclass
class InternalProviderExecutor:
    """Execute an @internal Provider read and POST the result.

    Parameters
    ----------
    base_url
        ZnDraw server URL (e.g. ``http://127.0.0.1:8000``).
    filebrowser_path
        Absolute or relative path rooting the ``filesystem`` handler.
        Relative paths are resolved against the taskiq-worker's cwd
        when the executor is constructed.
    _transport
        Optional httpx transport for testing. Not set in production.
    """

    base_url: str
    filebrowser_path: str
    _transport: Any = field(default=None, repr=False)

    async def __call__(
        self,
        provider_cls: type[Provider],
        params_json: str,
        provider_id: str,
        request_id: str,
        token: str,
    ) -> None:
        """Execute the provider and POST the result. Raises on upload failure."""
        base_url = self.base_url
        transport = self._transport

        def _run() -> None:
            handler = self._resolve_handler(provider_cls)
            params = json.loads(params_json) if params_json else {}
            instance = provider_cls(**params)
            result = instance.read(handler)

            if provider_cls.content_type == "application/json":
                content = json.dumps(result).encode()
            else:
                content = result  # type: ignore[assignment]

            client_kwargs: dict[str, Any] = {"timeout": 30.0}
            if transport is not None:
                client_kwargs["transport"] = transport

            with httpx.Client(**client_kwargs) as client:
                resp = client.post(
                    f"{base_url}/v1/joblib/providers/{provider_id}/results",
                    content=content,
                    headers={
                        "Authorization": f"Bearer {token}",
                        "X-Request-Hash": request_id,
                    },
                )
                resp.raise_for_status()

        await asyncio.to_thread(_run)

    def _resolve_handler(self, provider_cls: type[Provider]) -> Any:
        """Resolve a handler for the provider category."""
        category = provider_cls.category
        if category == "filesystem":
            import fsspec
            from fsspec.implementations.dirfs import DirFileSystem

            return DirFileSystem(
                path=self.filebrowser_path, fs=fsspec.filesystem("file")
            )
        raise ValueError(
            f"No internal handler configured for provider category '{category}'"
        )
```

- [ ] **Step 7.4: Add `BUNDLED_PROVIDERS` to the providers package.**

Overwrite `/Users/fzills/tools/zndraw-fastapi/src/zndraw/providers/__init__.py`:

```python
"""Built-in Provider subclasses shipped with the server.

These are registered at ``@internal`` scope on startup and dispatched
via the taskiq worker (see ``zndraw.providers.executor``).
"""

from zndraw.providers.filesystem import FilesystemRead

BUNDLED_PROVIDERS = [FilesystemRead]

__all__ = ["BUNDLED_PROVIDERS", "FilesystemRead"]
```

- [ ] **Step 7.5: Run the test — must pass.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw/test_providers_executor.py -v
```

Expected: 1 passed.

- [ ] **Step 7.6: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add src/zndraw/providers/executor.py src/zndraw/providers/__init__.py tests/zndraw/test_providers_executor.py && git commit -m "$(cat <<'EOF'
feat(providers): InternalProviderExecutor + BUNDLED_PROVIDERS

InternalProviderExecutor resolves an fsspec DirFileSystem handler
from filebrowser_path, runs Provider.read(), and POSTs the bytes to
/v1/joblib/providers/{id}/results using the internal worker JWT.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: `ensure_internal_worker_row` + host-app lifespan wiring

**Files:**
- Modify: `src/zndraw/database.py`
- Test: `tests/zndraw/test_lifespan.py` (append) — or a dedicated `test_internal_providers_lifespan.py` if the file is already large

> Check `tests/zndraw/test_lifespan.py` first; if it's >500 lines or structurally messy, create `tests/zndraw/test_internal_providers_lifespan.py` instead.

- [ ] **Step 8.1: Add `ensure_internal_worker_row` helper.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw/database.py`, beside `ensure_internal_worker` (~line 98), add:

```python
async def ensure_internal_worker_row(
    session: AsyncSession,
    user_id: uuid.UUID,
) -> uuid.UUID:
    """Create or reuse a Worker row owned by the internal worker user.

    Idempotent. Returns the worker id.
    """
    from sqlmodel import select as _select

    from zndraw_joblib.models import Worker

    result = await session.exec(
        _select(Worker).where(Worker.user_id == user_id).limit(1)
    )
    existing = result.one_or_none()
    if existing is not None:
        return existing.id
    worker = Worker(user_id=user_id)
    session.add(worker)
    await session.commit()
    await session.refresh(worker)
    return worker.id
```

- [ ] **Step 8.2: Add `_collect_providers`.**

In the same file, near `_collect_extensions` (~line 64), add:

```python
def _collect_providers() -> list[type]:
    """Collect all built-in provider classes to register at @internal."""
    from zndraw.providers import BUNDLED_PROVIDERS

    return list(BUNDLED_PROVIDERS)
```

- [ ] **Step 8.3: Wire `init_database` to seed provider rows.**

In `init_database` (~line 143), after the `ensure_internal_worker(...)` call, add:

```python
    # Seed internal worker row + @internal provider rows
    settings_ = settings  # local alias for clarity
    if settings_.filebrowser_path.lower() != "none":
        from zndraw_joblib.registry import ensure_internal_providers

        async with session_maker() as session:
            # Look up the internal worker user to get its id
            result = await session.exec(
                select(User).where(User.email == settings_.internal_worker_email)
            )
            internal_user = result.one()
            worker_id = await ensure_internal_worker_row(session, internal_user.id)

        await ensure_internal_providers(
            _collect_providers(),
            session_maker,
            user_id=internal_user.id,
            worker_id=worker_id,
        )
```

- [ ] **Step 8.4: Wire the lifespan to register taskiq tasks + attach registry to `app.state`.**

In `lifespan` (~line 192), locate the block after `register_internal_jobs` / `register_internal_tasks` (~line 295-307) and, in BOTH branches (`if settings.init_db_on_startup:` and `else:`), after the extensions registration, add:

```python
        # Register @internal provider tasks
        if settings.filebrowser_path.lower() != "none":
            from zndraw.providers.executor import InternalProviderExecutor
            from zndraw_joblib.registry import register_internal_providers

            provider_executor = InternalProviderExecutor(
                base_url=f"http://{executor_host}:{settings.port}",
                filebrowser_path=str(Path(settings.filebrowser_path).resolve()),
            )
            provider_registry = register_internal_providers(
                broker, _collect_providers(), provider_executor
            )
            app.state.internal_provider_registry = provider_registry
        else:
            app.state.internal_provider_registry = None
```

Also add `from pathlib import Path` at the top of the file if not already present.

- [ ] **Step 8.5: Update the external-worker broker.**

In `/Users/fzills/tools/zndraw-fastapi/src/zndraw/broker.py`, after the existing `register_internal_tasks(...)` call, add:

```python
from zndraw.providers.executor import InternalProviderExecutor
from zndraw_joblib.registry import register_internal_providers

if settings.filebrowser_path.lower() != "none":
    from pathlib import Path

    from zndraw.database import _collect_providers

    provider_executor = InternalProviderExecutor(
        base_url=server_url,
        filebrowser_path=str(Path(settings.filebrowser_path).resolve()),
    )
    register_internal_providers(broker, _collect_providers(), provider_executor)
```

- [ ] **Step 8.6: Write an integration test.**

Append to `/Users/fzills/tools/zndraw-fastapi/tests/zndraw/test_providers_filesystem.py`:

```python
# =============================================================================
# @internal default filesystem provider
# =============================================================================


def test_default_internal_filesystem_listed(server):
    """The default @internal:filesystem:FilesystemRead is listed in every room."""
    vis = ZnDraw(url=server)
    try:
        providers = _list_providers(vis)
        internal = [
            p
            for p in providers
            if p.full_name == "@internal:filesystem:FilesystemRead"
        ]
        assert len(internal) == 1
    finally:
        vis.disconnect()


def test_default_internal_filesystem_read(server, tmp_path, monkeypatch):
    """A read against @internal:filesystem:FilesystemRead returns local files."""
    # Seed a file in the server's cwd — but that's not this process's tmp_path,
    # so use a scoped ZNDRAW_SERVER_FILEBROWSER_PATH fixture instead.
    # This particular test relies on the cwd-based default; assert the endpoint
    # responds (200 with a list) rather than asserting specific filenames.
    vis = ZnDraw(url=server)
    try:
        items = _poll_until_200(
            vis, "@internal:filesystem:FilesystemRead", {"path": "/"}
        )
        assert isinstance(items, list)
    finally:
        vis.disconnect()
```

- [ ] **Step 8.7: Add a `filebrowser_path="none"` integration test using `server_factory`.**

Append to the same file:

```python
def test_filebrowser_path_none_disables_default(server_factory):
    """Setting ZNDRAW_SERVER_FILEBROWSER_PATH=none drops the @internal provider."""
    instance = server_factory({"ZNDRAW_SERVER_FILEBROWSER_PATH": "none"})
    vis = ZnDraw(url=instance.url)
    try:
        providers = _list_providers(vis)
        assert not any(p.room_id == "@internal" for p in providers)
    finally:
        vis.disconnect()
```

- [ ] **Step 8.8: Run the integration tests.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw/test_providers_filesystem.py -v
```

Expected: all pass (existing tests + 3 new ones).

- [ ] **Step 8.9: Run the full joblib + provider suites.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest tests/zndraw_joblib tests/zndraw/test_providers_filesystem.py tests/zndraw/test_providers_executor.py tests/zndraw/worker/test_internal.py -q
```

Expected: all pass.

- [ ] **Step 8.10: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add src/zndraw/database.py src/zndraw/broker.py tests/zndraw/test_providers_filesystem.py && git commit -m "$(cat <<'EOF'
feat(zndraw): register default @internal filesystem provider

On startup, the server seeds a ProviderRecord at
@internal:filesystem:FilesystemRead and registers the matching
taskiq task against the broker. Reads dispatch to the taskiq worker
which resolves an fsspec DirFileSystem rooted at filebrowser_path.

Setting ZNDRAW_SERVER_FILEBROWSER_PATH=none disables the default.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 5: Frontend — hide activity-bar icon when empty

### Task 9: `useHasFilesystemProviders` hook

**Files:**
- Create: `frontend/src/hooks/useHasFilesystemProviders.ts`

- [ ] **Step 9.1: Create the hook.**

Write `/Users/fzills/tools/zndraw-fastapi/frontend/src/hooks/useHasFilesystemProviders.ts`:

```typescript
import { useQuery } from "@tanstack/react-query";
import { useEffect } from "react";
import { listProviders } from "../myapi/client";
import { queryClient } from "../myapi/queryClient";
import { socket } from "../socket";
import { useAppStore } from "../store";

/**
 * Returns `true` when the current room has at least one provider of
 * category `"filesystem"`. Invalidates on the `providers_invalidate`
 * Socket.IO event so a newly-registered provider lights the icon up
 * without a page reload.
 */
export function useHasFilesystemProviders(): boolean {
	const roomId = useAppStore((s) => s.roomId);

	const { data } = useQuery({
		queryKey: ["filesystemProviders", roomId],
		queryFn: () => listProviders(roomId!, "filesystem"),
		enabled: !!roomId,
		retry: false,
		staleTime: 5_000,
	});

	useEffect(() => {
		const onInvalidate = () => {
			queryClient.invalidateQueries({
				queryKey: ["filesystemProviders", roomId],
			});
		};
		socket.on("providers_invalidate", onInvalidate);
		return () => {
			socket.off("providers_invalidate", onInvalidate);
		};
	}, [roomId]);

	return (data?.length ?? 0) > 0;
}
```

> Before writing, open `frontend/src/myapi/queryClient.ts` and `frontend/src/socket.ts` to confirm the exports used exist. If `queryClient` is exported from a different path (e.g. `frontend/src/App.tsx` or `frontend/src/main.tsx`), adjust the import. If `socket` is exported with a different event-bus API (e.g. typed `tsio.on(Event, handler)`), mirror the existing pattern used elsewhere in the codebase — e.g. look at how `FilesystemPanel.tsx` or `RoomsPanel.tsx` subscribe to invalidations for the right idiom.

- [ ] **Step 9.2: Typecheck.**

```bash
cd /Users/fzills/tools/zndraw-fastapi/frontend && bun --bun tsc --noEmit
```

Expected: zero errors. If `queryClient` is imported differently in the codebase, adapt.

- [ ] **Step 9.3: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add frontend/src/hooks/useHasFilesystemProviders.ts && git commit -m "$(cat <<'EOF'
feat(frontend): useHasFilesystemProviders hook

TanStack-Query hook that reports whether the current room has any
provider of category "filesystem". Invalidates on
providers_invalidate socket events.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 10: Hide `filesystem` icon in `ActivityBar`

**Files:**
- Modify: `frontend/src/panels/ActivityBar.tsx`

- [ ] **Step 10.1: Apply the filter.**

In `/Users/fzills/tools/zndraw-fastapi/frontend/src/panels/ActivityBar.tsx`:

1. Add the import near the top:

```typescript
import { useHasFilesystemProviders } from "../hooks/useHasFilesystemProviders";
```

2. Inside the `ActivityBar` function, after the existing `useAppStore` / `useCallback` hooks, add:

```typescript
	const hasFilesystemProviders = useHasFilesystemProviders();
	const filteredIcons = hasFilesystemProviders
		? icons
		: icons.filter((id) => id !== "filesystem");
```

3. Replace the two references to `icons` that drive rendering. The state-machine line (~144) stays using `icons.length` for drag-target sizing — actually it must use `filteredIcons.length` so a bar with only the filesystem icon reverts to its "empty" sliver state when the provider vanishes. Update:

```typescript
	const state: SliverState =
		filteredIcons.length === 0
			? isDragActive
				? isOverZone
					? "over-zone"
					: "hot"
				: "sliver"
			: "full";
```

And the render loop (`icons.map((id, idx) => {...})`) becomes:

```typescript
			{filteredIcons.map((id, idx) => {
```

- [ ] **Step 10.2: Typecheck + lint.**

```bash
cd /Users/fzills/tools/zndraw-fastapi/frontend && bun --bun tsc --noEmit && bun run lint
```

Expected: zero errors.

- [ ] **Step 10.3: Manual visual check.**

Start both servers:

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run zndraw --no-browser &
cd /Users/fzills/tools/zndraw-fastapi/frontend && bun run dev &
```

Using playwright-cli (see `skills/playwright-cli/SKILL.md` via `/playwright-cli`):

1. Navigate to `http://localhost:5173/`.
2. Take a screenshot of the left activity bar — the Folder icon should be visible (because the default filesystem provider is registered).
3. Stop the server, restart it with `ZNDRAW_SERVER_FILEBROWSER_PATH=none uv run zndraw --no-browser`.
4. Reload the page.
5. Take a screenshot — the Folder icon should be gone.

Check the screenshots! Compare "visible vs hidden" — do they match the expectation from step 10.1?

- [ ] **Step 10.4: Commit.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && git add frontend/src/panels/ActivityBar.tsx && git commit -m "$(cat <<'EOF'
feat(frontend): hide filesystem icon when no providers exist

ActivityBar filters the filesystem panel id out of the rendered
icons when the active room has no filesystem providers. Matches the
ZNDRAW_SERVER_FILEBROWSER_PATH=none behavior.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 6: End-to-end verification

### Task 11: Manual end-to-end sanity + final test sweep

**Files:** none — verification only.

- [ ] **Step 11.1: Start the real server in the repo root (so cwd = repo).**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run zndraw --no-browser &
cd /Users/fzills/tools/zndraw-fastapi/frontend && bun run dev &
```

- [ ] **Step 11.2: Open the UI and exercise the file browser.**

With playwright-cli:

1. Navigate to `http://localhost:5173/`.
2. Open the Files panel (click the folder icon in the left activity bar).
3. Confirm a listing appears rooted at the repo root — expect to see entries like `src/`, `tests/`, `frontend/`, `pyproject.toml`.
4. Click into `src/` — the panel should navigate one level deeper.
5. Screenshot the panel. Look at it: is the listing coherent? Any console errors?

If the panel shows "No filesystems registered for this room", check that the `@internal` provider is in fact listed via:

```bash
curl -s http://localhost:8000/v1/joblib/rooms/@global/providers | jq '.items[] | select(.room_id == "@internal")'
```

- [ ] **Step 11.3: Run the full backend suite.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uv run pytest -q
```

Expected: no new failures vs. the baseline from Task 0. Record anything that regresses.

- [ ] **Step 11.4: Typecheck frontend.**

```bash
cd /Users/fzills/tools/zndraw-fastapi/frontend && bun --bun tsc --noEmit && bun run lint
```

Expected: zero errors.

- [ ] **Step 11.5: Lint backend.**

```bash
cd /Users/fzills/tools/zndraw-fastapi && uvx prek run --all-files
```

Expected: zero errors.

- [ ] **Step 11.6: Final commit is unnecessary (everything is already committed).** If you touched anything during verification, either revert or commit it separately with a clear message.

---

## Spec-coverage self-review

| Spec requirement | Implemented in |
|---|---|
| `register_internal_providers(broker, providers, executor)` | Task 2 |
| Each `@internal` provider becomes `@internal:<category>:<ClassName>` taskiq task | Task 2 |
| Task handler resolves handler, runs `read()`, POSTs result | Task 7 |
| DB invariant: `ProviderRecord` rows per `@internal` provider, idempotent | Task 3 |
| Dispatch fork in `read_provider` | Task 6 |
| Visibility widening of `_room_provider_filter` to include `@internal` | Task 4 |
| `filebrowser_path` setting, default `"."`, `"none"` sentinel disables | Task 1 |
| `_collect_providers()` returning bundled classes (v1: `FilesystemRead`) | Task 8 |
| `broker.py` + `database.py` wire-up | Tasks 7–8 |
| `FilesystemRead` unchanged | (untouched) |
| Frontend hides filesystem icon when no filesystem providers | Tasks 9–10 |
| Re-queries on `providers_invalidate` | Task 9 |
| docker-compose unchanged | (no change needed) |
| Typecheck + lint zero errors | Task 0 + Task 11 |
| Unit: `register_internal_providers` creates ProviderRecord rows | Task 3 |
| Integration: list providers returns `@internal` by default | Task 8 |
| Integration: `FILEBROWSER_PATH=none` → empty list | Task 8 |
| Integration: reading via existing endpoint returns listings | Task 8 |
| Manual: `uv run zndraw`, open room, browse filesystem | Task 11 |
| `@internal` provider not user-deletable | Task 5 |
