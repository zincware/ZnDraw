# Worker Auth JWT Redesign — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the shared worker password with per-task JWT tokens minted at dispatch time, eliminating the well-known default superuser credential.

**Architecture:** The server auto-generates a random UUID password for the internal worker user on startup (write-only, never shared). When dispatching internal extension tasks to TaskIQ, the server mints a short-lived JWT and includes it in the job payload. The TaskIQ worker uses this token directly via `ZnDraw(token=...)` — no password or email needed.

**Tech Stack:** FastAPI, fastapi-users (JWTStrategy), TaskIQ, pydantic-settings

**Spec:** `docs/superpowers/specs/2026-04-01-worker-auth-jwt-design.md`

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `src/zndraw/config.py` | Modify | Add `internal_worker_email`, remove `worker_password` |
| `src/zndraw/database.py` | Modify | Auto-gen UUID password, remove `WORKER_EMAIL` constant, add `lookup_worker_user()` |
| `src/zndraw/executor.py` | Modify | Drop email/password fields, accept `token` param |
| `src/zndraw_joblib/registry.py` | Modify | Add `token` to protocol, task fn, and closure |
| `src/zndraw_joblib/dependencies.py` | Modify | Add `WorkerTokenDep` |
| `src/zndraw_joblib/router.py` | Modify | Inject `WorkerTokenDep`, pass to `kiq()` |
| `src/zndraw/broker.py` | Modify | Simplify executor — only `base_url` |
| `src/zndraw/routes/auth.py` | Modify | Add login guard for internal worker email |
| `docker/templates/.env` | Modify | Remove `ZNDRAW_SERVER_WORKER_PASSWORD` |
| `tests/zndraw/test_worker_auth.py` | Create | Login guard + worker token tests |

---

### Task 1: Config — add `internal_worker_email`, remove `worker_password`

**Files:**
- Modify: `src/zndraw/config.py:67-69`

- [ ] **Step 1: Remove `worker_password`, add `internal_worker_email`**

In `src/zndraw/config.py`, replace the `worker_password` field with `internal_worker_email`:

```python
    # Auth
    guest_password: SecretStr = SecretStr("zndraw")
    internal_worker_email: str = "worker@internal.user"
```

- [ ] **Step 2: Verify config loads**

Run: `uv run python -c "from zndraw.config import Settings; s = Settings(); print(s.internal_worker_email); assert not hasattr(s, 'worker_password')"`
Expected: prints `worker@internal.user`

- [ ] **Step 3: Commit**

```bash
git add src/zndraw/config.py
git commit -m "refactor: replace worker_password with internal_worker_email in Settings"
```

---

### Task 2: Database — auto-gen UUID password, add `lookup_worker_user`

**Files:**
- Modify: `src/zndraw/database.py:95-137` (WORKER_EMAIL constant + ensure_internal_worker)
- Modify: `src/zndraw/database.py:178` (init_database call)
- Modify: `src/zndraw/database.py:288-292` (lifespan executor creation)

- [ ] **Step 1: Update `ensure_internal_worker` — remove password param, auto-gen UUID**

In `src/zndraw/database.py`, remove the `WORKER_EMAIL` module constant. Update `ensure_internal_worker` to take `email` (from settings) and generate a random UUID password internally:

```python
import uuid

async def ensure_internal_worker(
    session: AsyncSession,
    email: str,
) -> None:
    """Create or update the internal worker superuser.

    Idempotent — safe to call on every startup. The password is a random
    UUID generated each time — it is never used for login (public login
    is blocked for this email).

    Parameters
    ----------
    session
        Async database session.
    email
        Internal worker email from ``Settings.internal_worker_email``.
    """
    password_helper = PasswordHelper()

    result = await session.execute(
        select(User).where(User.email == email)  # type: ignore[arg-type]
    )
    existing = result.scalar_one_or_none()

    hashed = password_helper.hash(str(uuid.uuid4()))

    if existing is None:
        worker = User(
            email=email,
            hashed_password=hashed,
            is_active=True,
            is_superuser=True,
            is_verified=True,
        )
        session.add(worker)
        await session.commit()
        log.info("Created internal worker user: %s", email)
    else:
        existing.hashed_password = hashed
        existing.is_superuser = True
        await session.commit()
        log.debug("Updated internal worker user: %s", email)
```

- [ ] **Step 2: Add `lookup_worker_user` helper**

Add this function after `ensure_internal_worker` in `database.py`:

```python
async def lookup_worker_user(session: AsyncSession, email: str) -> User:
    """Look up the internal worker user by email.

    Parameters
    ----------
    session
        Async database session.
    email
        Internal worker email from ``Settings.internal_worker_email``.

    Raises
    ------
    RuntimeError
        If the worker user does not exist (db-init not run).
    """
    result = await session.execute(
        select(User).where(User.email == email)  # type: ignore[arg-type]
    )
    user = result.scalar_one_or_none()
    if user is None:
        raise RuntimeError(
            f"Internal worker user '{email}' not found. "
            "Has the database been initialized (zndraw-db / init_db_on_startup)?"
        )
    return user
```

- [ ] **Step 3: Update `init_database` call**

In `init_database()` (around line 178), change:

```python
        await ensure_internal_worker(session, settings.worker_password)
```

to:

```python
        await ensure_internal_worker(session, settings.internal_worker_email)
```

- [ ] **Step 4: Update lifespan executor creation**

In `lifespan()` (around lines 288-292), simplify the executor — remove `worker_email` and `worker_password`:

```python
        executor = InternalExtensionExecutor(
            base_url=f"http://{executor_host}:{settings.port}",
        )
```

Also remove the `WORKER_EMAIL` import from the top of the file (line 95 constant is gone) and clean up the now-unused `from pydantic import SecretStr` import if `SecretStr` is no longer used in this file.

- [ ] **Step 5: Verify database init still works**

Run: `uv run python -c "import asyncio; from zndraw.database import init_database; asyncio.run(init_database())"`
Expected: completes without error

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/database.py
git commit -m "refactor: auto-gen worker password, add lookup_worker_user helper"
```

---

### Task 3: Executor + Protocol — accept `token` instead of credentials

**Files:**
- Modify: `src/zndraw/executor.py` (full file)
- Modify: `src/zndraw_joblib/registry.py:22-35` (InternalExecutor protocol)
- Modify: `src/zndraw_joblib/registry.py:63-72` (_make_task_fn closure)

- [ ] **Step 1: Simplify `InternalExtensionExecutor`**

Replace the entire `src/zndraw/executor.py` with:

```python
"""InternalExecutor implementation for zndraw-joblib.

Runs built-in extensions (modifiers, selections, analysis) in a thread
because the ZnDraw client uses synchronous HTTP and Socket.IO.

Status transitions use the same JobManager HTTP API as external workers.
"""

from __future__ import annotations

import asyncio
import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    pass

log = logging.getLogger(__name__)


@dataclass
class InternalExtensionExecutor:
    """Executes built-in extensions using a per-task JWT token.

    Only the server base URL is captured at creation time.
    The JWT token is passed per-task at call time.

    Parameters
    ----------
    base_url
        ZnDraw server URL.
    """

    base_url: str

    async def __call__(
        self,
        extension_cls: type,
        payload: dict[str, Any],
        room_id: str,
        task_id: str,
        token: str,
    ) -> None:
        """Execute extension via asyncio.to_thread."""
        base_url = self.base_url

        def _run() -> None:
            from zndraw.client import ZnDraw
            from zndraw_joblib.client import ClaimedTask

            vis = ZnDraw(
                url=base_url,
                room=room_id,
                token=token,
            )
            task: ClaimedTask | None = None
            try:
                instance = extension_cls(**payload)
                task = ClaimedTask(
                    task_id=task_id,
                    job_name="",
                    room_id=room_id,
                    extension=instance,
                )
                vis.jobs.start(task)
                instance.run(vis)
            except Exception as e:
                try:
                    if task is not None:
                        vis.jobs.fail(task, str(e))
                    else:
                        vis.jobs.fail_by_id(task_id, str(e))
                except Exception:
                    log.exception("Could not report failure for task %s", task_id)
                raise
            else:
                vis.jobs.complete(task)
            finally:
                vis.disconnect()

        await asyncio.to_thread(_run)
```

- [ ] **Step 2: Update `InternalExecutor` protocol**

In `src/zndraw_joblib/registry.py`, update the protocol (lines 22-35):

```python
class InternalExecutor(Protocol):
    """Protocol for the host-provided executor callback.

    The server base URL is captured at creation time.
    A per-task JWT token is passed at call time.
    """

    async def __call__(
        self,
        extension_cls: type[Extension],
        payload: dict[str, Any],
        room_id: str,
        task_id: str,
        token: str,
    ) -> None: ...
```

- [ ] **Step 3: Update `_make_task_fn` closure to pass `token`**

In `src/zndraw_joblib/registry.py`, update the closure (lines 63-72):

```python
        def _make_task_fn(
            cls: type[Extension] = ext_cls, ex: InternalExecutor = executor
        ):
            async def _execute(
                task_id: str, room_id: str, payload: dict[str, Any], token: str
            ) -> None:
                await ex(cls, payload, room_id, task_id, token)

            return _execute
```

- [ ] **Step 4: Verify import works**

Run: `uv run python -c "from zndraw.executor import InternalExtensionExecutor; e = InternalExtensionExecutor(base_url='http://localhost:8000'); print('OK')"`
Expected: prints `OK`

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/executor.py src/zndraw_joblib/registry.py
git commit -m "refactor: executor accepts per-task JWT token instead of static credentials"
```

---

### Task 4: WorkerTokenDep + dispatch wiring

**Files:**
- Modify: `src/zndraw_joblib/dependencies.py` (add `WorkerTokenDep`)
- Modify: `src/zndraw_joblib/router.py:583-626` (submit_task)

- [ ] **Step 1: Add `WorkerTokenDep` to dependencies.py**

Add the following at the end of `src/zndraw_joblib/dependencies.py`:

```python
from typing import Protocol as TypingProtocol


class WorkerTokenFactory(TypingProtocol):
    """Protocol for minting a JWT for the internal worker user."""

    async def __call__(self) -> str: ...


async def _no_worker_token() -> str:
    raise NotImplementedError(
        "WorkerTokenFactory not configured — host app must override get_worker_token"
    )


async def get_worker_token() -> str:
    """Return a fresh JWT for the internal worker user.

    Host apps must override this dependency via
    ``app.dependency_overrides[get_worker_token]``.
    """
    return await _no_worker_token()


WorkerTokenDep = Annotated[str, Depends(get_worker_token)]
```

- [ ] **Step 2: Wire `WorkerTokenDep` into `submit_task`**

In `src/zndraw_joblib/router.py`, add the import at the top (with the other dependency imports around line 24):

```python
from zndraw_joblib.dependencies import (
    FrameRoomCleanupDep,
    JobLibSettingsDep,
    ResultBackendDep,
    WorkerTokenDep,
    WritableRoomDep,
    get_internal_registry,
    get_tsio,
    request_hash,
    validate_room_id,
)
```

Then update the `submit_task` function signature (around line 583) to accept the new dependency:

```python
async def submit_task(
    room_id: WritableRoomDep,
    job_name: str,
    request: TaskSubmitRequest,
    response: Response,
    session: SessionDep,
    user: CurrentUserDep,
    internal_registry: InternalRegistryDep,
    tsio: TsioDep,
    worker_token: WorkerTokenDep,
):
```

Then update the `kiq()` call (around line 622) to pass `token`:

```python
            await internal_registry.tasks[job.full_name].kiq(
                task_id=str(task.id),
                room_id=room_id,
                payload=request.payload,
                token=worker_token,
            )
```

- [ ] **Step 3: Wire the dependency override in the lifespan**

In `src/zndraw/database.py`, inside the `lifespan()` function, after the broker setup (after line ~307), add the dependency override. First add the necessary imports at the top of `database.py`:

```python
from zndraw_joblib.dependencies import get_worker_token
```

Then in the lifespan, after `await broker.startup()` (around line 307):

```python
        # Wire WorkerTokenDep — mints JWTs for the internal worker user
        from fastapi_users.authentication import JWTStrategy

        async def _mint_worker_token() -> str:
            async with app.state.session_maker() as session:
                worker = await lookup_worker_user(
                    session, settings.internal_worker_email
                )
            strategy = JWTStrategy(
                secret=auth_settings.secret_key.get_secret_value(),
                lifetime_seconds=auth_settings.token_lifetime_seconds,
            )
            return await strategy.write_token(worker)

        app.dependency_overrides[get_worker_token] = _mint_worker_token
```

- [ ] **Step 4: Verify imports**

Run: `uv run python -c "from zndraw_joblib.dependencies import WorkerTokenDep, get_worker_token; print('OK')"`
Expected: prints `OK`

- [ ] **Step 5: Commit**

```bash
git add src/zndraw_joblib/dependencies.py src/zndraw_joblib/router.py src/zndraw/database.py
git commit -m "feat: mint per-task JWT for internal worker via WorkerTokenDep"
```

---

### Task 5: Simplify `broker.py`

**Files:**
- Modify: `src/zndraw/broker.py` (full file)

- [ ] **Step 1: Remove credential references from broker**

Replace `src/zndraw/broker.py`:

```python
"""Module-level broker for external TaskIQ workers.

Usage::

    taskiq worker zndraw.broker:broker

The broker connects to Redis (via ``ZNDRAW_SERVER_REDIS_URL``) and registers
all built-in extensions. The executor connects back to the FastAPI
server at ``ZNDRAW_SERVER_INTERNAL_URL``.
"""

from taskiq_redis import ListQueueBroker

from zndraw.config import Settings
from zndraw.database import _collect_extensions
from zndraw.executor import InternalExtensionExecutor
from zndraw_joblib import register_internal_tasks

settings = Settings()

if settings.redis_url is None:
    raise RuntimeError(
        "ZNDRAW_SERVER_REDIS_URL must be set for external TaskIQ workers. "
        "External workers cannot use fakeredis."
    )

broker = ListQueueBroker(settings.redis_url)

server_url = settings.internal_url or f"http://localhost:{settings.port}"

executor = InternalExtensionExecutor(
    base_url=server_url,
)

register_internal_tasks(broker, _collect_extensions(), executor)
```

- [ ] **Step 2: Verify broker module loads**

Run: `ZNDRAW_SERVER_REDIS_URL=redis://localhost:6379 uv run python -c "from zndraw.broker import broker; print('OK')"`
Expected: prints `OK` (may warn about Redis connection but module loads)

- [ ] **Step 3: Commit**

```bash
git add src/zndraw/broker.py
git commit -m "refactor: simplify broker — executor only needs base_url"
```

---

### Task 6: Block worker login

**Files:**
- Modify: `src/zndraw/routes/auth.py:46-50`
- Create: `tests/zndraw/test_worker_auth.py`

- [ ] **Step 1: Write failing test for login guard**

Create `tests/zndraw/test_worker_auth.py`:

```python
"""Tests for internal worker auth security."""

import pytest
from httpx import AsyncClient


@pytest.mark.anyio
async def test_worker_login_blocked(client: AsyncClient, settings) -> None:
    """POST /v1/auth/jwt/login with the worker email must return 403."""
    resp = await client.post(
        "/v1/auth/jwt/login",
        data={
            "username": settings.internal_worker_email,
            "password": "does-not-matter",
        },
    )
    assert resp.status_code == 403


@pytest.mark.anyio
async def test_regular_login_still_works(client: AsyncClient) -> None:
    """A non-worker email should not be blocked by the guard."""
    # Register a user first
    resp = await client.post(
        "/v1/auth/register",
        json={
            "email": "test-login@example.com",
            "password": "testpassword123",
        },
    )
    assert resp.status_code == 201

    # Login should work
    resp = await client.post(
        "/v1/auth/jwt/login",
        data={
            "username": "test-login@example.com",
            "password": "testpassword123",
        },
    )
    assert resp.status_code == 200
    assert "access_token" in resp.json()
```

- [ ] **Step 2: Run tests — verify they fail**

Run: `uv run pytest tests/zndraw/test_worker_auth.py -v`
Expected: `test_worker_login_blocked` FAILS (currently returns 400, not 403)

- [ ] **Step 3: Add login guard route**

In `src/zndraw/routes/auth.py`, add a custom login endpoint BEFORE the `router.include_router(fastapi_users.get_auth_router(...))` call. The custom route shadows the fastapi-users login route because FastAPI uses first-match routing.

Add these imports at the top:

```python
from fastapi import APIRouter, Depends, HTTPException
from fastapi.security import OAuth2PasswordRequestForm
from fastapi_users.authentication import JWTStrategy

from zndraw.config import Settings, SettingsDep, get_zndraw_settings
```

Then add the custom login endpoint before the `router.include_router` for JWT (before line 47):

```python
@router.post("/jwt/login")
async def login(
    settings: SettingsDep,
    auth_settings: AuthSettingsDep,
    user_manager: Annotated[UserManager, Depends(get_user_manager)],
    credentials: OAuth2PasswordRequestForm = Depends(),
) -> dict:
    """Login endpoint with internal worker email guard.

    Shadows the fastapi-users login route (first-match wins in FastAPI)
    to block the internal service account from public login.
    """
    if credentials.username == settings.internal_worker_email:
        raise HTTPException(
            status_code=403,
            detail="Internal service accounts cannot log in via the public API",
        )

    user = await user_manager.authenticate(credentials)
    if user is None or not user.is_active:
        raise HTTPException(status_code=400, detail="LOGIN_BAD_CREDENTIALS")

    strategy = JWTStrategy(
        secret=auth_settings.secret_key.get_secret_value(),
        lifetime_seconds=auth_settings.token_lifetime_seconds,
    )
    token = await strategy.write_token(user)
    return {"access_token": token, "token_type": "bearer"}
```

Note: `SettingsDep` is already defined in `config.py` — import it. The full import line becomes:

```python
from zndraw.config import Settings, SettingsDep, get_zndraw_settings
```

Remove the now-redundant standalone `Settings` and `get_zndraw_settings` imports if they were separate.

- [ ] **Step 4: Run tests — verify they pass**

Run: `uv run pytest tests/zndraw/test_worker_auth.py -v`
Expected: both tests PASS

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/routes/auth.py tests/zndraw/test_worker_auth.py
git commit -m "feat: block internal worker email from public login endpoint"
```

---

### Task 7: Docker + docs cleanup and grep verification

**Files:**
- Modify: `docker/templates/.env:10`
- Modify: `docs/superpowers/plans/2026-03-25-pydantic-settings-phase1.md:332` (if reference exists)

- [ ] **Step 1: Remove `ZNDRAW_SERVER_WORKER_PASSWORD` from Docker `.env` template**

In `docker/templates/.env`, remove this line:

```
ZNDRAW_SERVER_WORKER_PASSWORD=zndraw-worker
```

- [ ] **Step 2: Clean up docs references**

In `docs/superpowers/plans/2026-03-25-pydantic-settings-phase1.md`, find line 332:

```
- `ZNDRAW_WORKER_PASSWORD` → `ZNDRAW_SERVER_WORKER_PASSWORD`
```

Remove that line (the setting no longer exists).

- [ ] **Step 3: Run grep verification**

Run: `grep -r 'worker_password\|WORKER_PASSWORD' --include='*.py' --include='*.md' --include='*.yaml' --include='*.yml' --include='*.env' . | grep -v '.git/' | grep -v 'worker-auth-jwt-design.md' | grep -v 'worker-auth-jwt.md'`

Expected: **zero hits**. If any remain, fix them before proceeding.

- [ ] **Step 4: Commit**

```bash
git add docker/templates/.env docs/superpowers/plans/2026-03-25-pydantic-settings-phase1.md
git commit -m "chore: remove all worker_password references from Docker and docs"
```

---

### Task 8: Integration test — full dispatch with JWT

**Files:**
- Modify: `tests/zndraw/test_worker_auth.py`

- [ ] **Step 1: Write integration test for WorkerTokenDep**

Add to `tests/zndraw/test_worker_auth.py`:

```python
@pytest.mark.anyio
async def test_worker_token_dep_mints_valid_jwt(
    client: AsyncClient, settings
) -> None:
    """The WorkerTokenDep should mint a JWT for the internal worker user."""
    from zndraw_joblib.dependencies import get_worker_token

    # The dependency is overridden in lifespan — call it through the app
    app = client._transport.app  # type: ignore[union-attr]
    override = app.dependency_overrides.get(get_worker_token)
    assert override is not None, "get_worker_token override not configured"

    token = await override()
    assert isinstance(token, str)
    assert len(token) > 0

    # Verify token is valid by calling /v1/auth/users/me
    resp = await client.get(
        "/v1/auth/users/me",
        headers={"Authorization": f"Bearer {token}"},
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["email"] == settings.internal_worker_email
    assert data["is_superuser"] is True
```

- [ ] **Step 2: Run test**

Run: `uv run pytest tests/zndraw/test_worker_auth.py -v`
Expected: all tests PASS

- [ ] **Step 3: Run full test suite**

Run: `uv run pytest --timeout=60 -x -q`
Expected: all tests pass (no regressions)

- [ ] **Step 4: Commit**

```bash
git add tests/zndraw/test_worker_auth.py
git commit -m "test: add integration test for WorkerTokenDep JWT minting"
```
