# Progress Bar Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add progress tracking so long-running Python client operations (frame uploads, extensions) show real-time progress bars to all room clients.

**Architecture:** REST endpoints (POST/PATCH/DELETE) for progress mutations, Redis hash for ephemeral storage, Socket.IO broadcasts for real-time updates, tqdm subclass on the client for ergonomic usage.

**Tech Stack:** FastAPI, Redis (hash), Socket.IO (broadcast), tqdm (client-side), Pydantic, httpx

---

### Task 1: Add Redis Key Pattern

**Files:**
- Modify: `src/zndraw/redis.py:107-119` (add after Provider Dispatch Keys section)

**Step 1: Add the key method**

Add a new section at the bottom of `RedisKey`:

```python
# =========================================================================
# Progress Tracker Keys
# =========================================================================

@staticmethod
def room_progress(room_id: str) -> str:
    """Hash of active progress trackers in a room.

    Field: progress_id, Value: JSON ProgressResponse.
    """
    return f"room:{room_id}:progress"
```

**Step 2: Commit**

```bash
git add src/zndraw/redis.py
git commit -m "feat: add RedisKey.room_progress for progress tracking"
```

---

### Task 2: Add Pydantic Models (Schemas + Socket Events)

**Files:**
- Modify: `src/zndraw/schemas.py` (add after Screenshot schemas, ~line 507)
- Modify: `src/zndraw/socket_events.py` (add after Typing broadcast, ~line 227)

**Step 1: Add REST schemas to `schemas.py`**

Add at the end of the file:

```python
# =============================================================================
# Progress Schemas
# =============================================================================


class ProgressCreate(BaseModel):
    """Request to start a new progress tracker."""

    progress_id: str
    description: str


class ProgressPatch(BaseModel):
    """Request to update an existing progress tracker."""

    description: str | None = None
    progress: float | None = None


class ProgressResponse(BaseModel):
    """Response for a single progress tracker."""

    progress_id: str
    description: str
    progress: float | None = None
```

**Step 2: Add socket broadcast events to `socket_events.py`**

Add at the bottom of the Broadcast Models section:

```python
class ProgressStart(BaseModel):
    """Broadcast when a new progress tracker is created."""

    progress_id: str
    description: str


class ProgressUpdate(BaseModel):
    """Broadcast when a progress tracker is updated."""

    progress_id: str
    description: str | None = None
    progress: float | None = None


class ProgressComplete(BaseModel):
    """Broadcast when a progress tracker finishes."""

    progress_id: str
```

**Step 3: Add `progress_trackers` to `RoomJoinResponse`**

In `socket_events.py`, modify `RoomJoinResponse` (line 61-69) to add:

```python
from zndraw.schemas import ProgressResponse

class RoomJoinResponse(BaseModel):
    """Response for room join."""

    room_id: str
    session_id: str
    step: int
    frame_count: int
    locked: bool
    camera_key: str | None = None
    progress_trackers: dict[str, ProgressResponse] = {}
```

Note: This adds an import of `ProgressResponse` from `schemas.py`. Check for circular imports — `socket_events.py` already imports from `schemas.py` (`RoomResponse`), so adding `ProgressResponse` is fine.

**Step 4: Commit**

```bash
git add src/zndraw/schemas.py src/zndraw/socket_events.py
git commit -m "feat: add progress Pydantic models and socket events"
```

---

### Task 3: Add Problem Type for Missing Progress Tracker

**Files:**
- Modify: `src/zndraw/exceptions.py`

**Step 1: Add `ProgressNotFound` problem type**

Add after `MessageNotFound` (line 353):

```python
class ProgressNotFound(ProblemType):
    """The requested progress tracker does not exist.

    This error occurs when attempting to update or complete a
    progress tracker that does not exist or has already completed.
    """

    title: ClassVar[str] = "Not Found"
    status: ClassVar[int] = 404

    @classmethod
    def raise_for_client(cls, problem: "ProblemDetail") -> NoReturn:
        raise KeyError(problem.detail or problem.title)
```

**Step 2: Register it in `PROBLEM_TYPES`**

Add `ProgressNotFound` to the list in `PROBLEM_TYPES` (around line 491-522).

**Step 3: Commit**

```bash
git add src/zndraw/exceptions.py
git commit -m "feat: add ProgressNotFound problem type"
```

---

### Task 4: Write Progress Route Tests

**Files:**
- Create: `tests/test_progress.py`

**Step 1: Write the test file**

```python
"""Tests for Progress REST API endpoints."""

import json
from collections.abc import AsyncIterator
from unittest.mock import AsyncMock

import pytest
import pytest_asyncio
from conftest import MockSioServer, auth_header, create_test_room, create_test_user_in_db
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel

from zndraw.config import Settings
from zndraw.redis import RedisKey


@pytest_asyncio.fixture(name="progress_session")
async def progress_session_fixture() -> AsyncIterator[AsyncSession]:
    """Create a fresh in-memory async database session."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)
        factory = async_sessionmaker(
            bind=engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            yield session
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="mock_sio")
async def mock_sio_fixture() -> MockSioServer:
    return MockSioServer()


@pytest_asyncio.fixture(name="mock_redis")
async def mock_redis_fixture() -> AsyncMock:
    """Mock Redis with hash operation support for progress tracking."""
    redis = AsyncMock()
    redis.get = AsyncMock(return_value=None)  # no edit lock
    redis._progress_store: dict[str, dict[str, str]] = {}

    async def hset(key: str, field: str, value: str) -> int:
        if key not in redis._progress_store:
            redis._progress_store[key] = {}
        redis._progress_store[key][field] = value
        return 1

    async def hget(key: str, field: str) -> str | None:
        return redis._progress_store.get(key, {}).get(field)

    async def hdel(key: str, *fields: str) -> int:
        deleted = 0
        if key in redis._progress_store:
            for f in fields:
                if f in redis._progress_store[key]:
                    del redis._progress_store[key][f]
                    deleted += 1
        return deleted

    async def hgetall(key: str) -> dict[str, str]:
        return redis._progress_store.get(key, {})

    redis.hset = AsyncMock(side_effect=hset)
    redis.hget = AsyncMock(side_effect=hget)
    redis.hdel = AsyncMock(side_effect=hdel)
    redis.hgetall = AsyncMock(side_effect=hgetall)

    return redis


@pytest_asyncio.fixture(name="progress_client")
async def progress_client_fixture(
    progress_session: AsyncSession, mock_sio: MockSioServer, mock_redis: AsyncMock
) -> AsyncIterator[AsyncClient]:
    from zndraw_auth import get_session
    from zndraw_auth.settings import AuthSettings

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_tsio

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield progress_session

    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_tsio] = lambda: mock_sio
    app.dependency_overrides[get_redis] = lambda: mock_redis
    app.state.settings = Settings()
    app.state.auth_settings = AuthSettings()

    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        yield client

    app.dependency_overrides.clear()


# =============================================================================
# POST (create progress tracker)
# =============================================================================


@pytest.mark.asyncio
async def test_create_progress(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> None:
    """POST creates a progress tracker and broadcasts progress_start."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "test-123", "description": "Loading..."},
        headers=auth_header(token),
    )
    assert response.status_code == 201
    data = response.json()
    assert data["progress_id"] == "test-123"
    assert data["description"] == "Loading..."
    assert data["progress"] is None

    # Verify socket broadcast
    start_events = [e for e in mock_sio.emitted if e["event"] == "progress_start"]
    assert len(start_events) == 1
    assert start_events[0]["data"]["progress_id"] == "test-123"

    # Verify Redis storage
    key = RedisKey.room_progress(room.id)
    stored = await mock_redis.hget(key, "test-123")
    assert stored is not None
    assert json.loads(stored)["description"] == "Loading..."


@pytest.mark.asyncio
async def test_create_progress_requires_auth(
    progress_client: AsyncClient, progress_session: AsyncSession
) -> None:
    """POST without auth returns 401."""
    user, _ = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "test-123", "description": "Loading..."},
    )
    assert response.status_code == 401


# =============================================================================
# PATCH (update progress tracker)
# =============================================================================


@pytest.mark.asyncio
async def test_update_progress(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> None:
    """PATCH updates progress value and broadcasts progress_update."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    # Create first
    await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "test-123", "description": "Loading..."},
        headers=auth_header(token),
    )

    # Update
    response = await progress_client.patch(
        f"/v1/rooms/{room.id}/progress/test-123",
        json={"progress": 50.0},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["progress"] == 50.0
    assert data["description"] == "Loading..."

    # Verify socket broadcast
    update_events = [e for e in mock_sio.emitted if e["event"] == "progress_update"]
    assert len(update_events) == 1
    assert update_events[0]["data"]["progress"] == 50.0


@pytest.mark.asyncio
async def test_update_progress_not_found(
    progress_client: AsyncClient, progress_session: AsyncSession
) -> None:
    """PATCH for non-existent tracker returns 404."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.patch(
        f"/v1/rooms/{room.id}/progress/nonexistent",
        json={"progress": 50.0},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "progress-not-found" in response.json()["type"]


@pytest.mark.asyncio
async def test_update_progress_description(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_redis: AsyncMock,
) -> None:
    """PATCH can update description."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "test-123", "description": "Step 1"},
        headers=auth_header(token),
    )

    response = await progress_client.patch(
        f"/v1/rooms/{room.id}/progress/test-123",
        json={"description": "Step 2", "progress": 75.0},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    data = response.json()
    assert data["description"] == "Step 2"
    assert data["progress"] == 75.0


# =============================================================================
# DELETE (complete progress tracker)
# =============================================================================


@pytest.mark.asyncio
async def test_delete_progress(
    progress_client: AsyncClient,
    progress_session: AsyncSession,
    mock_sio: MockSioServer,
    mock_redis: AsyncMock,
) -> None:
    """DELETE removes tracker and broadcasts progress_complete."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    # Create first
    await progress_client.post(
        f"/v1/rooms/{room.id}/progress",
        json={"progress_id": "test-123", "description": "Loading..."},
        headers=auth_header(token),
    )

    # Delete
    response = await progress_client.delete(
        f"/v1/rooms/{room.id}/progress/test-123",
        headers=auth_header(token),
    )
    assert response.status_code == 204

    # Verify socket broadcast
    complete_events = [
        e for e in mock_sio.emitted if e["event"] == "progress_complete"
    ]
    assert len(complete_events) == 1
    assert complete_events[0]["data"]["progress_id"] == "test-123"

    # Verify removed from Redis
    key = RedisKey.room_progress(room.id)
    stored = await mock_redis.hget(key, "test-123")
    assert stored is None


@pytest.mark.asyncio
async def test_delete_progress_not_found(
    progress_client: AsyncClient, progress_session: AsyncSession
) -> None:
    """DELETE for non-existent tracker returns 404."""
    user, token = await create_test_user_in_db(progress_session)
    room = await create_test_room(progress_session, user)

    response = await progress_client.delete(
        f"/v1/rooms/{room.id}/progress/nonexistent",
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "progress-not-found" in response.json()["type"]


# =============================================================================
# Room not found
# =============================================================================


@pytest.mark.asyncio
async def test_progress_room_not_found(
    progress_client: AsyncClient, progress_session: AsyncSession
) -> None:
    """Operations on non-existent room return 404."""
    _, token = await create_test_user_in_db(progress_session)

    response = await progress_client.post(
        "/v1/rooms/nonexistent/progress",
        json={"progress_id": "x", "description": "y"},
        headers=auth_header(token),
    )
    assert response.status_code == 404
    assert "room-not-found" in response.json()["type"]
```

**Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_progress.py -v`
Expected: FAIL (routes don't exist yet)

**Step 3: Commit**

```bash
git add tests/test_progress.py
git commit -m "test: add progress REST endpoint tests (red)"
```

---

### Task 5: Implement Progress REST Endpoints

**Files:**
- Create: `src/zndraw/routes/progress.py`
- Modify: `src/zndraw/app.py` (register router)

**Step 1: Create the route file**

```python
"""Progress tracking REST API endpoints.

Ephemeral progress trackers stored in Redis (hash), broadcast via Socket.IO.
"""

import json

from fastapi import APIRouter, Response, status

from zndraw.dependencies import (
    CurrentUserDep,
    RedisDep,
    SessionDep,
    SioDep,
    room_channel,
    verify_room,
)
from zndraw.exceptions import (
    NotAuthenticated,
    ProgressNotFound,
    RoomNotFound,
    problem_responses,
)
from zndraw.redis import RedisKey
from zndraw.schemas import ProgressCreate, ProgressPatch, ProgressResponse
from zndraw.socket_events import ProgressComplete, ProgressStart, ProgressUpdate

router = APIRouter(prefix="/v1/rooms/{room_id}/progress", tags=["progress"])


@router.post(
    "",
    status_code=status.HTTP_201_CREATED,
    responses=problem_responses(NotAuthenticated, RoomNotFound),
)
async def create_progress(
    session: SessionDep,
    sio: SioDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
    request: ProgressCreate,
) -> ProgressResponse:
    """Create a new progress tracker in the room."""
    await verify_room(session, room_id)

    tracker = ProgressResponse(
        progress_id=request.progress_id,
        description=request.description,
    )
    await redis.hset(  # type: ignore[misc]
        RedisKey.room_progress(room_id),
        request.progress_id,
        tracker.model_dump_json(),
    )

    await sio.emit(
        ProgressStart(
            progress_id=request.progress_id,
            description=request.description,
        ),
        room=room_channel(room_id),
    )

    return tracker


@router.patch(
    "/{progress_id}",
    responses=problem_responses(NotAuthenticated, RoomNotFound, ProgressNotFound),
)
async def update_progress(
    session: SessionDep,
    sio: SioDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
    progress_id: str,
    request: ProgressPatch,
) -> ProgressResponse:
    """Update an existing progress tracker."""
    await verify_room(session, room_id)

    raw = await redis.hget(RedisKey.room_progress(room_id), progress_id)  # type: ignore[misc]
    if raw is None:
        raise ProgressNotFound.exception(
            f"Progress tracker {progress_id} not found"
        )

    current = json.loads(raw)
    if request.description is not None:
        current["description"] = request.description
    if request.progress is not None:
        current["progress"] = request.progress

    await redis.hset(  # type: ignore[misc]
        RedisKey.room_progress(room_id),
        progress_id,
        json.dumps(current),
    )

    await sio.emit(
        ProgressUpdate(
            progress_id=progress_id,
            description=request.description,
            progress=request.progress,
        ),
        room=room_channel(room_id),
    )

    return ProgressResponse(**current)


@router.delete(
    "/{progress_id}",
    status_code=status.HTTP_204_NO_CONTENT,
    responses=problem_responses(NotAuthenticated, RoomNotFound, ProgressNotFound),
)
async def delete_progress(
    session: SessionDep,
    sio: SioDep,
    redis: RedisDep,
    _current_user: CurrentUserDep,
    room_id: str,
    progress_id: str,
) -> Response:
    """Complete and remove a progress tracker."""
    await verify_room(session, room_id)

    deleted = await redis.hdel(RedisKey.room_progress(room_id), progress_id)  # type: ignore[misc]
    if not deleted:
        raise ProgressNotFound.exception(
            f"Progress tracker {progress_id} not found"
        )

    await sio.emit(
        ProgressComplete(progress_id=progress_id),
        room=room_channel(room_id),
    )

    return Response(status_code=status.HTTP_204_NO_CONTENT)
```

**Step 2: Register the router in `app.py`**

Add to imports (around line 24-40):
```python
from zndraw.routes.progress import router as progress_router
```

Add to the router list (around line 87):
```python
app.include_router(progress_router)
```

**Step 3: Run tests to verify they pass**

Run: `uv run pytest tests/test_progress.py -v`
Expected: ALL PASS

**Step 4: Run full test suite**

Run: `uv run pytest tests/ -x`
Expected: ALL PASS (no regressions)

**Step 5: Commit**

```bash
git add src/zndraw/routes/progress.py src/zndraw/app.py
git commit -m "feat: implement progress REST endpoints"
```

---

### Task 6: Include Progress Trackers in RoomJoinResponse

**Files:**
- Modify: `src/zndraw/socketio.py` (room_join handler, ~line 318-327)

**Step 1: Add Redis read of active trackers to `room_join`**

In the `room_join` handler, before the final `return RoomJoinResponse(...)` (around line 320), add:

```python
# Load active progress trackers from Redis
from zndraw.schemas import ProgressResponse
progress_raw = await redis.hgetall(RedisKey.room_progress(data.room_id))  # type: ignore[misc]
progress_trackers = {
    pid: ProgressResponse(**json.loads(pdata))
    for pid, pdata in progress_raw.items()
}
```

Then add `progress_trackers=progress_trackers` to the `RoomJoinResponse(...)` constructor.

Note: The `ProgressResponse` import should go at the top of the file with other imports from `zndraw.schemas`. The `json` module is already imported.

**Step 2: Run full tests**

Run: `uv run pytest tests/ -x`
Expected: ALL PASS

**Step 3: Commit**

```bash
git add src/zndraw/socketio.py
git commit -m "feat: include progress trackers in room join response"
```

---

### Task 7: Add Progress Methods to Python Client APIManager

**Files:**
- Modify: `src/zndraw/client.py` (add after Screenshot Operations, ~line 928)

**Step 1: Add progress methods to `APIManager`**

After the Screenshot Operations section (line 928), add:

```python
    # -------------------------------------------------------------------------
    # Progress Operations
    # -------------------------------------------------------------------------

    def progress_start(self, progress_id: str, description: str) -> dict[str, Any]:
        """Create a new progress tracker in the room."""
        response = self.http.post(
            f"/v1/rooms/{self.room_id}/progress",
            json={"progress_id": progress_id, "description": description},
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def progress_update(
        self,
        progress_id: str,
        description: str | None = None,
        progress: float | None = None,
    ) -> dict[str, Any]:
        """Update an existing progress tracker."""
        payload: dict[str, Any] = {}
        if description is not None:
            payload["description"] = description
        if progress is not None:
            payload["progress"] = progress
        response = self.http.patch(
            f"/v1/rooms/{self.room_id}/progress/{progress_id}",
            json=payload,
            headers=self._headers(),
        )
        self.raise_for_status(response)
        return response.json()

    def progress_complete(self, progress_id: str) -> None:
        """Complete and remove a progress tracker."""
        response = self.http.delete(
            f"/v1/rooms/{self.room_id}/progress/{progress_id}",
            headers=self._headers(),
        )
        self.raise_for_status(response)
```

**Step 2: Commit**

```bash
git add src/zndraw/client.py
git commit -m "feat: add progress methods to APIManager"
```

---

### Task 8: Implement ZnDrawTqdm

**Files:**
- Create: `src/zndraw/tqdm.py`
- Modify: `src/zndraw/__init__.py` (export)
- Modify: `pyproject.toml` (add tqdm dependency)

**Step 1: Add tqdm dependency**

Run: `uv add tqdm`

**Step 2: Create `src/zndraw/tqdm.py`**

```python
"""tqdm subclass that reports progress to a ZnDraw room."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any
from uuid import uuid4

from tqdm.auto import tqdm as std_tqdm

if TYPE_CHECKING:
    from zndraw.client import ZnDraw


class ZnDrawTqdm(std_tqdm):  # type: ignore[type-arg]
    """tqdm subclass that sends progress updates to a ZnDraw room.

    Uses tqdm's built-in ``mininterval`` throttling to limit REST calls.
    Terminal output is suppressed (headless).

    Parameters
    ----------
    *args : Any
        Positional arguments forwarded to tqdm.
    vis : ZnDraw
        The ZnDraw client instance (must have ``vis.api``).
    description : str
        Human-readable label shown in the UI.
    mininterval : float
        Minimum seconds between REST updates (default 0.5).
    **kwargs : Any
        Keyword arguments forwarded to tqdm.

    Examples
    --------
    >>> for atoms in ZnDrawTqdm(frames, vis=vis, description="Loading"):
    ...     vis.append(atoms)
    """

    def __init__(
        self,
        *args: Any,
        vis: ZnDraw,
        description: str = "Processing...",
        mininterval: float = 0.5,
        **kwargs: Any,
    ) -> None:
        kwargs["mininterval"] = mininterval
        super().__init__(*args, **kwargs)
        self._vis_api = vis.api
        self._progress_id = str(uuid4())
        self._vis_api.progress_start(self._progress_id, description)

    def display(self, *args: Any, **kwargs: Any) -> None:
        """Suppress terminal output (headless mode)."""

    def update(self, n: int = 1) -> bool | None:
        """Advance the progress bar and send a throttled REST update."""
        displayed = super().update(n)
        if displayed:
            d = self.format_dict
            total = d.get("total")
            progress = (d["n"] / total * 100) if total else None
            self._vis_api.progress_update(self._progress_id, progress=progress)
        return displayed

    def close(self) -> None:
        """Complete the progress tracker and clean up."""
        if not self.disable:
            self._vis_api.progress_complete(self._progress_id)
        super().close()
```

**Step 3: Export from `__init__.py`**

Modify `src/zndraw/__init__.py`:

```python
"""ZnDraw - Interactive visualization for atomistic simulations."""

from zndraw.client import ZnDraw
from zndraw.tqdm import ZnDrawTqdm

try:
    from zndraw._version import __version__, __version_tuple__
except ImportError:
    __version__ = "0.0.0-dev"
    __version_tuple__ = (0, 0, 0, "dev")

__all__ = ["ZnDraw", "ZnDrawTqdm", "__version__", "__version_tuple__"]
```

**Step 4: Run type checker**

Run: `uv run pyright .`
Expected: No new errors

**Step 5: Commit**

```bash
git add src/zndraw/tqdm.py src/zndraw/__init__.py pyproject.toml uv.lock
git commit -m "feat: add ZnDrawTqdm progress bar subclass"
```

---

### Task 9: Adapt Frontend Socket Handlers

**Files:**
- Modify: `frontend/src/hooks/useSocketManager.ts` (~lines 736-756, 802-805, 891-894)
- Modify: `frontend/src/components/ProgressNotifications.tsx` (~line 23-26)

**Step 1: Update `onProgressStarted` to not require `roomId`**

The backend no longer sends `roomId` in the socket event (it's room-scoped).
Update the frontend handler at ~line 743:

```typescript
function onProgressStarted(data: any) {
    const { progress_id, description } = data;
    addProgressTracker(progress_id, description, null, roomId);
}
```

Note: We use the `roomId` from the hook's closure (the current room) instead of `data.roomId`.

Also update `onProgressUpdate` (~line 748):
```typescript
function onProgressUpdate(data: any) {
    const { progress_id, description, progress } = data;
    updateProgressTracker(progress_id, description, progress);
}
```

And `onProgressComplete` (~line 753):
```typescript
function onProgressComplete(data: any) {
    const { progress_id } = data;
    removeProgressTracker(progress_id);
}
```

**Step 2: Handle `progress_trackers` from room join response**

In the `onConnect` function (or wherever `room_join` response is handled), add extraction of `progress_trackers` from the join response. Find where `room_join` emits a callback and the response is processed — add:

```typescript
if (data.progress_trackers) {
    setProgressTrackers(data.progress_trackers);
}
```

**Step 3: Remove `progress_init` handler**

Remove these lines (~802, 891):
```typescript
// REMOVE: socket.on("progress_init", onProgressInitial);
// REMOVE: socket.off("progress_init", onProgressInitial);
```

And remove the `onProgressInitial` function (~line 736-741).

**Step 4: Remove `roomId` filter from `ProgressNotifications.tsx`**

Since events are room-scoped, we don't need to filter by `roomId`. Simplify the filter at line 23-26:

```typescript
const progressList = Object.values(progressTrackers).filter(
    (tracker) => !hiddenProgressIds.has(tracker.progressId),
);
```

Also remove the `roomId` import from the store if it's only used here. Check first.

**Step 5: Commit**

```bash
cd frontend && git add src/hooks/useSocketManager.ts src/components/ProgressNotifications.tsx
git commit -m "feat: adapt frontend for progress bar socket events"
```

---

### Task 10: Run Full Test Suite + Format + Lint

**Step 1: Format**

Run: `uv run ruff format .`

**Step 2: Fix imports**

Run: `uv run ruff check --select I --fix .`

**Step 3: Type check**

Run: `uv run pyright .`
Expected: No new errors (Redis async stubs may show known false positives)

**Step 4: Run full test suite**

Run: `uv run pytest tests/ -x`
Expected: ALL PASS

**Step 5: Commit any formatting/lint fixes**

```bash
git add -A && git commit -m "chore: format and lint"
```
