# Storage Refactor Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `AsebytesStorage` + `StorageRouter` with thin `FrameStorage` registry; all I/O through AsyncBlobIO pandas-like subscript API.

**Architecture:** `FrameStorage` is a room-scoped `AsyncBlobIO` registry + provider frame count metadata (Redis). No read/write wrappers. Routes call `storage[room_id]` to get an `AsyncBlobIO` and use its subscript API (`io[index]`, `io[start:stop].to_list()`, `io[b"key"][index]`, etc.) directly. Write guards move to a FastAPI dependency (`RequireWritableDep`).

**Tech Stack:** Python 3.12+, FastAPI, asebytes (AsyncBlobIO), Redis, pytest-asyncio

**Spec:** `docs/superpowers/specs/2026-03-17-storage-refactor-design.md`

---

## File Structure

| File | Role | Action |
|------|------|--------|
| `src/zndraw/storage/frame_storage.py` | `FrameStorage` class + `_create_blob_backend` + `to_raw_frame` + `RawFrame` | **Create** |
| `src/zndraw/storage/asebytes_backend.py` | Old `AsebytesStorage` | **Delete** |
| `src/zndraw/storage/router.py` | Old `StorageRouter` | **Delete** |
| `src/zndraw/storage/__init__.py` | Re-exports | **Rewrite** |
| `src/zndraw/dependencies.py` | `FrameStorageDep`, `RequireWritableDep` | **Modify** |
| `src/zndraw/database.py` | Lifespan init | **Modify** |
| `src/zndraw/result_backends.py` | `StorageResultBackend` | **Modify** |
| `src/zndraw/routes/frames.py` | Frame CRUD | **Modify** |
| `src/zndraw/routes/isosurface.py` | Isosurface extraction | **Modify** |
| `src/zndraw/routes/trajectory.py` | Trajectory download/upload | **Modify** |
| `src/zndraw/routes/rooms.py` | Room CRUD | **Modify** |
| `src/zndraw/routes/step.py` | Frame navigation | **Modify** |
| `src/zndraw/routes/server_settings.py` | Server settings | **Modify** |
| `src/zndraw/socketio.py` | Socket.IO handlers | **Modify** |
| `tests/test_storage_asebytes.py` | Unit tests for storage | **Rewrite** as `test_frame_storage.py` |
| `tests/test_storage_router.py` | Unit tests for router | **Rewrite** into `test_frame_storage.py` |
| All other test files | Integration tests | **Modify** fixtures |

---

## Chunk 1: Core — FrameStorage class + dependencies + init

### Task 1: Create `FrameStorage` class

**Files:**
- Create: `src/zndraw/storage/frame_storage.py`

- [ ] **Step 1: Create `frame_storage.py` with `RawFrame`, `to_raw_frame`, `_create_blob_backend`, and `FrameStorage`**

Move `RawFrame`, `to_raw_frame()`, and `_create_blob_backend()` from `asebytes_backend.py` verbatim.
Write `FrameStorage` per the spec — `__getitem__`, `get_length`, `has_mount`, `set_frame_count`, `clear_frame_count`, `clear`, `close`.

```python
"""Room-scoped AsyncBlobIO registry with provider frame count support."""

import base64
from typing import Any

import msgpack
from asebytes import AsyncBlobIO
from asebytes._async_adapters import AsyncObjectToBlobReadWriteAdapter
from asebytes._async_backends import AsyncReadWriteBackend
from redis.asyncio import Redis as AsyncRedis

from zndraw.redis import RedisKey

# Type alias for raw frame data: dict with bytes keys and bytes values
# Keys are like b"arrays.positions", b"cell", etc.
# Values are msgpack-numpy encoded bytes
RawFrame = dict[bytes, bytes]


def to_raw_frame(frame: dict[str, Any] | RawFrame) -> RawFrame:
    """Convert input frame dict to raw bytes format.

    Handles two input formats:
    1. Already raw: dict[bytes, bytes] -- pass through
    2. Base64 encoded: dict[str, str] with keys like "b64:..." -- decode

    Parameters
    ----------
    frame
        Input frame data (dict[str, Any] or dict[bytes, bytes]).
    """
    if not frame:
        return {}

    first_key = next(iter(frame.keys()))

    if isinstance(first_key, bytes):
        return frame  # type: ignore[return-value]

    result: RawFrame = {}
    for k, v in frame.items():
        key_bytes: bytes
        if isinstance(k, str) and k.startswith("b64:"):
            key_bytes = base64.b64decode(k[4:])
        elif isinstance(k, str):
            key_bytes = k.encode()
        elif isinstance(k, bytes):
            key_bytes = k
        else:
            key_bytes = str(k).encode()

        val_bytes: bytes
        if isinstance(v, str):
            val_bytes = base64.b64decode(v)
        elif isinstance(v, bytes):
            val_bytes = v
        else:
            packed = msgpack.packb(v)
            val_bytes = packed if packed is not None else b""

        result[key_bytes] = val_bytes

    return result


def _create_blob_backend(
    uri: str, group: str | None = None
) -> AsyncReadWriteBackend[bytes, bytes]:
    """Create an async blob backend from a URI string.

    All backends from asebytes registries are object-level (str keys).
    This function wraps them to blob-level (bytes keys) via
    AsyncObjectToBlobReadWriteAdapter, converting sync backends to
    async first when needed.

    Parameters
    ----------
    uri
        Storage backend URI string.
    group
        Optional group/namespace for data isolation (e.g., room_id).
    """
    from asebytes._registry import get_async_backend_cls

    cls = get_async_backend_cls(uri, readonly=False)
    if hasattr(cls, "from_uri"):
        backend = cls.from_uri(uri, group=group)
    else:
        backend = cls(uri, group=group)

    # Ensure async
    if not isinstance(backend, AsyncReadWriteBackend):
        from asebytes._async_backends import sync_to_async

        backend = sync_to_async(backend)  # type: ignore[arg-type]

    # All registry backends are object-level (str keys) -> wrap to blob
    return AsyncObjectToBlobReadWriteAdapter(backend)  # type: ignore[arg-type]


class FrameStorage:
    """Room-scoped AsyncBlobIO registry with provider frame count support.

    No I/O wrappers. Routes use AsyncBlobIO directly via ``storage[room_id]``.

    Parameters
    ----------
    uri
        Storage backend URI (``memory://``, ``path.lmdb``, ``mongodb://...``).
    redis
        Async Redis client for provider frame count metadata.
    """

    def __init__(self, uri: str, redis: AsyncRedis) -> None:  # type: ignore[type-arg]
        self._uri = uri
        self._redis = redis
        self._rooms: dict[str, AsyncBlobIO] = {}

    def __getitem__(self, room_id: str) -> AsyncBlobIO:
        """Get or create the AsyncBlobIO for a room."""
        if room_id not in self._rooms:
            backend = _create_blob_backend(self._uri, group=room_id)
            self._rooms[room_id] = AsyncBlobIO(backend)
        return self._rooms[room_id]

    # -- Provider frame count (Redis) -----------------------------------------

    async def get_length(self, room_id: str) -> int:
        """Frame count: storage first, then Redis provider fallback."""
        io = self[room_id]
        length = await io.len()
        if length > 0:
            return length
        cached = await self._redis.get(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id)
        )
        return int(cached) if cached else 0

    async def has_mount(self, room_id: str) -> bool:
        """Check if a room has a provider-backed frame count."""
        return await self._redis.exists(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id)
        ) > 0

    async def set_frame_count(self, room_id: str, count: int) -> None:
        """Store provider frame count in Redis."""
        await self._redis.set(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id), count
        )

    async def clear_frame_count(self, room_id: str) -> None:
        """Remove provider frame count from Redis."""
        await self._redis.delete(  # type: ignore[misc]
            RedisKey.provider_frame_count(room_id)
        )

    # -- Lifecycle -------------------------------------------------------------

    async def clear(self, room_id: str) -> None:
        """Clear all frames AND provider frame count for a room."""
        io = self[room_id]
        await io.clear()
        await self.clear_frame_count(room_id)

    async def close(self) -> None:
        """Release in-memory handles."""
        self._rooms.clear()
```

- [ ] **Step 2: Verify file has no syntax errors**

Run: `uv run python -c "from zndraw.storage.frame_storage import FrameStorage, RawFrame, to_raw_frame"`
Expected: No output (clean import)

- [ ] **Step 3: Commit**

```bash
git add src/zndraw/storage/frame_storage.py
git commit -m "feat: add FrameStorage class (room registry + provider metadata)"
```

---

### Task 2: Update `storage/__init__.py` and dependencies

**Files:**
- Modify: `src/zndraw/storage/__init__.py`
- Modify: `src/zndraw/dependencies.py`

- [ ] **Step 1: Rewrite `storage/__init__.py`**

```python
"""Frame storage backends."""

from .frame_storage import FrameStorage, RawFrame, to_raw_frame

__all__ = ["FrameStorage", "RawFrame", "to_raw_frame"]
```

- [ ] **Step 2: Update `dependencies.py`**

Replace `StorageDep`, `StorageRouterDep`, `get_storage`, `get_storage_router` with `FrameStorageDep`, `RequireWritableDep`, `get_frame_storage`:

- Remove imports of `AsebytesStorage` and `StorageRouter`
- Add import of `FrameStorage`
- Replace `get_storage()` with `get_frame_storage()` returning `FrameStorage`
- Delete `get_storage_router()` entirely
- Add `require_writable_room()` dependency
- Add `RequireWritableDep` type alias
- Import `RoomReadOnly` from exceptions and `Path` from fastapi

The new dependency section (replacing lines 32-34 and 56-72):

```python
from zndraw.storage import FrameStorage
# ... (remove AsebytesStorage and StorageRouter imports)

def get_frame_storage(request: Request) -> FrameStorage:
    """Get frame storage registry from app.state."""
    return request.app.state.frame_storage


FrameStorageDep = Annotated[FrameStorage, Depends(get_frame_storage)]


async def require_writable_room(
    storage: FrameStorageDep,
    room_id: str = Path(),
) -> None:
    """Raise RoomReadOnly if the room has a provider mount."""
    if await storage.has_mount(room_id):
        raise RoomReadOnly.exception("Room is provider-backed (read-only)")


RequireWritableDep = Annotated[None, Depends(require_writable_room)]
```

- [ ] **Step 3: Verify imports work**

Run: `uv run python -c "from zndraw.dependencies import FrameStorageDep, RequireWritableDep"`
Expected: No output (clean import)

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/storage/__init__.py src/zndraw/dependencies.py
git commit -m "feat: replace StorageDep/StorageRouterDep with FrameStorageDep + RequireWritableDep"
```

---

### Task 3: Update `database.py` lifespan init

**Files:**
- Modify: `src/zndraw/database.py`

- [ ] **Step 1: Update lifespan storage initialization**

Find the section (around lines 239-248) that creates `AsebytesStorage` and `StorageRouter`. Replace with:

```python
# Before:
from zndraw.storage.router import StorageRouter
default_storage = AsebytesStorage(uri=settings.storage)
app.state.frame_storage = StorageRouter(
    default=default_storage,
    redis=app.state.redis,
)

# After:
from zndraw.storage import FrameStorage
app.state.frame_storage = FrameStorage(
    uri=settings.storage,
    redis=app.state.redis,
)
```

Also find where `default_storage` is passed to `StorageResultBackend` (around line 320) and update to pass `app.state.frame_storage` instead.

- [ ] **Step 2: Verify app starts without errors**

Run: `uv run python -c "from zndraw.database import create_app; create_app()"`
Expected: No import errors (app won't fully start without Redis, but the import should work)

- [ ] **Step 3: Commit**

```bash
git add src/zndraw/database.py
git commit -m "feat: init FrameStorage directly in lifespan (replaces AsebytesStorage + StorageRouter)"
```

---

### Task 4: Update `result_backends.py`

**Files:**
- Modify: `src/zndraw/result_backends.py`

- [ ] **Step 1: Update `StorageResultBackend` to use `FrameStorage` + subscript API**

Change the type annotation and all method bodies:

```python
# TYPE_CHECKING import changes:
# Before:
from zndraw.storage import AsebytesStorage
# After:
from zndraw.storage import FrameStorage

# Class changes:
class StorageResultBackend:
    """Adapt ``FrameStorage`` to the ``ResultBackend`` protocol.

    Uses cache keys as room_ids and stores raw bytes as a single-entry
    frame (``{b"_": data}``).  Uses column view ``io[b"_"][0]`` for
    direct scalar access.
    """

    def __init__(self, storage: FrameStorage) -> None:
        self._storage = storage

    async def store(self, key: str, data: bytes, ttl: int) -> None:  # noqa: ARG002
        io = self._storage[key]
        await io.clear()
        packed = msgpack.packb(data)
        assert packed is not None
        await io.extend([{b"_": packed}])

    async def get(self, key: str) -> bytes | None:
        io = self._storage[key]
        if await io.len() == 0:
            return None
        packed = await io[b"_"][0]
        if packed is None:
            return None
        return msgpack.unpackb(packed)

    async def delete(self, key: str) -> None:
        await self._storage[key].clear()

    # ... (rest unchanged)
```

- [ ] **Step 2: Commit**

```bash
git add src/zndraw/result_backends.py
git commit -m "refactor: StorageResultBackend uses FrameStorage + io[b'_'][0] column view"
```

---

### Task 5: Delete old storage files

**Files:**
- Delete: `src/zndraw/storage/asebytes_backend.py`
- Delete: `src/zndraw/storage/router.py`

- [ ] **Step 1: Delete old files**

```bash
git rm src/zndraw/storage/asebytes_backend.py src/zndraw/storage/router.py
```

- [ ] **Step 2: Verify no dangling imports**

Run: `uv run python -c "from zndraw.storage import FrameStorage, RawFrame, to_raw_frame; from zndraw.dependencies import FrameStorageDep"`
Expected: Clean import

- [ ] **Step 3: Commit**

```bash
git commit -m "chore: delete AsebytesStorage and StorageRouter (replaced by FrameStorage)"
```

---

## Chunk 2: Route Migration — frames, isosurface, trajectory

### Task 6: Migrate `routes/frames.py`

**Files:**
- Modify: `src/zndraw/routes/frames.py`

This is the largest migration. Changes:

- [ ] **Step 1: Update imports**

```python
# Before:
from zndraw.dependencies import StorageDep, ...
# After:
from zndraw.dependencies import FrameStorageDep, RequireWritableDep, ...
```

- [ ] **Step 2: Update `list_frames` (lines ~180-279)**

Replace `storage: StorageDep` with `storage: FrameStorageDep` in signature.

Replace read calls:
- `storage.get_length(room_id)` stays as-is (method exists on `FrameStorage`)
- `storage.get_many(room_id, requested_indices)` becomes `await storage[room_id][requested_indices].to_list()`
- `storage.get_range(room_id, effective_start, effective_stop)` becomes `await storage[room_id][effective_start:effective_stop].to_list()`

Add provider padding after range/many reads:
```python
# After fetching frames via subscript:
if len(frames_or_none) < expected and await storage.has_mount(room_id):
    frames_or_none.extend([None] * (expected - len(frames_or_none)))
```

- [ ] **Step 3: Update `get_frame` (lines ~287-338)**

Replace `storage: StorageDep` with `storage: FrameStorageDep`.

Storage-hit path with keys:
```python
# Before:
frame = await storage.get(room_id, index)
# ...
if keys:
    key_bytes = frozenset(k.strip().encode() for k in keys.split(","))
    frames = _filter_frames_by_keys(frames, key_bytes)

# After:
io = storage[room_id]
if keys:
    key_bytes_list = [k.strip().encode() for k in keys.split(",")]
    frame = await io.get(index, keys=key_bytes_list)
else:
    frame = await io[index]
# ... provider fallback path still uses _filter_frames_by_keys on provider frame
```

- [ ] **Step 4: Update `get_frame_metadata` (lines ~373-424)**

Replace `storage: StorageDep` with `storage: FrameStorageDep`.

```python
# Before:
frame = await storage.get(room_id, index)

# After:
frame = await storage[room_id][index]
```

- [ ] **Step 5: Update `append_frames` (lines ~430-470)**

Replace `storage: StorageDep` with `storage: FrameStorageDep`. Add `_: RequireWritableDep` parameter.

```python
# Before:
old_total = await storage.get_length(room_id)
raw_frames = [to_raw_frame(f) for f in request.frames]
new_total = await storage.extend(room_id, raw_frames)

# After:
old_total = await storage.get_length(room_id)
raw_frames = [to_raw_frame(f) for f in request.frames]
io = storage[room_id]
new_total = await io.extend(raw_frames)
```

- [ ] **Step 6: Update `update_frame` (lines ~479-510)**

Replace `storage: StorageDep` with `storage: FrameStorageDep`. Add `_: RequireWritableDep`.

```python
# Before:
await storage.set_item(room_id, index, raw_frame)

# After:
await storage[room_id][index].set(raw_frame)
```

- [ ] **Step 7: Update `merge_frame` (lines ~514-560)**

Replace `storage: StorageDep` with `storage: FrameStorageDep`. Add `_: RequireWritableDep`.

```python
# Before:
await storage.merge_item(room_id, index, partial)

# After:
await storage[room_id][index].update(partial)
```

- [ ] **Step 8: Update `delete_frame` (lines ~565-590)**

Replace `storage: StorageDep` with `storage: FrameStorageDep`. Add `_: RequireWritableDep`.

```python
# Before:
await storage.delete_range(room_id, index, index + 1)

# After:
await storage[room_id][index:index + 1].delete()
```

- [ ] **Step 9: Commit**

```bash
git add src/zndraw/routes/frames.py
git commit -m "refactor: frames routes use FrameStorage subscript API"
```

---

### Task 7: Migrate `routes/isosurface.py`

**Files:**
- Modify: `src/zndraw/routes/isosurface.py`

- [ ] **Step 1: Update imports and endpoint**

Replace `StorageDep` with `FrameStorageDep`.

Replace frame loading with column view for selective key load:

```python
# Before:
frame = await storage.get(room_id, index)
if frame is None:
    raise FrameNotFound.exception(...)
key_bytes = cube_key.encode()
if key_bytes not in frame:
    raise UnprocessableContent.exception(...)
cube_dict = msgpack.unpackb(frame[key_bytes], ...)

# After:
io = storage[room_id]
key_bytes = cube_key.encode()
try:
    cube_raw = await io[key_bytes][index]
except (IndexError, KeyError):
    cube_raw = None
if cube_raw is None:
    raise UnprocessableContent.exception(
        f"Key '{cube_key}' not found in frame {index}"
    )
cube_dict = msgpack.unpackb(cube_raw, ...)
```

Note: We still need to check that the frame exists first (for proper FrameNotFound vs UnprocessableContent errors). Read the full endpoint to determine the exact error flow.

- [ ] **Step 2: Commit**

```bash
git add src/zndraw/routes/isosurface.py
git commit -m "refactor: isosurface uses io[cube_key][index] column view for selective load"
```

---

### Task 8: Migrate `routes/trajectory.py`

**Files:**
- Modify: `src/zndraw/routes/trajectory.py`

- [ ] **Step 1: Update imports**

Remove `StorageRouter` import. Replace `StorageDep` with `FrameStorageDep`.

- [ ] **Step 2: Update download endpoint**

```python
# Before:
if isinstance(storage, StorageRouter) and await storage.has_mount(room_id):
    raise InvalidPayload.exception(...)

# After:
if await storage.has_mount(room_id):
    raise InvalidPayload.exception(...)
```

```python
# Before:
first_frame = await storage.get(room_id, index_list[0])

# After:
first_frame = await storage[room_id][index_list[0]]
```

```python
# Before:
raw_frames = await storage.get_many(room_id, batch_indices)

# After:
raw_frames = await storage[room_id][batch_indices].to_list()
```

- [ ] **Step 3: Update upload endpoint**

```python
# Before:
old_total = await storage.get_length(room_id)
new_total = await storage.extend(room_id, frames)

# After:
old_total = await storage.get_length(room_id)
new_total = await storage[room_id].extend(frames)
```

- [ ] **Step 4: Commit**

```bash
git add src/zndraw/routes/trajectory.py
git commit -m "refactor: trajectory routes use FrameStorage subscript API"
```

---

## Chunk 3: Route Migration — rooms, step, server_settings, socketio

### Task 9: Migrate `routes/rooms.py`

**Files:**
- Modify: `src/zndraw/routes/rooms.py`

- [ ] **Step 1: Update imports**

Remove `StorageRouterDep` import. Replace `StorageDep` with `FrameStorageDep`.

- [ ] **Step 2: Update `build_room_update` helper**

Change type annotation from `storage: Any` (or whatever it is) to `storage: FrameStorage`.

- [ ] **Step 3: Update `create_room`**

Remove `storage_router: StorageRouterDep` parameter. Replace with single `storage: FrameStorageDep`.

```python
# Before:
if source_room is not None and await storage_router.has_mount(copy_from):
    raise RoomReadOnly.exception(...)
source_frames_or_none = await storage.get_range(copy_from, 0, None)
source_frames = [f for f in source_frames_or_none if f is not None]
if source_frames:
    await storage.extend(room_id, source_frames)

# After:
if source_room is not None and await storage.has_mount(copy_from):
    raise RoomReadOnly.exception(...)
source_frames_or_none = await storage[copy_from][0:].to_list()
source_frames = [f for f in source_frames_or_none if f is not None]
if source_frames:
    await storage[room_id].extend(source_frames)
```

Also update `@empty` and fallback paths:
```python
# Before:
await storage.extend(room_id, [{}])

# After:
await storage[room_id].extend([{}])
```

- [ ] **Step 4: Update `update_room`**

Remove `try/except NotImplementedError` around `set_frame_count`/`clear_frame_count`:

```python
# Before:
try:
    if count > 0:
        await storage.set_frame_count(room.id, count)
    else:
        await storage.clear_frame_count(room.id)
except NotImplementedError as err:
    raise Forbidden.exception(...)

# After:
if count > 0:
    await storage.set_frame_count(room.id, count)
else:
    await storage.clear_frame_count(room.id)
```

- [ ] **Step 5: Update all other `StorageDep` references to `FrameStorageDep`**

Every function with `storage: StorageDep` changes to `storage: FrameStorageDep`.
`storage.get_length(...)` calls stay as-is.

- [ ] **Step 6: Commit**

```bash
git add src/zndraw/routes/rooms.py
git commit -m "refactor: rooms routes use FrameStorage (single dep, subscript API)"
```

---

### Task 10: Migrate `routes/step.py`, `routes/server_settings.py`, `socketio.py`

**Files:**
- Modify: `src/zndraw/routes/step.py`
- Modify: `src/zndraw/routes/server_settings.py`
- Modify: `src/zndraw/socketio.py`

These are trivial — only `get_length` calls. Just swap `StorageDep` → `FrameStorageDep` in imports and signatures.

- [ ] **Step 1: Update `step.py`**

Replace import and both function signatures. All `storage.get_length(room_id)` calls stay as-is.

- [ ] **Step 2: Update `server_settings.py`**

Replace import and function signatures. The `storage` param passed to `build_room_update()` is already handled.

- [ ] **Step 3: Update `socketio.py`**

Replace `StorageDep` import with `FrameStorageDep`. Update `room_join` handler signature.

- [ ] **Step 4: Verify all routes import cleanly**

Run: `uv run python -c "from zndraw.routes import frames, isosurface, trajectory, rooms, step, server_settings; from zndraw.socketio import register_events"`
Expected: No import errors

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/routes/step.py src/zndraw/routes/server_settings.py src/zndraw/socketio.py
git commit -m "refactor: step, server_settings, socketio use FrameStorageDep"
```

---

## Chunk 4: Test Migration

### Task 11: Rewrite `test_storage_asebytes.py` as `test_frame_storage.py`

**Files:**
- Create: `tests/test_frame_storage.py`
- Delete: `tests/test_storage_asebytes.py`

- [ ] **Step 1: Create `test_frame_storage.py`**

Port all tests from `test_storage_asebytes.py`, changing:
- `AsebytesStorage("memory://")` → need `FrameStorage("memory://", redis)` but `FrameStorage` requires Redis. For unit tests that don't need provider features, create a mock Redis or use a real Redis fixture.
- `storage._get_io(room)` → `storage[room]`
- `storage.get(room, index)` → `await storage[room][index]`
- `storage.get_range(room, start, stop)` → `await storage[room][start:stop].to_list()`
- `storage.get_many(room, indices)` → `await storage[room][indices].to_list()`
- `storage.extend(room, frames)` → `await storage[room].extend(frames)`
- `storage.set_item(room, index, frame)` → `await storage[room][index].set(frame)`
- `storage.merge_item(room, index, partial)` → `await storage[room][index].update(partial)`
- `storage.delete_range(room, start, stop)` → `await storage[room][start:stop].delete()`
- `storage.clear(room)` → `await storage.clear(room)`
- `storage.get_length(room)` → `await storage.get_length(room)`
- Remove `test_reserve` and `test_remove_items` (dead code)

- [ ] **Step 2: Run tests**

Run: `uv run pytest tests/test_frame_storage.py -v`
Expected: All tests pass

- [ ] **Step 3: Delete old test file**

```bash
git rm tests/test_storage_asebytes.py
```

- [ ] **Step 4: Commit**

```bash
git add tests/test_frame_storage.py
git commit -m "test: rewrite storage tests for FrameStorage subscript API"
```

---

### Task 12: Rewrite `test_storage_router.py`

**Files:**
- Modify or rewrite: `tests/test_storage_router.py`

Port all tests that test provider mount functionality. These now test `FrameStorage` directly:
- `test_set_frame_count_makes_room_mounted` → test `storage.has_mount()` after `storage.set_frame_count()`
- `test_get_length_returns_provider_count` → test `storage.get_length()` fallback
- `test_clear_frame_count_removes_mount` → test `storage.clear_frame_count()`
- `test_clear_removes_provider_count` → test `storage.clear()` removes Redis key
- Provider padding tests move to route-level tests (since padding logic moved to `list_frames`)
- Write-guard tests (extend/set/delete raises on provider room) → these now test `RequireWritableDep` in route tests, not storage

- [ ] **Step 1: Rewrite tests**

- [ ] **Step 2: Run tests**

Run: `uv run pytest tests/test_storage_router.py -v`
Expected: All tests pass

- [ ] **Step 3: Commit**

```bash
git add tests/test_storage_router.py
git commit -m "test: rewrite router tests for FrameStorage provider metadata"
```

---

### Task 13: Update integration test fixtures

**Files:**
- Modify: `tests/test_routes_frames.py`
- Modify: `tests/test_trajectory.py`
- Modify: `tests/test_isosurface.py`
- Modify: `tests/test_routes_step.py`
- Modify: `tests/test_frames_provider_dispatch.py`
- Modify: `tests/test_result_backends.py`
- Modify: `tests/test_constraints.py`
- Modify: `tests/test_lifespan.py`

All integration tests create `AsebytesStorage("memory://")` in fixtures and override `get_storage`.
Change to `FrameStorage("memory://", redis)` and override `get_frame_storage`.

- [ ] **Step 1: Update each test file's fixture**

Pattern for each file:
```python
# Before:
from zndraw.storage import AsebytesStorage
from zndraw.dependencies import get_storage

storage = AsebytesStorage("memory://")
app.dependency_overrides[get_storage] = lambda: storage

# After:
from zndraw.storage import FrameStorage
from zndraw.dependencies import get_frame_storage

# For tests that don't need Redis (most integration tests):
# Use a mock Redis or the redis fixture
storage = FrameStorage("memory://", redis_client)
app.dependency_overrides[get_frame_storage] = lambda: storage
```

For tests without Redis fixture, create a minimal mock:
```python
from unittest.mock import AsyncMock
mock_redis = AsyncMock()
mock_redis.get = AsyncMock(return_value=None)
mock_redis.exists = AsyncMock(return_value=0)
storage = FrameStorage("memory://", mock_redis)
```

- [ ] **Step 2: Update storage method calls in tests**

Tests that call storage methods directly (e.g., `frame_storage.extend(room.id, ...)`) need to use subscript API:
```python
# Before:
await frame_storage.extend(room.id, [{"a": 1}])
await frame_storage.get(room.id, 0)
await frame_storage.get_length(room.id)

# After:
await frame_storage[room.id].extend([{"a": 1}])
await frame_storage[room.id][0]
await frame_storage.get_length(room.id)  # this stays
```

- [ ] **Step 3: Update `test_frames_provider_dispatch.py`**

This file uses `reserve()` to simulate provider mounts. Replace with `set_frame_count()`:
```python
# Before:
await prov_storage.reserve(room.id, 3)

# After:
await prov_storage.set_frame_count(room.id, 3)
```

- [ ] **Step 4: Update `test_result_backends.py`**

Change fixture to use `FrameStorage`:
```python
# Before:
storage_backend = AsebytesStorage("memory://")
backend = StorageResultBackend(storage_backend)

# After:
storage = FrameStorage("memory://", mock_redis)
backend = StorageResultBackend(storage)
```

- [ ] **Step 5: Update `test_constraints.py`**

Tests use `AsebytesStorage` directly for storage roundtrip testing. Either:
- Change to `FrameStorage` with mock Redis, or
- Use `AsyncBlobIO` directly (since these tests don't need room routing)

- [ ] **Step 6: Update `test_lifespan.py`**

```python
# Before:
backend = AsebytesStorage(uri="memory://")
assert isinstance(backend, AsebytesStorage)

# After:
# Test FrameStorage creation or just remove this test
# (it tests a trivial constructor)
```

- [ ] **Step 7: Run full test suite**

Run: `uv run pytest tests/ -v --timeout=900`
Expected: All tests pass

- [ ] **Step 8: Commit**

```bash
git add tests/
git commit -m "test: migrate all integration tests to FrameStorage + subscript API"
```

---

## Chunk 5: Cleanup and Verification

### Task 14: Final verification

- [ ] **Step 1: Verify zero SLF001 violations**

Run: `uv run ruff check src/ --select SLF001`
Expected: No violations (or only pre-existing ones unrelated to storage)

- [ ] **Step 2: Verify no remaining references to old classes**

```bash
grep -r "AsebytesStorage\|StorageRouter\|StorageDep\|StorageRouterDep\|get_storage_router\|_get_io\|_backend\." src/zndraw/ --include="*.py" | grep -v __pycache__
```
Expected: No matches

- [ ] **Step 3: Run type checker**

Run: `uv run pyright .`
Expected: No new errors related to storage types

- [ ] **Step 4: Run formatter and import sorter**

```bash
uv run ruff format .
uv run ruff check --select I --fix .
```

- [ ] **Step 5: Run full test suite one more time**

Run: `uv run pytest tests/ -v --timeout=900`
Expected: All tests pass

- [ ] **Step 6: Final commit (if any formatting changes)**

```bash
git add -A
git commit -m "chore: format and lint cleanup after storage refactor"
```
