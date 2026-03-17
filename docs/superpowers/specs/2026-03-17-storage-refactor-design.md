# Storage Refactor: Replace Custom Wrappers with Direct AsyncBlobIO

**Date:** 2026-03-17
**Status:** Draft
**Scope:** `src/zndraw/storage/`, `src/zndraw/dependencies.py`, `src/zndraw/database.py`,
all routes, `result_backends.py`, tests

## Problem

The current storage layer wraps `AsyncBlobIO` in two classes (`AsebytesStorage` and
`StorageRouter`) that:

1. **Duplicate AsyncBlobIO's public API** with custom method signatures (`get`, `get_range`,
   `get_many`, `set_item`, `merge_item`, `delete_range`, `reserve`, `remove_items`).
2. **Access `_backend` directly** in 6 places (`noqa: SLF001`) despite public API equivalents
   existing for every operation.
3. **Block selective key loading** — every read loads the full frame, even when a single key
   is needed (e.g., isosurface loads entire frame for one cube key).
4. **Have broken typing** — `StorageDep` is typed as `AsebytesStorage` but the runtime object
   is `StorageRouter`; routes call `StorageRouter`-only methods on the wrong type.
5. **Contain dead code** — `reserve()` and `remove_items()` are never called from routes.
6. **Perform redundant conversion** — `to_raw_frame()` is called by both routes and storage
   methods; the storage call is a no-op because routes already convert.

## Goal

Replace `AsebytesStorage` + `StorageRouter` with a thin `FrameStorage` registry class.
All frame I/O goes through `AsyncBlobIO` public API directly. Zero `_backend` access.
Selective key loading at call sites that benefit from it.

## Design

### Architecture

```
Before:
  Routes -> StorageRouter -> AsebytesStorage -> io._backend (private!)

After:
  Routes -> FrameStorage[room_id] -> AsyncBlobIO public API (reads + writes)
  Routes -> FrameStorage.get_length()  (provider-aware frame count)
  Routes -> RequireWritableDep          (write guard via FastAPI dependency)
```

### `FrameStorage` class

Lives in `storage/frame_storage.py` (replaces both `asebytes_backend.py` class and
`router.py`).

```python
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

    def __init__(self, uri: str, redis: AsyncRedis) -> None:
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
        cached = await self._redis.get(
            RedisKey.provider_frame_count(room_id)
        )
        return int(cached) if cached else 0

    async def has_mount(self, room_id: str) -> bool:
        """Check if a room has a provider-backed frame count."""
        return await self._redis.exists(
            RedisKey.provider_frame_count(room_id)
        ) > 0

    async def set_frame_count(self, room_id: str, count: int) -> None:
        """Store provider frame count in Redis."""
        await self._redis.set(
            RedisKey.provider_frame_count(room_id), count
        )

    async def clear_frame_count(self, room_id: str) -> None:
        """Remove provider frame count from Redis."""
        await self._redis.delete(
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

**What it does NOT have:** `get()`, `get_range()`, `get_many()`, `set_item()`,
`merge_item()`, `extend()`, `delete_range()`, `reserve()`, `remove_items()`.
These are all `AsyncBlobIO` methods — routes call them directly.

### Write guards as FastAPI dependency

Provider-backed rooms are read-only. Instead of checking in every storage write
method, a FastAPI dependency enforces this at the route level:

```python
async def require_writable_room(
    storage: FrameStorageDep,
    room_id: str = Path(),
) -> None:
    """Raise RoomReadOnly if the room has a provider mount."""
    if await storage.has_mount(room_id):
        raise RoomReadOnly.exception("Room is provider-backed (read-only)")

RequireWritableDep = Annotated[None, Depends(require_writable_room)]
```

Added to write endpoints: `append_frames`, `update_frame`, `merge_frame`,
`delete_frame`.

**Note:** `create_room` checks `has_mount(copy_from)` (the *source* room), not
`has_mount(room_id)` (the new room). This cannot use `RequireWritableDep` (which
checks the path's `room_id`). Keep this as an inline `storage.has_mount(copy_from)`
check in the route handler.

### Dependency injection

```python
# Before (broken typing, two deps for same object):
StorageDep = Annotated[AsebytesStorage, Depends(get_storage)]
StorageRouterDep = Annotated[StorageRouter, Depends(get_storage_router)]

# After (one dep, correct type):
FrameStorageDep = Annotated[FrameStorage, Depends(get_frame_storage)]
```

### Initialization (`database.py`)

```python
# Before:
default_storage = AsebytesStorage(uri=settings.storage)
app.state.frame_storage = StorageRouter(default=default_storage, redis=app.state.redis)

# After:
app.state.frame_storage = FrameStorage(uri=settings.storage, redis=app.state.redis)
```

## Call Site Migration

### Reads — use AsyncBlobIO public API

| Before | After |
|--------|-------|
| `storage.get(room_id, index)` | `storage[room_id].get(index)` |
| `storage.get(room_id, index)` then `_filter_frames_by_keys` | `storage[room_id].get(index, keys=[...])` |
| `storage.get_range(room_id, start, stop)` | `await storage[room_id][start:stop].to_list()` |
| `storage.get_many(room_id, indices)` | `await storage[room_id][indices].to_list()` |
| `storage.get_length(room_id)` | `storage.get_length(room_id)` (unchanged) |
| `storage.has_mount(room_id)` | `storage.has_mount(room_id)` (unchanged) |
| `storage.set_frame_count(room_id, n)` | `storage.set_frame_count(room_id, n)` (unchanged) |
| `storage.clear_frame_count(room_id)` | `storage.clear_frame_count(room_id)` (unchanged) |

### Writes — use AsyncBlobIO public API + write guard dependency

| Before | After |
|--------|-------|
| `storage.extend(room_id, frames)` | `await storage[room_id].extend(frames)` |
| `storage.set_item(room_id, index, frame)` | `await storage[room_id][index].set(frame)` |
| `storage.merge_item(room_id, index, partial)` | `await storage[room_id].update(index, partial)` |
| `storage.delete_range(room_id, start, stop)` | `await storage[room_id][start:stop].delete()` |

Write endpoints add `RequireWritableDep` to their signatures. The guard runs before
the handler body, so no write code executes for provider-backed rooms.

### High-impact optimizations

**Isosurface** — biggest win, cube data can be megabytes:
```python
# Before: loads ENTIRE frame
frame = await storage.get(room_id, index)
cube_data = frame[cube_key.encode()]

# After: loads ONLY the requested key
io = storage[room_id]
frame = await io.get(index, keys=[cube_key.encode()])
cube_data = frame[cube_key.encode()]
```

**Single frame with `keys` query param** — selective load on storage-hit path:
```python
# Before: load all, filter in Python
frame = await storage.get(room_id, index)
frames = _filter_frames_by_keys([frame], key_bytes)

# After: selective load at backend (storage-hit path)
io = storage[room_id]
frame = await io.get(index, keys=list(key_bytes))
```

**Note:** When the storage returns `None` and the frame is fetched via provider
dispatch, the provider returns a full frame. Post-load `_filter_frames_by_keys`
is still needed for the provider-fallback path in `get_frame`. The selective
`io.get(keys=...)` optimization applies only when the frame is already in storage.

**Result backend** — full migration of `StorageResultBackend`:
```python
# Before (store):
await self._storage.clear(key)
packed = msgpack.packb(data)
await self._storage.extend(key, [{b"_": packed}])

# After (store):
io = self._storage[key]
await io.clear()
packed = msgpack.packb(data)
await io.extend([{b"_": packed}])

# Before (get):
frame = await self._storage.get(key, 0)
packed = frame.get(b"_")

# After (get):
io = self._storage[key]
frame = await io.get(0, keys=[b"_"])
packed = frame.get(b"_") if frame else None

# Before (delete):
await self._storage.clear(key)

# After (delete):
await self._storage[key].clear()
```

`StorageResultBackend` takes `FrameStorage` but only uses `__getitem__` to get
`AsyncBlobIO` instances. It never calls `get_length`, `has_mount`, or provider
methods. The Redis client in `FrameStorage` is irrelevant for result backend
keys — `clear_frame_count` on non-existent Redis keys is a harmless no-op.

**Note:** Result backend calls `io.clear()` directly (not `storage.clear(room_id)`)
because it does not need provider frame count cleanup.

### Batch reads with key filtering

For `list_frames` with the `keys` query parameter in batch mode, keep post-load
filtering with `_filter_frames_by_keys()`. The batch+keys case is rare and the
Python-side filtering cost is negligible compared to I/O. The big wins are all
single-row selective loads.

### Provider frame padding in `list_frames`

`StorageRouter.get_range()` currently pads results with `None` for provider-backed
rooms. This logic moves into `list_frames` itself — it already handles `None`
entries by dispatching provider requests. The padding becomes:

```python
io = storage[room_id]
total = await storage.get_length(room_id)  # provider-aware count
actual_stop = min(stop, total)
frames = await io[start:actual_stop].to_list()
# The slice resolves against io.len() (storage length), which may be less than
# the provider-declared total from get_length(). E.g. provider declares 100
# frames but only 30 are materialized: io[0:50] returns 30 rows, we pad 20.
if len(frames) < (actual_stop - start) and await storage.has_mount(room_id):
    frames.extend([None] * ((actual_stop - start) - len(frames)))
```

## Files Changed

| File | Change |
|------|--------|
| `storage/frame_storage.py` | **New** — `FrameStorage` class, `_create_blob_backend()`, `to_raw_frame()` |
| `storage/asebytes_backend.py` | **Delete** (or keep only `to_raw_frame` + `_create_blob_backend` if preferred) |
| `storage/router.py` | **Delete** |
| `storage/__init__.py` | Export `FrameStorage`, `RawFrame`, `to_raw_frame` |
| `dependencies.py` | `FrameStorageDep`, `RequireWritableDep`; remove `StorageDep`, `StorageRouterDep` |
| `database.py` | Direct `FrameStorage(uri, redis)` init |
| `routes/frames.py` | Use `storage[room_id]` for reads/writes; selective keys; add `RequireWritableDep` |
| `routes/isosurface.py` | Selective key load |
| `routes/trajectory.py` | Use view API; remove `isinstance(storage, StorageRouter)` check |
| `routes/rooms.py` | Single dep; `storage.has_mount()` directly |
| `routes/step.py` | Trivial — only uses `get_length` |
| `routes/server_settings.py` | Trivial — only passes storage to `build_room_update()` |
| `socketio.py` | Trivial — only uses `get_length` |
| `result_backends.py` | Use `storage[key].get(0, keys=[b"_"])` |
| `tests/` | Update all storage tests to new API |

## What Gets Deleted

- `AsebytesStorage` class
- `StorageRouter` class
- `storage/router.py` file
- `StorageDep` and `StorageRouterDep` type aliases
- `get_storage()` and `get_storage_router()` dependency functions
- `_filter_frames_by_keys()` for single-frame storage-hit reads (kept for provider-fallback
  path and batch+keys case)
- `reserve()` and `remove_items()` (dead code, never called from routes)
- All `noqa: SLF001` suppressions (zero `_backend` access)
- Redundant `to_raw_frame()` calls inside storage methods
- `try: ... except NotImplementedError` around `set_frame_count`/`clear_frame_count` in
  `update_room` (dead code — `FrameStorage` always supports these methods)

## What Stays

- `to_raw_frame()` — standalone utility, called by routes before writes
- `RawFrame` type alias — `dict[bytes, bytes]`
- `_create_blob_backend()` — creates backend from URI
- `_filter_frames_by_keys()` — kept for batch+keys case and provider-fallback path

## Migration Notes

- `storage[room_id]` (`__getitem__`) is **synchronous** — it creates `AsyncBlobIO`
  lazily and can be used in non-async setup code. All I/O methods on the returned
  `AsyncBlobIO` are async.
- Tests using `storage._get_io(room_id)` migrate to `storage[room_id]`.
- `build_room_update()` type annotation changes from `storage: Any` to
  `storage: FrameStorage`.
- `_rooms` dict concurrency: `__getitem__` does check-then-write which is safe in a
  single-event-loop async context (no preemption between check and write).

## Success Criteria

1. All existing tests pass (adapted to new API).
2. Zero `_backend` access / `SLF001` suppressions.
3. `FrameStorage` has no read or write wrapper methods — only registry + provider metadata.
4. Isosurface endpoint loads only the requested cube key, not the full frame.
5. Single-frame `keys` query parameter uses `io.get(keys=...)` instead of post-load filter.
6. Single `FrameStorageDep` replaces both `StorageDep` and `StorageRouterDep`.
7. Write guard enforced via `RequireWritableDep` FastAPI dependency.

## Non-Goals

- Changing the `AsyncBlobIO` or asebytes library.
- Changing the REST API contract (all endpoints keep the same request/response format).
- Optimizing batch+keys reads (kept as post-load filter; minor impact).
- Changing the provider dispatch mechanism.
