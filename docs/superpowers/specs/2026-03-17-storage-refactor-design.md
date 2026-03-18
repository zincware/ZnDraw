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
5. **Contain dead code** — `reserve()` and `remove_items()` are never called from
   production code. Provider mounts use `set_frame_count()` in Redis, not
   pre-allocated storage slots. Test-only usage simulates provider scenarios.
6. **Perform redundant conversion** — `to_raw_frame()` is called by both routes and storage
   methods; the storage call is a no-op because routes already convert.

## Goal

Replace `AsebytesStorage` + `StorageRouter` with a thin `FrameStorage` registry class.
All frame I/O goes through `AsyncBlobIO`'s **pandas-like subscript API** directly.
Zero `_backend` access. Selective key loading at call sites that benefit from it.

### Design Principle: Prefer the Subscript API

`AsyncBlobIO` provides a pandas/xarray-style `__getitem__` API that returns lazy
views. **This is the idiomatic asebytes interface and MUST be preferred throughout
the codebase.** The subscript API is expressive, chainable, and maps directly to
efficient backend operations.

```python
# PREFERRED — subscript/view API (pandas-like)
frame   = await io[index]                   # single row
frames  = await io[start:stop].to_list()    # range of rows
frames  = await io[[0, 5, 9]].to_list()     # sparse rows
value   = await io[b"cube_key"][index]       # single column value
values  = await io[b"energy"].to_list()      # full column
await io[index].set(data)                    # write row
await io[start:stop].delete()                # delete range
await io[index].update(partial)              # merge into row

# ACCEPTABLE — .get() when you need a filtered dict back
frame = await io.get(index, keys=[...])      # partial row as dict

# FORBIDDEN — private backend access
rows = await io._backend.get_many(indices)   # NEVER
```

**When to use `.get(index, keys=[...])`:** Only when the caller needs a
`dict[bytes, bytes]` with a subset of keys (e.g., the `?keys=` query parameter
on frame endpoints, where the response format is a frame dict). For single-value
extraction, use the column view: `await io[b"key"][index]`.

## AsyncBlobIO Subscript API Reference

`AsyncBlobIO` uses a pandas-like `__getitem__` API. Subscript is always **sync**
and returns a lazy view. The backend is only hit on `await` / `.to_list()` /
`async for`.

### Row access (subscript with int, slice, or list[int])

| Syntax | View type | Materialized result |
|--------|-----------|-------------------|
| `io[i]` | `AsyncSingleRowView` | `await` -> `dict[bytes, bytes] \| None` |
| `io[a:b]` | `AsyncRowView` | `.to_list()` -> `list[dict \| None]` |
| `io[[0, 5, 9]]` | `AsyncRowView` | `.to_list()` -> `list[dict \| None]` |

### Column access — single key (subscript with `bytes`)

| Syntax | View type | Materialized result |
|--------|-----------|-------------------|
| `io[b"key"]` | `AsyncColumnView` | `.to_list()` -> `list[bytes]` |
| `io[b"key"][i]` | `AsyncSingleColumnView` | `await` -> `bytes` (scalar) |
| `io[b"key"][a:b]` | `AsyncColumnView` | `.to_list()` -> `list[bytes]` |

### Multi-key access (subscript with `list[bytes]`)

| Syntax | View type | Materialized result |
|--------|-----------|-------------------|
| `io[[b"k1", b"k2"]][i]` | `AsyncSingleColumnView` | `await` -> `[val1, val2]` |
| `io[[b"k1", b"k2"]][a:b]` | `AsyncColumnView` | `.to_list()` -> `[[v1, v2], ...]` |
| `io[[b"k1", b"k2"]]` | `AsyncColumnView` | `.to_dict()` -> `{b"k1": [...], b"k2": [...]}` |

Multi-key subscript calls `_read_rows(indices, keys=...)` on the backend —
**only the requested keys are loaded from storage**. This is the most efficient
way to load a known subset of keys across multiple rows.

`.to_dict()` returns a **column-oriented** dict (key -> list of values across rows).
`.to_list()` returns **row-oriented** lists of values (one list per row).

### Mutations (on views)

| Syntax | Effect |
|--------|--------|
| `await io[i].set(data)` | Replace row at index |
| `await io[i].set(None)` | Soft-delete (placeholder) |
| `await io[i].update(partial)` | Merge partial dict into row |
| `await io[a:b].delete()` | Delete range with index shifting |
| `await io.update(i, partial)` | Same as `io[i].update(partial)` |
| `await io.extend(frames)` | Append rows |

### Efficiency by access pattern

| Pattern | Best approach | Backend calls |
|---------|--------------|---------------|
| 1 key, 1 row | `await io[b"key"][i]` | 1 (`get_column`) |
| N keys, 1 row (values) | `await io[[b"k1", b"k2"]][i]` | 1 (`get`) |
| N keys, 1 row (as dict) | `await io.get(i, keys=[...])` | 1 (`get`) |
| 1 key, N rows | `await io[b"key"].to_list()` | 1 (`get_column`) |
| N keys, N rows | `await io[[b"k1", b"k2"]][a:b].to_list()` | 1 (`get_many`) |
| N rows (full) | `await io[a:b].to_list()` | 1 (`get_many`) |
| Sparse rows | `await io[[0, 5, 9]].to_list()` | 1 (`get_many`) |

### When to use which

- **`io[b"key"][i]`** — extracting a single value (isosurface cube data, result backend)
- **`io[[b"k1", b"k2"]][i]`** — extracting multiple known values from one row
- **`io[[b"k1", b"k2"]][a:b].to_list()`** — sparse batch load of specific keys
- **`io.get(i, keys=[...])`** — only when caller needs a `dict[bytes, bytes]` back
  (e.g., `?keys=` query param where response format is a frame dict)
- **`io[a:b].to_list()`** — full row batch (rendering, trajectory export)

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

### Reads — subscript API

| Before | After (subscript) |
|--------|-------------------|
| `storage.get(room_id, index)` | `await storage[room_id][index]` |
| `storage.get(room_id, index)` + extract 1 key | `await storage[room_id][b"key"][index]` |
| `storage.get(room_id, index)` + filter N keys as dict | `await storage[room_id].get(index, keys=[...])` |
| `storage.get_range(room_id, start, stop)` | `await storage[room_id][start:stop].to_list()` |
| `storage.get_many(room_id, indices)` | `await storage[room_id][indices].to_list()` |
| `storage.get_length(room_id)` | `storage.get_length(room_id)` (unchanged) |
| `storage.has_mount(room_id)` | `storage.has_mount(room_id)` (unchanged) |
| `storage.set_frame_count(room_id, n)` | `storage.set_frame_count(room_id, n)` (unchanged) |
| `storage.clear_frame_count(room_id)` | `storage.clear_frame_count(room_id)` (unchanged) |

### Writes — subscript API + write guard dependency

| Before | After (subscript) |
|--------|-------------------|
| `storage.extend(room_id, frames)` | `await storage[room_id].extend(frames)` |
| `storage.set_item(room_id, index, frame)` | `await storage[room_id][index].set(frame)` |
| `storage.merge_item(room_id, index, partial)` | `await storage[room_id][index].update(partial)` |
| `storage.delete_range(room_id, start, stop)` | `await storage[room_id][start:stop].delete()` |

Write endpoints add `RequireWritableDep` to their signatures. The guard runs before
the handler body, so no write code executes for provider-backed rooms.

### High-impact optimizations

**Isosurface** — biggest win, cube data can be megabytes:
```python
# Before: loads ENTIRE frame for one key
frame = await storage.get(room_id, index)
cube_data = msgpack.unpackb(frame[cube_key.encode()], ...)

# After: column view extracts the single value directly
cube_raw = await storage[room_id][cube_key.encode()][index]
cube_data = msgpack.unpackb(cube_raw, ...)
```

**Single frame with `keys` query param** — selective load on storage-hit path:
```python
# Before: load all, filter in Python
frame = await storage.get(room_id, index)
frames = _filter_frames_by_keys([frame], key_bytes)

# After: filtered get returns partial dict
frame = await storage[room_id].get(index, keys=list(key_bytes))
```

**Trajectory validation** — only needs atom count, not full frame:
```python
# Before: loads entire frame to count atoms
first_frame = await storage.get(room_id, index_list[0])
n_atoms = len(decode(first_frame))

# After: subscript row view
first_frame = await storage[room_id][index_list[0]]
n_atoms = len(decode(first_frame))
```

**Note:** When the storage returns `None` and the frame is fetched via provider
dispatch, the provider returns a full frame (it unconditionally enriches with
colors, radii, connectivity before caching). Passing `keys` through the dispatch
chain is not worth it: it would create separate cache entries per key subset,
and the provider still loads the full source object. Post-load
`_filter_frames_by_keys` is kept for the provider-fallback path.

**Result backend** — use column view API for direct scalar access:
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

# After (get) — column view returns scalar directly:
packed = await self._storage[key][b"_"][0]

# Before (delete):
await self._storage.clear(key)

# After (delete):
await self._storage[key].clear()
```

Uses `io[b"_"][0]` (column view → `AsyncSingleColumnView`) which returns the
raw bytes value directly via `_read_column`. No dict unpacking needed.

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
| `result_backends.py` | Use `storage[key][b"_"][0]` column view for get; `storage[key].extend/clear` for store/delete |
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
4. **All call sites use the pandas-like subscript API** (`io[...]`, `.to_list()`,
   `.to_dict()`) as the primary interface. `.get(keys=...)` only where a filtered
   dict is explicitly required.
5. Isosurface uses `io[cube_key][index]` — column view, single value.
6. Result backend uses `io[b"_"][0]` — column view, single value.
7. Batch reads use `io[start:stop].to_list()` / `io[indices].to_list()`.
8. Single `FrameStorageDep` replaces both `StorageDep` and `StorageRouterDep`.
9. Write guard enforced via `RequireWritableDep` FastAPI dependency.

## Non-Goals

- Changing the `AsyncBlobIO` or asebytes library.
- Changing the REST API contract (all endpoints keep the same request/response format).
- Optimizing batch+keys reads (kept as post-load filter; minor impact).
- Changing the provider dispatch mechanism.
