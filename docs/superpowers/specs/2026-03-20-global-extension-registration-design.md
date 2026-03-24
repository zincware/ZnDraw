# Global Extension Registration

## Summary

Enable global (`@global`) and per-room extension registration through the public
`ZnDraw` client API with proper admin enforcement, and provide backwards-compatible
deprecation shims for the v0.6.0 `register_extension()` method.

## Motivation

Users migrating from v0.6.0 to v0.7.0 hit `AttributeError` on
`vis.register_extension(cls, public=True)` because the method was removed.
The replacement `register_job()` only supports room-scoped registration — there
is no public API path to register a global extension visible in all rooms.

The server already supports `@global` job registration (admin-only, enforced at
`PUT /v1/joblib/rooms/@global/jobs` with 403 for non-superusers), but the client
never exposes it.

## Design

### Constant: `GLOBAL_ROOM`

```python
# src/zndraw/__init__.py
GLOBAL_ROOM = "@global"
```

Exported from the package so users don't need to remember the string. Rooms
cannot start with `@`, so there is no collision risk.

### `register_job` signature change

```python
def register_job(
    self,
    cls: type,
    *,
    room: Literal["@global"] | str | None = None,
    public: Annotated[
        bool | None,
        deprecated("Use room='@global' instead of public=True"),
    ] = None,
) -> None:
```

Logic:

1. If `public is True` and `room` is also set, raise `ValueError`.
2. If `public is True`, emit `DeprecationWarning` via `warnings.warn` and set
   `room = GLOBAL_ROOM`.
3. If `room == GLOBAL_ROOM`, pass it directly to `self.jobs.register()`.
4. Otherwise, use `self._resolve_room(room)` (falls back to `self.room`).

Note: `public=False` (explicit) is treated the same as `public=None` (default) —
only `public=True` triggers the conflict check and deprecation path.

```python
def register_job(self, cls, *, room=None, public=None):
    if public and room is not None:
        raise ValueError("Cannot specify both 'room' and 'public'")
    if public:
        warnings.warn(
            "public=True is deprecated, use room='@global' instead",
            DeprecationWarning,
            stacklevel=2,
        )
        room = GLOBAL_ROOM
    elif room != GLOBAL_ROOM:
        room = self._resolve_room(room)

    self._ensure_socket_connected()
    self.jobs.register(cls, room=room)
```

### Deprecated `register_extension` method

```python
@deprecated("Use register_job(cls, room='@global') for global, or register_job(cls) for room-scoped")
def register_extension(self, cls: type, *, public: bool = False, **kwargs) -> None:
    room = GLOBAL_ROOM if public else kwargs.get("room")
    self.register_job(cls, room=room)
```

### No server changes

Admin enforcement for `@global` registration already exists in
`zndraw_joblib/router.py:269`. The client surfaces the HTTP 403 as-is.

## Tests

1. **Admin registers `@global` extension** — succeeds, extension is registered.
2. **Non-admin tries `@global` registration** — gets 403 from server.
3. **Global extension visible in all rooms** — register `@global` extension from
   one room, verify it appears in `vis.extensions` from a different room.
4. **Deprecated `register_extension(cls, public=True)`** — works, emits
   deprecation warning, delegates correctly.
5. **`room` and `public` both set** — raises `ValueError`.

## Files to change

- `src/zndraw/__init__.py` — export `GLOBAL_ROOM`
- `src/zndraw/client/core.py` — modify `register_job`, add `register_extension`
- `tests/` — new test file or extend existing integration tests

## Out of scope

- Server-side changes (already enforced).
- Client-side pre-check of admin status (simple try-and-fail approach).
- Migration guide / changelog (separate task).
