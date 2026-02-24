# Progress Bar Design

## Summary

Migrate the progress bar feature from the old ZnDraw implementation.
REST endpoints for mutations, Socket.IO broadcasts for real-time updates,
tqdm subclass on the Python client for ergonomic usage with built-in debouncing.

## Decision Record

- **Transport:** REST + socket broadcast (matches chat, frames pattern)
- **Storage:** Redis hash (ephemeral, no SQL)
- **Client API:** `ZnDrawTqdm` tqdm subclass with `mininterval=0.5` throttling
- **Initial sync:** Progress trackers included in `RoomJoinResponse`
- **Chat integration:** None — progress is ephemeral, chat is persistent

## Data Model

### Socket Broadcast Events (`socket_events.py`)

No `room_id` — events are broadcast to the room, so listeners know which room.

```python
class ProgressStart(BaseModel):
    """Broadcast when a new tracker is created."""
    progress_id: str
    description: str

class ProgressUpdate(BaseModel):
    """Broadcast when a tracker is updated."""
    progress_id: str
    description: str | None = None
    progress: float | None = None

class ProgressComplete(BaseModel):
    """Broadcast when a tracker finishes."""
    progress_id: str
```

### REST Schemas (`schemas.py`)

```python
class ProgressCreate(BaseModel):
    progress_id: str
    description: str

class ProgressPatch(BaseModel):
    description: str | None = None
    progress: float | None = None

class ProgressResponse(BaseModel):
    progress_id: str
    description: str
    progress: float | None = None
```

### Redis Storage

Hash key: `room:{room_id}:progress`
Field: `{progress_id}`, Value: JSON `ProgressResponse`.

```python
# RedisKey
@staticmethod
def room_progress(room_id: str) -> str:
    return f"room:{room_id}:progress"
```

### RoomJoinResponse Extension

```python
class RoomJoinResponse(BaseModel):
    ...  # existing fields
    progress_trackers: dict[str, ProgressResponse] = {}
```

## REST Endpoints

File: `src/zndraw/routes/progress.py`

| Method | Path | Body | Action |
|--------|------|------|--------|
| POST   | `/v1/rooms/{room_id}/progress` | `ProgressCreate` | Store in Redis, broadcast `ProgressStart` |
| PATCH  | `/v1/rooms/{room_id}/progress/{progress_id}` | `ProgressPatch` | Update Redis, broadcast `ProgressUpdate` |
| DELETE | `/v1/rooms/{room_id}/progress/{progress_id}` | — | Remove from Redis, broadcast `ProgressComplete` |

Response codes: 201, 200, 204.

## Python Client: `ZnDrawTqdm`

tqdm subclass with built-in debouncing via `mininterval` (default 0.5s).
Suppresses terminal output. Fires REST calls only on tqdm refresh.

```python
class ZnDrawTqdm(std_tqdm):
    def __init__(self, *args, vis, description="Processing...", mininterval=0.5, **kwargs):
        kwargs["mininterval"] = mininterval
        super().__init__(*args, **kwargs)
        self._vis_api = vis.api
        self._progress_id = str(uuid4())
        self._vis_api.progress_start(self._progress_id, description)

    def display(self, *args, **kwargs):
        pass

    def update(self, n=1):
        displayed = super().update(n)
        if displayed:
            d = self.format_dict
            progress = (d["n"] / d["total"] * 100) if d["total"] else None
            self._vis_api.progress_update(self._progress_id, progress=progress)
        return displayed

    def close(self):
        self._vis_api.progress_complete(self._progress_id)
        super().close()
```

**Throttling:** With `mininterval=0.5`, max 2 REST calls/second regardless
of iteration speed. `list(ZnDrawTqdm(range(50000)))` sends ~1 update total.

## Frontend Changes

Minimal — UI and state management already exist.

1. **`RoomJoinResponse`:** Extract `progress_trackers` and call `setProgressTrackers()`
2. **Remove `room_id` filtering** from `ProgressNotifications.tsx` — events are room-scoped
3. **Remove `progress_init` handler** — data comes via `RoomJoinResponse`

## Testing

- Unit: Progress REST endpoints (create, update, delete, 404 on missing)
- Integration: REST create → verify socket broadcast to room
- E2E: `ZnDrawTqdm` → progress updates reach room clients
