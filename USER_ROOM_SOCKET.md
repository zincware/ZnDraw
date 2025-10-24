# Refactoring user / room / socket management

## Overview

Consolidate socket management to use a single connection per client with room-based broadcasting and unified room metadata events.

## Architecture Goals

1. **Single Socket Connection** - One WebSocket connection per client
2. **Room-Based Broadcasting** - Use Flask-SocketIO's `join_room`/`leave_room`
3. **Consolidated Events** - Replace fragmented events with unified `room:update` and `room:delete`
4. **Typed Data Model** - Pydantic model shared between REST API and WebSocket
5. **Efficient Updates** - Send only changed fields, not full objects
6. **No Backward Compatibility** - Clean removal of legacy code

## Socket Connection Flow

### Client Connection Lifecycle

1. **On Page Load** - Client connects to WebSocket (ONE connection for entire session)
2. **On `/rooms` Page** - Client joins `overview:public` room
3. **On `/rooms/:roomId/:userId` Page** - Client:
   - Leaves `overview:public` room (if joined)
   - Joins `room:<room_id>` room
4. **On Navigation Back to `/rooms`** - Client:
   - Leaves `room:<room_id>` room
   - Re-joins `overview:public` room

### Initial Data Loading

- **On `/rooms` page load** - ONE-TIME REST API call to `GET /api/rooms` (full list)
- **On specific room page load** - ONE-TIME REST API call to `GET /api/rooms/:roomId` (room detail)
- **Subsequent updates** - Via WebSocket events only

## Room Metadata Model

### Pydantic Model Definition

```python
from pydantic import BaseModel, Field
from typing import Optional

class LockMetadata(BaseModel):
    """Lock metadata for temporary locks (e.g., uploads)."""
    msg: Optional[str] = None
    userName: Optional[str] = None
    timestamp: Optional[float] = None  # Unix timestamp

class RoomMetadata(BaseModel):
    """
    Room metadata shared between REST API and WebSocket events.
    This is the single source of truth for room state.
    """
    id: str = Field(..., description="Unique room identifier")
    description: Optional[str] = Field(None, description="Human-readable description")
    frameCount: int = Field(0, description="Number of trajectory frames")
    locked: bool = Field(False, description="Permanent lock (immutable)")
    hidden: bool = Field(False, description="Hidden from room list")
    isDefault: bool = Field(False, description="Set as default room")
    metadataLocked: Optional[LockMetadata] = Field(None, description="Temporary lock metadata")
    presenterSid: Optional[str] = Field(None, description="Socket ID of current presenter")

    # Additional metadata (file paths, etc.) stored separately via RoomMetadataManager
    # Not included in socket events
```

### Fields Explanation

| Field | Type | Description | Redis Key |
|-------|------|-------------|-----------|
| `id` | `str` | Room identifier | Key prefix `room:{id}:*` |
| `description` | `str \| None` | Room description | `room:{id}:description` |
| `frameCount` | `int` | Number of frames | Computed from `room:{id}:trajectory:indices` |
| `locked` | `bool` | Permanent lock | `room:{id}:locked` |
| `hidden` | `bool` | Hidden from list | `room:{id}:hidden` |
| `isDefault` | `bool` | Is default room | `default_room` global key |
| `metadataLocked` | `LockMetadata \| None` | Temporary lock info | `lock:trajectory:meta:{room_id}` (ttl-based) |
| `presenterSid` | `str \| None` | Current presenter socket ID | `presenter:{room_id}` (in-memory or Redis) |

## Event Consolidation

### Events to REMOVE (Legacy)

| Old Event | Replaced By | Reason |
|-----------|-------------|--------|
| `len_frames` | `room:update` | Consolidate frame count updates |
| `lock:metadata:updated` | `room:update` | Consolidate lock state broadcast |
| `lock:released` | `room:update` | Consolidate lock state broadcast |
| `presenter_update` | `room:update` | Consolidate presenter state |
| Current `room:updated` | `room:update` (standardized) | Unify naming and structure |

### Events to KEEP (Room-Specific Operations)

These events are for specific room operations, NOT room metadata:

- `frame_update` - Current frame index changed (viewer state)
- `invalidate:*` - Cache invalidation (settings, geometry, figure, selection, etc.)
- `frames:invalidate` - Frame cache invalidation
- `bookmarks:invalidate` - Bookmarks changed
- `queue:update` - Job queue updates
- `frame_selection:update` - Frame selection changed
- `chat:message:new` - New chat message
- `chat:message:updated` - Updated chat message
- **Lock operation handlers** - `lock:acquire`, `lock:release`, `lock:refresh`, `lock:msg` (use `socket.call` for validation, then emit `room:update`)

### New Consolidated Events

#### 1. `room:update` Event

**Purpose:** All room metadata changes (frames, locks, description, etc.)

**Payload Structure:**
```python
{
    "roomId": str,           # Required: which room changed
    "created": bool,         # Optional: true if this is a new room
    # ... only changed fields below ...
    "description": str | None,
    "frameCount": int,
    "locked": bool,
    "hidden": bool,
    "isDefault": bool,
    "metadataLocked": LockMetadata | None,
    "presenterSid": str | None
}
```

**Examples:**
```python
# Frame count updated
{"roomId": "my-room", "frameCount": 120}

# Lock acquired
{"roomId": "my-room", "metadataLocked": {"msg": "Uploading", "userName": "alice", "timestamp": 1234567890.0}}

# Lock released
{"roomId": "my-room", "metadataLocked": None}

# Presenter changed
{"roomId": "my-room", "presenterSid": "socket-123-abc"}

# Presenter cleared
{"roomId": "my-room", "presenterSid": None}

# Room created
{"roomId": "new-room", "created": True, "description": "New Room", "frameCount": 0, "locked": False, "hidden": False, "isDefault": False}

# Multiple fields updated
{"roomId": "my-room", "description": "Updated desc", "locked": True}
```

**Broadcasting:**
- **To `overview:public`** - Always broadcast (for room list updates)
- **To `room:<room_id>`** - Always broadcast (for room viewers)

#### 2. `room:delete` Event

**Purpose:** Room deletion

**Payload Structure:**
```python
{
    "roomId": str  # Which room was deleted
}
```

**Broadcasting:**
- **To `overview:public`** - Always broadcast
- **To `room:<room_id>`** - Always broadcast (clients should navigate away)

## Implementation Plan

### Phase 1: Backend - Data Model (1-2 hours)

**File:** `src/zndraw/app/models.py` (new file)

```python
from pydantic import BaseModel, Field
from typing import Optional

class LockMetadata(BaseModel):
    msg: Optional[str] = None
    userName: Optional[str] = None
    timestamp: Optional[float] = None

class RoomMetadata(BaseModel):
    id: str
    description: Optional[str] = None
    frameCount: int = 0
    locked: bool = False
    hidden: bool = False
    isDefault: bool = False
    metadataLocked: Optional[LockMetadata] = None
    presenterSid: Optional[str] = None
```

**File:** `src/zndraw/app/room_manager.py` (new file)

Create helper functions:
```python
def get_room_metadata(redis_client, room_id: str) -> RoomMetadata:
    """Fetch complete room metadata from Redis."""
    pass

def emit_room_update(socketio, room_id: str, **changes):
    """
    Emit room:update event to both overview:public and room:<room_id>.
    Only sends changed fields.
    """
    payload = {"roomId": room_id, **changes}
    socketio.emit("room:update", payload, to="overview:public", namespace="/")
    socketio.emit("room:update", payload, to=f"room:{room_id}", namespace="/")
```

### Phase 2: Backend - Socket Event Handlers (2-3 hours)

**File:** `src/zndraw/app/events.py`

**Add new handlers:**
```python
@socketio.on("join:overview")
def on_join_overview():
    """Client joining /rooms page."""
    join_room("overview:public")
    emit({"status": "joined", "room": "overview:public"})

@socketio.on("leave:overview")
def on_leave_overview():
    """Client leaving /rooms page."""
    leave_room("overview:public")
    emit({"status": "left", "room": "overview:public"})

@socketio.on("join:room")
def on_join_room(data):
    """Client joining specific room page."""
    room_id = data.get("roomId")

    # Leave overview if joined
    leave_room("overview:public")

    # Join specific room
    join_room(f"room:{room_id}")
    emit({"status": "joined", "room": f"room:{room_id}"})

@socketio.on("leave:room")
def on_leave_room(data):
    """Client leaving specific room page."""
    room_id = data.get("roomId")
    leave_room(f"room:{room_id}")
    emit({"status": "left", "room": f"room:{room_id}"})
```

**Update existing handlers:**

1. **Lock operations** - KEEP `lock:acquire`, `lock:release`, `lock:refresh`, `lock:msg` handlers (they use `socket.call` for validation)
   - After successful validation, emit `room:update` with lock state:
   ```python
   # In lock:msg handler (around line 560):
   # After validation, broadcast state change:
   emit_room_update(
       socketio,
       room=room,
       metadataLocked=metadata_with_user
   )

   # In lock:release handler (around line 475):
   # After validation, broadcast state change:
   emit_room_update(
       socketio,
       room=room,
       metadataLocked=None
   )
   ```

2. **Presenter mode** (around line 170):
```python
# BEFORE: emit("presenter_update", {"presenterSid": ...}, to=f"room:{room_name}")
# AFTER:
emit_room_update(
    socketio,
    room_id=room_name,
    presenterSid=presenter_sid  # or None to clear
)
```

### Phase 3: Backend - REST Endpoints (2-3 hours)

**File:** `src/zndraw/app/routes.py`

**Update all `socketio.emit("room:updated", ...)` calls:**

1. **Frame count updates** (around line 245-251):
```python
# REMOVE: socketio.emit("len_frames", ...)
# CHANGE: socketio.emit("room:updated", ...) → emit_room_update(...)
emit_room_update(socketio, room_id, frameCount=frame_count)
```

2. **Room updates** (around line 1239):
```python
# CHANGE: socketio.emit("room:updated", event_data, namespace="/")
# TO:
emit_room_update(socketio, room_id, **changes)
```

3. **Default room updates** (around line 1284-1316):
```python
# CHANGE: socketio.emit("room:updated", {"roomId": ..., "isDefault": ...})
# TO:
emit_room_update(socketio, room_id, isDefault=True)
emit_room_update(socketio, previous_default, isDefault=False)
```

4. **Room creation** (around line 2542):
```python
# CHANGE: socketio.emit("room:updated", {..., "created": True})
# TO:
emit_room_update(
    socketio,
    room_id,
    created=True,
    description=description,
    frameCount=frame_count,
    locked=False,
    hidden=False,
    isDefault=False
)
```

### Phase 4: Frontend - Single Socket Manager (3-4 hours)

**DELETE:** `app/src/hooks/useGlobalSocketListener.ts` (entire file)

**File:** `app/src/hooks/useSocketManager.ts`

**Add room management:**
```typescript
interface SocketManagerOptions {
  roomId?: string;  // If on /rooms/:roomId page
  isOverview?: boolean;  // If on /rooms page
}

export const useSocketManager = (options: SocketManagerOptions = {}) => {
  const { roomId, isOverview } = options;

  useEffect(() => {
    // ... existing socket creation ...

    // Handle room joining/leaving
    if (isOverview) {
      socket.emit("join:overview");
    } else if (roomId) {
      socket.emit("join:room", { roomId });
    }

    return () => {
      if (isOverview) {
        socket.emit("leave:overview");
      } else if (roomId) {
        socket.emit("leave:room", { roomId });
      }
    };
  }, [roomId, isOverview]);

  // ... rest of socket manager ...
}
```

**Update event handlers:**
```typescript
// REMOVE: socket.on("len_frames", ...)
// REMOVE: socket.on("lock:metadata:updated", ...)
// REMOVE: socket.on("lock:released", ...)

// CHANGE: socket.on("room:updated", ...) → socket.on("room:update", ...)
function onRoomUpdate(data: any) {
  const { roomId, created, ...updates } = data;

  if (created) {
    // New room created
    useRoomsStore.getState().setRoom(roomId, {
      id: roomId,
      frameCount: 0,
      locked: false,
      hidden: false,
      isDefault: false,
      ...updates
    });
  } else {
    // Update existing room
    useRoomsStore.getState().updateRoom(roomId, updates);
  }
}

socket.on("room:update", onRoomUpdate);

// NEW: Handle room deletion
function onRoomDelete(data: any) {
  useRoomsStore.getState().removeRoom(data.roomId);
}

socket.on("room:delete", onRoomDelete);
```

### Phase 5: Frontend - Component Updates (1-2 hours)

**File:** `app/src/pages/roomList.tsx`

```typescript
export default function RoomListPage() {
  // REMOVE: useGlobalSocketListener();

  // ADD: Use main socket manager with overview flag
  useSocketManager({ isOverview: true });

  // ... rest remains the same ...
}
```

**File:** `app/src/pages/landingPage.tsx` (main room viewer)

```typescript
export default function LandingPage() {
  const { roomId, userId } = useParams();

  // UPDATE: Pass roomId to socket manager
  useSocketManager({ roomId });

  // ... rest remains the same ...
}
```

**File:** `app/src/store.tsx`

Update lock metadata type:
```typescript
interface LockMetadata {
  msg?: string;
  userName?: string;
  timestamp?: number;
}

interface AppState {
  lockMetadata: LockMetadata | null;
  // ...
}
```

### Phase 6: Testing & Cleanup (1-2 hours)

1. **Test room list updates** - Upload file, verify frame count updates in real-time
2. **Test lock metadata** - Verify lock pills show/hide correctly
3. **Test room creation** - New rooms appear without refresh
4. **Test navigation** - Verify join/leave room events fire correctly
5. **Test default room** - Setting default updates all clients
6. **Remove legacy code** - Search for and delete:
   - All `len_frames` references
   - All `lock:metadata:updated` references
   - All `lock:released` references
   - `useGlobalSocketListener.ts` file

## Migration Checklist

### Backend Changes

- [ ] Create `src/zndraw/app/models.py` with Pydantic models
- [ ] Create `src/zndraw/app/room_manager.py` with helper functions
- [ ] Add join/leave room handlers in `events.py`
- [ ] Update lock metadata emit in `events.py` (line 587)
- [ ] Update lock release emit in `events.py` (line 475)
- [ ] Update frame count emits in `routes.py` (line 245, 251)
- [ ] Update room update emit in `routes.py` (line 1239)
- [ ] Update default room emits in `routes.py` (line 1284-1316)
- [ ] Update room creation emit in `routes.py` (line 2542)
- [ ] Remove all legacy event emits

### Frontend Changes

- [ ] Delete `app/src/hooks/useGlobalSocketListener.ts`
- [ ] Update `useSocketManager.ts` with room join/leave logic
- [ ] Remove `len_frames` event handler
- [ ] Remove `lock:metadata:updated` event handler
- [ ] Remove `lock:released` event handler
- [ ] Rename `room:updated` → `room:update` handler
- [ ] Add `room:delete` handler
- [ ] Update `roomList.tsx` to use `useSocketManager({ isOverview: true })`
- [ ] Update `landingPage.tsx` to use `useSocketManager({ roomId })`
- [ ] Update lock metadata TypeScript types in `store.tsx`
- [ ] Update Room interface in `myapi/client.ts` (align with Pydantic model)

### Testing

- [ ] Test `/rooms` page receives updates (frame count, locks, new rooms)
- [ ] Test room viewer receives updates
- [ ] Test navigation between `/rooms` and specific room
- [ ] Test lock metadata display (yellow pills)
- [ ] Test permanent lock display (red pills)
- [ ] Test default room setting
- [ ] Test room creation
- [ ] Verify no duplicate socket connections in browser DevTools

## Notes

- **No debouncing** - Not implementing for now, can add later if needed
- **Same data structure** - `room:update` event uses identical structure for overview and room broadcasts
- **Msgpack** - Can be added later for optimization, start with JSON
- **No backward compatibility** - Clean break, remove all legacy code
- **Lock operations flow** - Client sends `lock:*` via `socket.call` → Server validates → Server responds to client → Server emits `room:update` with lock state to all other clients
- **Presenter mode** - State changes emit `room:update` with `presenterSid` field
