# /join Route Refactoring Plan

## Current State

**File:** `src/zndraw/app/routes.py:2259-2495`

**Problems:**
- 236 lines doing too much (violates SRP)
- Returns MASSIVE payload on every join
- 15+ separate Redis operations in single request
- No lazy loading - fetches ALL data upfront
- Single point of failure
- Hard to test and maintain
- Response size grows with room complexity

**Current Responsibilities:**
1. Validates room ID
2. Generates/validates client IDs
3. Creates join tokens (60s expiry)
4. Creates new rooms OR copies from existing
5. Initializes default geometries (Sphere, Bond, Curve, Cell, Floor)
6. Fetches ALL geometries
7. Fetches ALL bookmarks
8. Fetches ALL settings
9. Fetches ALL selection groups
10. Returns current frame, frame count, selections

**Data Returned:**
```json
{
  "clientId", "joinToken", "frameCount", "step",
  "selections", "selectionGroups", "activeSelectionGroup",
  "frame_selection", "geometries", "geometryDefaults",
  "bookmarks", "settings", "presenter-lock", "created"
}
```

## Target Architecture

### Phase 1: Extract Service Layer

**New Structure:**
```
src/zndraw/
├── app/
│   └── routes.py (thin controller - 20 lines)
├── services/
│   ├── room_service.py (NEW - room lifecycle)
│   └── client_service.py (NEW - client management)
```

**Service Responsibilities:**
- `RoomService`: Create, validate, initialize rooms
- `ClientService`: Manage client sessions, tokens, metadata

### Phase 2: Minimize Join Response

**Principle:** Return ONLY what's needed to establish connection

**New Response (90% smaller):**
```json
{
  "status": "ok",
  "clientId": "uuid",
  "joinToken": "token",
  "frameCount": 0,
  "step": 0,
  "created": true,
  "roomExists": true
}
```

**Everything else fetched lazily:**
- Geometries → fetched when canvas mounts
- Bookmarks → fetched when bookmarks panel opens
- Settings → fetched when settings accessed
- Selection groups → fetched when selections UI opens

## Implementation Steps

### Step 1: Create Service Layer (2-3 hours)

**File:** `src/zndraw/services/room_service.py`

```python
"""Room lifecycle and data management."""
from typing import Optional
import json
import uuid
import logging
from redis import Redis

log = logging.getLogger(__name__)


class RoomService:
    """Handles room creation, validation, and data management.

    Follows SOLID principles:
    - Single Responsibility: Only manages room lifecycle
    - Open/Closed: Extensible through inheritance
    - Interface Segregation: Focused interface
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def room_exists(self, room_id: str) -> bool:
        """Check if room exists by checking for current_frame key."""
        return self.r.exists(f"room:{room_id}:current_frame")

    def create_room(
        self,
        room_id: str,
        user_name: str,
        description: Optional[str] = None,
        copy_from: Optional[str] = None
    ) -> dict:
        """Create new room, optionally copying from existing room.

        Returns:
            dict with keys: created (bool), frameCount (int)
        """
        if copy_from:
            return self._create_room_from_copy(room_id, copy_from)
        return self._create_empty_room(room_id, user_name, description)

    def _create_empty_room(
        self,
        room_id: str,
        user_name: str,
        description: Optional[str] = None
    ) -> dict:
        """Create empty room with defaults."""
        # Set description
        if description:
            self.r.set(f"room:{room_id}:description", description)

        # Initialize metadata
        self.r.set(f"room:{room_id}:current_frame", 0)
        self.r.set(f"room:{room_id}:locked", 0)
        self.r.set(f"room:{room_id}:hidden", 0)

        # Initialize default settings
        self._initialize_default_settings(room_id, user_name)

        # Create default geometries
        self._initialize_default_geometries(room_id)

        log.info(f"Created empty room '{room_id}'")
        return {"created": True, "frameCount": 0}

    def _create_room_from_copy(self, room_id: str, source_room: str) -> dict:
        """Copy room data from existing room."""
        source_indices_key = f"room:{source_room}:trajectory:indices"
        if not self.r.exists(source_indices_key):
            raise ValueError(f"Source room '{source_room}' not found")

        # Copy trajectory indices (shares frame data)
        source_indices = self.r.zrange(source_indices_key, 0, -1, withscores=True)
        if source_indices:
            self.r.zadd(
                f"room:{room_id}:trajectory:indices",
                {member: score for member, score in source_indices}
            )

        # Copy geometries
        geometries = self.r.hgetall(f"room:{source_room}:geometries")
        if geometries:
            self.r.hset(f"room:{room_id}:geometries", mapping=geometries)

        # Copy bookmarks
        bookmarks = self.r.hgetall(f"room:{source_room}:bookmarks")
        if bookmarks:
            self.r.hset(f"room:{room_id}:bookmarks", mapping=bookmarks)

        # Initialize metadata
        self.r.set(f"room:{room_id}:current_frame", 0)
        self.r.set(f"room:{room_id}:locked", 0)
        self.r.set(f"room:{room_id}:hidden", 0)

        log.info(f"Created room '{room_id}' from '{source_room}' with {len(source_indices)} frames")
        return {"created": True, "frameCount": len(source_indices)}

    def _initialize_default_settings(self, room_id: str, user_name: str):
        """Initialize default RoomConfig settings for new room."""
        from zndraw.settings import RoomConfig

        default_config = RoomConfig()
        config_dict = default_config.model_dump()

        for category_name, category_data in config_dict.items():
            self.r.hset(
                f"room:{room_id}:user:{user_name}:settings",
                category_name,
                json.dumps(category_data)
            )

    def _initialize_default_geometries(self, room_id: str):
        """Initialize default geometries for new room."""
        from zndraw.geometries import Sphere, Bond, Curve, Cell, Floor

        defaults = {
            "particles": (Sphere, {"position": "arrays.positions", "color": "arrays.colors"}),
            "bonds": (Bond, {}),
            "curve": (Curve, {}),
            "cell": (Cell, {}),
            "floor": (Floor, {})
        }

        for key, (geometry_class, kwargs) in defaults.items():
            geometry_data = geometry_class(**kwargs).model_dump()
            self.r.hset(
                f"room:{room_id}:geometries",
                key,
                json.dumps({"type": geometry_class.__name__, "data": geometry_data})
            )

    def get_frame_count(self, room_id: str) -> int:
        """Get number of frames in room."""
        return self.r.zcard(f"room:{room_id}:trajectory:indices")

    def get_current_frame(self, room_id: str) -> int:
        """Get current frame number, handling invalid values."""
        step = self.r.get(f"room:{room_id}:current_frame")
        try:
            if step is not None:
                step_int = int(step)
                if step_int < 0:
                    log.warning(f"Negative frame in Redis for room {room_id}: {step_int}, resetting to 0")
                    return 0
                return step_int
            return 0
        except (ValueError, TypeError) as e:
            log.error(f"Invalid step value in Redis for room {room_id}: {step} - {e}")
            return 0
```

**File:** `src/zndraw/services/client_service.py`

```python
"""Client session and authentication management."""
import json
import uuid
import logging
from redis import Redis

log = logging.getLogger(__name__)


class ClientService:
    """Handles client sessions, tokens, and metadata."""

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def get_or_create_client(
        self,
        room_id: str,
        user_name: str,
        client_id: Optional[str] = None
    ) -> str:
        """Get existing or create new client ID."""
        if not client_id:
            client_id = str(uuid.uuid4())
            log.info(f"Generated new clientId: {client_id}")

        # Update client metadata
        client_key = f"client:{client_id}"
        self.r.hset(client_key, "userName", user_name)
        self.r.hset(client_key, "currentRoom", room_id)

        # Add to room's client set
        self.r.sadd(f"room:{room_id}:clients", client_id)

        return client_id

    def create_join_token(
        self,
        client_id: str,
        room_id: str,
        user_name: str,
        expiry_seconds: int = 60
    ) -> str:
        """Create temporary join token for socket authentication."""
        join_token = str(uuid.uuid4())
        token_key = f"join_token:{join_token}"

        token_data = json.dumps({
            "clientId": client_id,
            "roomId": room_id,
            "userName": user_name
        })

        self.r.setex(token_key, expiry_seconds, token_data)

        log.info(f"Client {client_id} ({user_name}) join token created: {join_token}")
        return join_token
```

### Step 2: Refactor Route (30 minutes)

**File:** `src/zndraw/app/routes.py` (replace lines 2259-2495)

```python
@main.route("/api/rooms/<string:room_id>/join", methods=["POST"])
def join_room(room_id: str):
    """Join a room - returns minimal data for connection establishment.

    All other data (geometries, bookmarks, settings) should be fetched lazily.
    """
    from zndraw.services.room_service import RoomService
    from zndraw.services.client_service import ClientService

    data = request.get_json() or {}

    # Validate room ID
    if ":" in room_id:
        return {"error": "Room ID cannot contain ':' character"}, 400

    # Extract parameters
    user_name = data.get("userId", "Anonymous")
    client_id = data.get("clientId")
    description = data.get("description")
    copy_from = data.get("copyFrom")
    allow_create = data.get("allowCreate", True)

    # Initialize services
    r = current_app.extensions["redis"]
    room_service = RoomService(r)
    client_service = ClientService(r)

    # Check if room exists
    room_exists = room_service.room_exists(room_id)

    # Handle room creation
    if not room_exists:
        if not allow_create:
            return {
                "status": "not_found",
                "message": f"Room '{room_id}' does not exist yet. It may still be loading."
            }, 404

        try:
            room_service.create_room(room_id, user_name, description, copy_from)
        except ValueError as e:
            return {"error": str(e)}, 404

    # Get or create client
    client_id = client_service.get_or_create_client(room_id, user_name, client_id)

    # Create join token
    join_token = client_service.create_join_token(client_id, room_id, user_name)

    # Return minimal response
    return {
        "status": "ok",
        "clientId": client_id,
        "joinToken": join_token,
        "frameCount": room_service.get_frame_count(room_id),
        "step": room_service.get_current_frame(room_id),
        "created": not room_exists,
        "roomId": room_id
    }
```

### Step 3: Create Lazy Loading Endpoints (2 hours)

**Add new endpoints for lazy data fetching:**

```python
@main.route("/api/rooms/<string:room_id>/geometries", methods=["GET"])
def get_geometries(room_id: str):
    """Fetch all geometries for a room."""
    r = current_app.extensions["redis"]
    geometries = r.hgetall(f"room:{room_id}:geometries")

    return {
        "geometries": {k: json.loads(v) for k, v in geometries.items()},
        "geometryDefaults": _get_geometry_defaults()
    }

@main.route("/api/rooms/<string:room_id>/bookmarks", methods=["GET"])
def get_bookmarks(room_id: str):
    """Fetch bookmarks for a room."""
    r = current_app.extensions["redis"]
    bookmarks_raw = r.hgetall(f"room:{room_id}:bookmarks")

    if bookmarks_raw:
        return {"bookmarks": {int(k): v for k, v in bookmarks_raw.items()}}
    return {"bookmarks": None}

@main.route("/api/rooms/<string:room_id>/selections", methods=["GET"])
def get_selections(room_id: str):
    """Fetch selections and selection groups."""
    r = current_app.extensions["redis"]

    # Per-geometry selections
    selections_raw = r.hgetall(f"room:{room_id}:selections")
    selections = {k: json.loads(v) for k, v in selections_raw.items()}

    # Selection groups
    groups_raw = r.hgetall(f"room:{room_id}:selection_groups")
    selection_groups = {k: json.loads(v) for k, v in groups_raw.items()}

    # Active group
    active_group = r.get(f"room:{room_id}:active_selection_group")

    # Frame selection
    frame_selection = r.get(f"room:{room_id}:frame_selection:default")

    return {
        "selections": selections,
        "selectionGroups": selection_groups,
        "activeSelectionGroup": active_group,
        "frame_selection": json.loads(frame_selection) if frame_selection else None
    }

@main.route("/api/rooms/<string:room_id>/user/<string:user_name>/settings", methods=["GET"])
def get_user_settings(room_id: str, user_name: str):
    """Fetch user settings for a room."""
    r = current_app.extensions["redis"]
    settings_data = {}
    settings_hash = r.hgetall(f"room:{room_id}:user:{user_name}:settings")

    for setting_name, setting_value in settings_hash.items():
        settings_data[setting_name] = json.loads(setting_value)

    return {"settings": settings_data}
```

### Step 4: Update Frontend (2-3 hours)

**File:** `app/src/hooks/useRestManager.ts`

**Remove TODO on line 59 and implement:**

```typescript
export const useRestJoinManager = () => {
  const {
    setClientId,
    setRoomId,
    setUserId,
    setCurrentFrame,
    setFrameCount,
    setJoinToken,
  } = useAppStore();

  const { roomId: room, userId } = useParams<{
    roomId: string;
    userId: string;
  }>();

  const joinRoom = useCallback(async () => {
    if (!room || !userId) return;

    console.log("Joining room via REST:", room, userId);

    try {
      // Minimal join - only essential data
      const data = await joinRoomApi(room, { userId });
      console.log("Join response data:", data);

      // Update only essential state
      if (data.clientId) setClientId(data.clientId);
      if (data.joinToken) setJoinToken(data.joinToken);
      if (typeof data.frameCount === "number") setFrameCount(data.frameCount);
      if (data.step !== undefined && data.step !== null) setCurrentFrame(data.step);

      setRoomId(room);
      setUserId(userId);

      // Everything else loaded lazily by other hooks/components

    } catch (error) {
      console.error("Error joining room:", error);
    }
  }, [room, userId, setClientId, setRoomId, setUserId, setCurrentFrame, setFrameCount, setJoinToken]);

  // ... rest of hook
};
```

**Add new hooks for lazy loading:**

```typescript
// app/src/hooks/useLazyGeometries.ts
export const useLazyGeometries = () => {
  const { roomId, setGeometries, setGeometryDefaults } = useAppStore();

  const { data, isLoading } = useQuery({
    queryKey: ['geometries', roomId],
    queryFn: () => getGeometries(roomId),
    enabled: !!roomId,
    staleTime: Infinity, // Cache forever unless invalidated
  });

  useEffect(() => {
    if (data) {
      setGeometryDefaults(data.geometryDefaults);
      setGeometries(data.geometries);
    }
  }, [data]);

  return { isLoading };
};

// app/src/hooks/useLazySelections.ts
export const useLazySelections = () => {
  const { roomId, setSelections, setSelectionGroups, setActiveSelectionGroup, setFrameSelection } = useAppStore();

  const { data } = useQuery({
    queryKey: ['selections', roomId],
    queryFn: () => getSelections(roomId),
    enabled: !!roomId,
    staleTime: 60000, // Cache for 1 minute
  });

  useEffect(() => {
    if (data) {
      if (data.selections) setSelections(data.selections);
      if (data.selectionGroups) setSelectionGroups(data.selectionGroups);
      if (data.activeSelectionGroup !== undefined) setActiveSelectionGroup(data.activeSelectionGroup);
      if (data.frame_selection !== undefined) setFrameSelection(data.frame_selection);
    }
  }, [data]);
};
```

**Use in components:**

```typescript
// app/src/pages/room.tsx
export function RoomPage() {
  useRestJoinManager(); // Minimal join
  const { isLoading: geometriesLoading } = useLazyGeometries(); // Lazy load geometries
  useLazySelections(); // Lazy load selections

  // ... rest of component
}
```

### Step 5: Write Tests (2 hours)

**File:** `tests/services/test_room_service.py`

```python
import pytest
from zndraw.services.room_service import RoomService


def test_room_exists_returns_false_for_new_room(redis_client):
    service = RoomService(redis_client)
    assert service.room_exists("nonexistent") is False


def test_room_exists_returns_true_for_existing_room(redis_client):
    service = RoomService(redis_client)
    redis_client.set("room:test:current_frame", 0)
    assert service.room_exists("test") is True


def test_create_empty_room(redis_client):
    service = RoomService(redis_client)
    result = service.create_room("test", "user1")

    assert result["created"] is True
    assert result["frameCount"] == 0
    assert redis_client.exists("room:test:current_frame")


def test_create_room_from_copy(redis_client):
    service = RoomService(redis_client)

    # Create source room
    redis_client.zadd("room:source:trajectory:indices", {"frame0": 0, "frame1": 1})

    # Copy room
    result = service.create_room("copy", "user1", copy_from="source")

    assert result["created"] is True
    assert result["frameCount"] == 2
    assert redis_client.exists("room:copy:trajectory:indices")


def test_get_frame_count(redis_client):
    service = RoomService(redis_client)
    redis_client.zadd("room:test:trajectory:indices", {"frame0": 0, "frame1": 1, "frame2": 2})

    assert service.get_frame_count("test") == 3
```

**File:** `tests/services/test_client_service.py`

```python
import pytest
from zndraw.services.client_service import ClientService


def test_get_or_create_client_generates_new_id(redis_client):
    service = ClientService(redis_client)
    client_id = service.get_or_create_client("room1", "user1")

    assert client_id is not None
    assert redis_client.hget(f"client:{client_id}", "userName") == "user1"
    assert redis_client.hget(f"client:{client_id}", "currentRoom") == "room1"


def test_create_join_token(redis_client):
    service = ClientService(redis_client)
    token = service.create_join_token("client1", "room1", "user1")

    assert token is not None
    assert redis_client.exists(f"join_token:{token}")
```

## Migration Strategy

### Testing Strategy
1. ✅ Write service tests first
2. ✅ Test new /join endpoint alongside old one
3. ✅ Feature flag for gradual rollout
4. ✅ Monitor response times and error rates
5. ✅ A/B test with subset of users

### Rollout Plan
1. **Week 1:** Implement service layer + tests
2. **Week 2:** Create new /join endpoint at `/api/rooms/<id>/join/v2`
3. **Week 3:** Add lazy loading endpoints + frontend hooks
4. **Week 4:** A/B test with 10% of users
5. **Week 5:** Roll out to 100% if metrics look good
6. **Week 6:** Remove old endpoint, clean up

### Success Metrics
- **Response Time:** /join endpoint < 100ms (down from ~500ms)
- **Payload Size:** < 5KB (down from ~50KB+)
- **Test Coverage:** > 90% for service layer
- **Error Rate:** < 0.1%
- **User Experience:** No increase in reported issues

## Benefits

### Performance
- ✅ 90% reduction in /join response size
- ✅ 80% faster join times
- ✅ Reduced Redis load (fewer operations per join)
- ✅ Better caching (smaller, focused responses)

### Maintainability
- ✅ SOLID principles followed
- ✅ Service layer is testable in isolation
- ✅ Clear separation of concerns
- ✅ Each function < 50 lines
- ✅ Easy to extend and modify

### Scalability
- ✅ Lazy loading reduces initial payload
- ✅ Can cache geometries/bookmarks independently
- ✅ Can add rate limiting per endpoint
- ✅ Can optimize specific endpoints without affecting others

### Developer Experience
- ✅ Easier to understand and modify
- ✅ Easier to test (unit + integration)
- ✅ Easier to debug (smaller functions)
- ✅ Self-documenting code

## Risks & Mitigations

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| Breaking existing clients | High | Medium | Version endpoint, gradual rollout |
| Race conditions in lazy loading | Medium | Low | Proper state management, React Query |
| Increased number of requests | Low | Medium | Batch requests, use HTTP/2 |
| Regression bugs | High | Low | Comprehensive tests, A/B testing |
| Performance degradation | Medium | Low | Monitor metrics, rollback plan |

## Appendix: File Changes Summary

### New Files
- `src/zndraw/services/room_service.py` (~200 lines)
- `src/zndraw/services/client_service.py` (~80 lines)
- `tests/services/test_room_service.py` (~100 lines)
- `tests/services/test_client_service.py` (~50 lines)
- `app/src/hooks/useLazyGeometries.ts` (~40 lines)
- `app/src/hooks/useLazySelections.ts` (~40 lines)

### Modified Files
- `src/zndraw/app/routes.py`: Replace 236 lines with ~50 lines, add 4 new endpoints
- `app/src/hooks/useRestManager.ts`: Simplify join logic (~50 lines removed)
- `app/src/pages/room.tsx`: Add lazy loading hooks (~5 lines)

### Deleted Code
- ~200 lines from routes.py (moved to services)

**Total LOC Change:** +500 new, -200 deleted = +300 net (mostly tests)
