## Current State

**File:** `src/zndraw/app/routes.py:2491+`

**Current Flow:**
1. Client calls `/api/login` with `userName` → gets JWT token + clientId
2. Client stores JWT in localStorage
3. Client sends JWT in `Authorization: Bearer <token>` header to `/join`
4. Server extracts `clientId` and `userName` from JWT
5. Client connects to socket with JWT in auth payload

**Problems (Still Apply):**
- 200+ lines doing too much (violates SRP)
- Returns MASSIVE payload on every join
- 15+ separate Redis operations in single request
- No lazy loading - fetches ALL data upfront
- Hard to test and maintain

## Target Architecture (JWT-Aware)

### Phase 1: Extract Service Layer

**New Structure:**
```
src/zndraw/
├── app/
│   └── routes.py (thin controller - 20 lines)
├── services/
│   ├── room_service.py (NEW - room lifecycle)
│   └── client_service.py (SIMPLIFIED - no token management)
├── auth.py (EXISTING - JWT utilities)
```

**Service Responsibilities:**
- `RoomService`: Create, validate, initialize rooms
- `ClientService`: Manage client metadata in Redis (simplified - no token creation)
- `auth.py`: JWT creation, validation, extraction (already exists)

### Phase 2: Minimize Join Response

**Principle:** Return ONLY what's needed to establish connection

**New Response (90% smaller):**
```json
{
  "status": "ok",
  "roomId": "room-name",
  "clientId": "uuid-from-jwt",
  "frameCount": 0,
  "step": 0,
  "created": true
}
```

**No More:**
- ❌ `joinToken` - JWT handles authentication
- ❌ `geometries` - fetch lazily
- ❌ `bookmarks` - fetch lazily
- ❌ `settings` - fetch lazily
- ❌ `selections` - fetch lazily
- ❌ `selectionGroups` - fetch lazily

## Implementation Steps

### Step 1: Create Service Layer (2 hours)

**File:** `src/zndraw/services/room_service.py`

```python
"""Room lifecycle and data management."""

import json
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
        description: str | None = None,
        copy_from: str | None = None,
    ) -> dict:
        """Create new room, optionally copying from existing room.

        Parameters
        ----------
        room_id : str
            Unique room identifier
        user_name : str
            Username creating the room (for default settings)
        description : str | None
            Optional room description
        copy_from : str | None
            Optional source room to copy from

        Returns
        -------
        dict
            {"created": bool, "frameCount": int}
        """
            # Validate room ID
        if ":" in room_id:
            return {"error": "Room ID cannot contain ':' character", "code": "INVALID_ROOM_ID"}, 400
        # reject any non-string/int room IDs
        if not re.match(r"^[a-zA-Z0-9_\-]+$", room_id):
            return {"error": "Room ID contains invalid characters", "code": "INVALID_ROOM_ID"}, 400

        if copy_from:
            return self._create_room_from_copy(room_id, copy_from)
        return self._create_empty_room(room_id, user_name, description)

    def _create_empty_room(
        self, room_id: str, user_name: str, description: str | None = None
    ) -> dict:
        """Create empty room with defaults.

        Note: User-specific settings are NOT initialized here (YAGNI).
        Settings will be initialized lazily when first accessed.
        """
        # Set description
        if description:
            self.r.set(f"room:{room_id}:description", description)

        # Initialize metadata
        self.r.set(f"room:{room_id}:current_frame", 0)
        self.r.set(f"room:{room_id}:locked", 0)
        self.r.set(f"room:{room_id}:hidden", 0)

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
                {member: score for member, score in source_indices},
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

        log.info(
            f"Created room '{room_id}' from '{source_room}' with {len(source_indices)} frames"
        )
        return {"created": True, "frameCount": len(source_indices)}

    def _initialize_default_geometries(self, room_id: str):
        """Initialize default geometries for new room."""
        from zndraw.geometries import Bond, Cell, Curve, Floor, Sphere

        defaults = {
            "particles": (
                Sphere,
                {"position": "arrays.positions", "color": "arrays.colors"},
            ),
            "bonds": (Bond, {}),
            "curve": (Curve, {}),
            "cell": (Cell, {}),
            "floor": (Floor, {}),
        }

        for key, (geometry_class, kwargs) in defaults.items():
            geometry_data = geometry_class(**kwargs).model_dump()
            self.r.hset(
                f"room:{room_id}:geometries",
                key,
                json.dumps({"type": geometry_class.__name__, "data": geometry_data}),
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
                    log.warning(
                        f"Negative frame in Redis for room {room_id}: {step_int}, resetting to 0"
                    )
                    return 0
                return step_int
            return 0
        except (ValueError, TypeError) as e:
            log.error(f"Invalid step value in Redis for room {room_id}: {step} - {e}")
            return 0
```

**File:** `src/zndraw/services/client_service.py`

```python
"""Client session metadata management (JWT-aware).

Note: JWT tokens are created/validated by auth.py.
This service only manages client metadata in Redis.
"""

import logging

from redis import Redis

log = logging.getLogger(__name__)


class ClientService:
    """Handles client session metadata in Redis.

    JWT authentication is handled separately by auth.py.
    This service focuses on storing client metadata like currentRoom.
    """

    def __init__(self, redis_client: Redis):
        self.r = redis_client

    def update_client_room(self, client_id: str, room_id: str) -> None:
        """Update client's current room in Redis.

        Parameters
        ----------
        client_id : str
            Client identifier (from JWT sub claim)
        room_id : str
            Room the client is joining
        """
        client_key = f"client:{client_id}"
        self.r.hset(client_key, "currentRoom", room_id)
        log.info(f"Client {client_id} updated room to {room_id}")

    def add_client_to_room(self, room_id: str, client_id: str) -> None:
        """Add client to room's client set.

        Parameters
        ----------
        room_id : str
            Room identifier
        client_id : str
            Client identifier
        """
        self.r.sadd(f"room:{room_id}:clients", client_id)

    def remove_client_from_room(self, room_id: str, client_id: str) -> None:
        """Remove client from room's client set.

        Parameters
        ----------
        room_id : str
            Room identifier
        client_id : str
            Client identifier
        """
        self.r.srem(f"room:{room_id}:clients", client_id)

    def get_room_clients(self, room_id: str) -> set[str]:
        """Get all client IDs currently in a room.

        Parameters
        ----------
        room_id : str
            Room identifier

        Returns
        -------
        set[str]
            Set of client IDs
        """
        members = self.r.smembers(f"room:{room_id}:clients")
        return {m for m in members}
```

### Step 2: Refactor Join Route (30 minutes)

**File:** `src/zndraw/app/routes.py`

Replace existing join_room function with:

```python
from flask import current_app, request

from zndraw.auth import AuthError, get_current_client
from zndraw.services.client_service import ClientService
from zndraw.services.room_service import RoomService


@main.route("/api/rooms/<string:room_id>/join", methods=["POST"])
def join_room(room_id: str):
    """Join a room - returns minimal data for connection establishment.

    JWT authentication required via Authorization header.
    All other data (geometries, bookmarks, settings) should be fetched lazily.

    Headers
    -------
    Authorization: Bearer <jwt-token> (required)

    Request
    -------
    {
        "description": "optional room description",
        "copyFrom": "optional-source-room",
        "allowCreate": true
    }

    Response (Success)
    ------------------
    {
        "status": "ok",
        "roomId": "room-name",
        "clientId": "uuid-from-jwt",
        "frameCount": 0,
        "step": 0,
        "created": true
    }

    Response (Error)
    ----------------
    {
        "error": "error message",
        "code": "ERROR_CODE"
    }
    """
    data = request.get_json() or {}

    # Extract client info from JWT (authentication happens here)
    try:
        client = get_current_client()
        client_id = client["clientId"]
        user_name = client["userName"]
    except AuthError as e:
        return {"error": e.message, "code": "AUTH_ERROR"}, e.status_code

    # Extract parameters
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
                "error": f"Room '{room_id}' does not exist",
                "code": "ROOM_NOT_FOUND",
            }, 404

        try:
            room_service.create_room(room_id, user_name, description, copy_from)
        except ValueError as e:
            return {"error": str(e), "code": "ROOM_CREATION_FAILED"}, 400

    # Update client's room membership
    client_service.update_client_room(client_id, room_id)
    client_service.add_client_to_room(room_id, client_id)

    # Return minimal response
    return {
        "status": "ok",
        "roomId": room_id,
        "clientId": client_id,
        "frameCount": room_service.get_frame_count(room_id),
        "step": room_service.get_current_frame(room_id),
        "created": not room_exists,
    }
```

### Step 3: Create Lazy Loading Endpoints (2 hours)

**Add new JWT-authenticated endpoints for lazy data fetching:**

**Note:** Since `@require_auth` may not exist yet, we'll use `get_current_client()` directly for authentication.

```python
import json

from flask import current_app

from zndraw.auth import AuthError, get_current_client
from zndraw.services.client_service import ClientService
from zndraw.services.room_service import RoomService


def _verify_room_access(room_id: str, client_id: str) -> tuple[dict, int] | None:
    """Verify client has access to the room.

    Returns
    -------
    tuple[dict, int] | None
        Error response tuple if access denied, None if access granted
    """
    r = current_app.extensions["redis"]
    room_service = RoomService(r)
    client_service = ClientService(r)

    # Check room exists
    if not room_service.room_exists(room_id):
        return {"error": "Room not found", "code": "ROOM_NOT_FOUND"}, 404

    # Check client is in room
    room_clients = client_service.get_room_clients(room_id)
    if client_id not in room_clients:
        return {"error": "Access denied - not a member of this room", "code": "ACCESS_DENIED"}, 403

    return None


@main.route("/api/rooms/<string:room_id>/geometries", methods=["GET"])
def get_geometries(room_id: str):
    """Fetch all geometries for a room (requires JWT).

    Headers
    -------
    Authorization: Bearer <jwt-token>

    Response (Success)
    ------------------
    {
        "geometries": {...},
        "geometryDefaults": {...}
    }

    Response (Error)
    ----------------
    {
        "error": "error message",
        "code": "ERROR_CODE"
    }
    """
    # Authenticate
    try:
        client = get_current_client()
        client_id = client["clientId"]
    except AuthError as e:
        return {"error": e.message, "code": "AUTH_ERROR"}, e.status_code

    # Verify room access
    access_error = _verify_room_access(room_id, client_id)
    if access_error:
        return access_error

    r = current_app.extensions["redis"]
    geometries_raw = r.hgetall(f"room:{room_id}:geometries")

    return {
        "geometries": {k: json.loads(v) for k, v in geometries_raw.items()} if geometries_raw else {},
    }


@main.route("/api/rooms/<string:room_id>/bookmarks", methods=["GET"])
def get_bookmarks(room_id: str):
    """Fetch bookmarks for a room (requires JWT).

    Headers
    -------
    Authorization: Bearer <jwt-token>
    """
    # Authenticate
    try:
        client = get_current_client()
        client_id = client["clientId"]
    except AuthError as e:
        return {"error": e.message, "code": "AUTH_ERROR"}, e.status_code

    # Verify room access
    access_error = _verify_room_access(room_id, client_id)
    if access_error:
        return access_error

    r = current_app.extensions["redis"]
    bookmarks_raw = r.hgetall(f"room:{room_id}:bookmarks")

    return {
        "bookmarks": {int(k): v for k, v in bookmarks_raw.items()} if bookmarks_raw else {},
    }


@main.route("/api/rooms/<string:room_id>/selections", methods=["GET"])
def get_selections(room_id: str):
    """Fetch selections and selection groups (requires JWT).

    Headers
    -------
    Authorization: Bearer <jwt-token>
    """
    # Authenticate
    try:
        client = get_current_client()
        client_id = client["clientId"]
    except AuthError as e:
        return {"error": e.message, "code": "AUTH_ERROR"}, e.status_code

    # Verify room access
    access_error = _verify_room_access(room_id, client_id)
    if access_error:
        return access_error

    r = current_app.extensions["redis"]

    # Per-geometry selections
    selections_raw = r.hgetall(f"room:{room_id}:selections")
    selections = {k: json.loads(v) for k, v in selections_raw.items()} if selections_raw else {}

    # Selection groups
    groups_raw = r.hgetall(f"room:{room_id}:selection_groups")
    selection_groups = {k: json.loads(v) for k, v in groups_raw.items()} if groups_raw else {}

    # Active group
    active_group = r.get(f"room:{room_id}:active_selection_group")

    # Frame selection
    frame_selection_raw = r.get(f"room:{room_id}:frame_selection:default")
    frame_selection = json.loads(frame_selection_raw) if frame_selection_raw else None

    return {
        "selections": selections,
        "selectionGroups": selection_groups,
        "activeSelectionGroup": active_group,
        "frameSelection": frame_selection,
    }


@main.route("/api/rooms/<string:room_id>/settings", methods=["GET"])
def get_settings(room_id: str):
    """Fetch user settings for a room (requires JWT).

    If settings don't exist, initializes them with defaults (lazy initialization).

    Headers
    -------
    Authorization: Bearer <jwt-token>
    """
    # Authenticate
    try:
        client = get_current_client()
        client_id = client["clientId"]
        user_name = client["userName"]
    except AuthError as e:
        return {"error": e.message, "code": "AUTH_ERROR"}, e.status_code

    # Verify room access
    access_error = _verify_room_access(room_id, client_id)
    if access_error:
        return access_error

    r = current_app.extensions["redis"]
    settings_key = f"room:{room_id}:user:{user_name}:settings"
    settings_hash = r.hgetall(settings_key)

    # Lazy initialization: create defaults if not exists
    if not settings_hash:
        from zndraw.settings import RoomConfig

        default_config = RoomConfig()
        config_dict = default_config.model_dump()

        for category_name, category_data in config_dict.items():
            r.hset(settings_key, category_name, json.dumps(category_data))

        settings_hash = r.hgetall(settings_key)

    settings_data = {k: json.loads(v) for k, v in settings_hash.items()}
    return {"settings": settings_data}
```

### Step 4: Update Frontend (2-3 hours)

**File:** `app/src/hooks/useRestManager.ts`

Update to ensure JWT token is included in requests:

```typescript
import { useCallback, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import { useAppStore } from '../store';
import { getToken } from '../utils/auth';

export const useRestJoinManager = () => {
  const {
    setClientId,
    setRoomId,
    setUserId,
    setCurrentFrame,
    setFrameCount,
  } = useAppStore();

  const { roomId: room, userId } = useParams<{
    roomId: string;
    userId: string;
  }>();

  const joinRoom = useCallback(async () => {
    if (!room || !userId) return;

    console.log('Joining room via REST:', room, userId);

    try {
      const token = getToken();
      if (!token) {
        throw new Error('No authentication token available');
      }

      // Minimal join - only essential data
      const response = await fetch(`/api/rooms/${room}/join`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${token}`,
        },
        body: JSON.stringify({
          allowCreate: true,
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.error || `Join failed: ${response.statusText}`);
      }

      const data = await response.json();
      console.log('Join response data:', data);

      // Update only essential state
      if (data.clientId) setClientId(data.clientId);
      if (typeof data.frameCount === 'number') setFrameCount(data.frameCount);
      if (data.step !== undefined && data.step !== null) setCurrentFrame(data.step);

      setRoomId(room);
      setUserId(userId);

      // Everything else loaded lazily by other hooks/components

    } catch (error) {
      console.error('Error joining room:', error);
      throw error;
    }
  }, [room, userId, setClientId, setRoomId, setUserId, setCurrentFrame, setFrameCount]);

  useEffect(() => {
    joinRoom();
  }, [joinRoom]);

  return { joinRoom };
};
```

**Add new hooks for lazy loading with JWT:**

```typescript
// app/src/hooks/useLazyGeometries.ts
import { useEffect } from 'react';
import { useQuery } from '@tanstack/react-query';
import { useAppStore } from '../store';
import { getToken } from '../utils/auth';

async function fetchGeometries(roomId: string) {
  const token = getToken();
  if (!token) {
    throw new Error('No authentication token available');
  }

  const response = await fetch(`/api/rooms/${roomId}/geometries`, {
    headers: {
      'Authorization': `Bearer ${token}`,
    },
  });

  if (!response.ok) {
    const errorData = await response.json();
    throw new Error(errorData.error || 'Failed to fetch geometries');
  }

  return response.json();
}

export const useLazyGeometries = () => {
  const { roomId, setGeometries } = useAppStore();

  const { data, isLoading, error } = useQuery({
    queryKey: ['geometries', roomId],
    queryFn: () => fetchGeometries(roomId),
    enabled: !!roomId,
    staleTime: Infinity, // Cache forever unless invalidated
  });

  useEffect(() => {
    if (data?.geometries) {
      setGeometries(data.geometries);
    }
  }, [data, setGeometries]);

  return { isLoading, error };
};

// app/src/hooks/useLazySelections.ts
import { useEffect } from 'react';
import { useQuery } from '@tanstack/react-query';
import { useAppStore } from '../store';
import { getToken } from '../utils/auth';

async function fetchSelections(roomId: string) {
  const token = getToken();
  if (!token) {
    throw new Error('No authentication token available');
  }

  const response = await fetch(`/api/rooms/${roomId}/selections`, {
    headers: {
      'Authorization': `Bearer ${token}`,
    },
  });

  if (!response.ok) {
    const errorData = await response.json();
    throw new Error(errorData.error || 'Failed to fetch selections');
  }

  return response.json();
}

export const useLazySelections = () => {
  const {
    roomId,
    setSelections,
    setSelectionGroups,
    setActiveSelectionGroup,
    setFrameSelection,
  } = useAppStore();

  const { data, error } = useQuery({
    queryKey: ['selections', roomId],
    queryFn: () => fetchSelections(roomId),
    enabled: !!roomId,
    staleTime: 60000, // Cache for 1 minute
  });

  useEffect(() => {
    if (data) {
      if (data.selections) setSelections(data.selections);
      if (data.selectionGroups) setSelectionGroups(data.selectionGroups);
      if (data.activeSelectionGroup !== undefined)
        setActiveSelectionGroup(data.activeSelectionGroup);
      if (data.frameSelection !== undefined)
        setFrameSelection(data.frameSelection);
    }
  }, [data, setSelections, setSelectionGroups, setActiveSelectionGroup, setFrameSelection]);

  return { error };
};

// app/src/hooks/useLazySettings.ts
import { useEffect } from 'react';
import { useQuery } from '@tanstack/react-query';
import { useAppStore } from '../store';
import { getToken } from '../utils/auth';

async function fetchSettings(roomId: string) {
  const token = getToken();
  if (!token) {
    throw new Error('No authentication token available');
  }

  const response = await fetch(`/api/rooms/${roomId}/settings`, {
    headers: {
      'Authorization': `Bearer ${token}`,
    },
  });

  if (!response.ok) {
    const errorData = await response.json();
    throw new Error(errorData.error || 'Failed to fetch settings');
  }

  return response.json();
}

export const useLazySettings = () => {
  const { roomId, setSettings } = useAppStore();

  const { data, error } = useQuery({
    queryKey: ['settings', roomId],
    queryFn: () => fetchSettings(roomId),
    enabled: !!roomId,
    staleTime: 300000, // Cache for 5 minutes
  });

  useEffect(() => {
    if (data?.settings) {
      setSettings(data.settings);
    }
  }, [data, setSettings]);

  return { error };
};
```

### Step 5: Write Tests (2 hours)

**File:** `tests/conftest.py` (Add fixture if not exists)

```python
import pytest
from redis import Redis


@pytest.fixture
def redis_client():
    """Provide a Redis client for testing.

    Note: Assumes Redis is running on localhost:6379.
    Consider using fakeredis for unit tests if needed.
    """
    client = Redis(host="localhost", port=6379, db=15, decode_responses=True)
    yield client
    # Cleanup: flush test database
    client.flushdb()
```

**File:** `tests/services/test_room_service.py`

```python
import pytest

from zndraw.services.room_service import RoomService


def test_room_exists_returns_false_for_new_room(redis_client):
    """Test room_exists returns False for nonexistent room."""
    service = RoomService(redis_client)
    assert service.room_exists("nonexistent") is False


def test_room_exists_returns_true_for_existing_room(redis_client):
    """Test room_exists returns True for existing room."""
    service = RoomService(redis_client)
    redis_client.set("room:test:current_frame", 0)
    assert service.room_exists("test") is True


def test_create_empty_room(redis_client):
    """Test creating an empty room with defaults."""
    service = RoomService(redis_client)
    result = service.create_room("test", "user1")

    assert result["created"] is True
    assert result["frameCount"] == 0
    assert redis_client.exists("room:test:current_frame")
    # Verify default geometries created
    assert redis_client.exists("room:test:geometries")


def test_create_empty_room_with_description(redis_client):
    """Test creating room with description."""
    service = RoomService(redis_client)
    result = service.create_room("test", "user1", description="Test Room")

    assert result["created"] is True
    assert redis_client.get("room:test:description") == "Test Room"


def test_create_room_from_copy(redis_client):
    """Test creating room by copying from existing room."""
    service = RoomService(redis_client)

    # Create source room
    redis_client.zadd("room:source:trajectory:indices", {"frame0": 0, "frame1": 1})

    # Copy room
    result = service.create_room("copy", "user1", copy_from="source")

    assert result["created"] is True
    assert result["frameCount"] == 2
    assert redis_client.exists("room:copy:trajectory:indices")


def test_create_room_from_nonexistent_source_raises_error(redis_client):
    """Test copying from nonexistent room raises ValueError."""
    service = RoomService(redis_client)

    with pytest.raises(ValueError, match="Source room 'nonexistent' not found"):
        service.create_room("test", "user1", copy_from="nonexistent")


def test_get_frame_count(redis_client):
    """Test getting frame count from room."""
    service = RoomService(redis_client)
    redis_client.zadd(
        "room:test:trajectory:indices", {"frame0": 0, "frame1": 1, "frame2": 2}
    )

    assert service.get_frame_count("test") == 3


def test_get_frame_count_empty_room(redis_client):
    """Test getting frame count from room with no frames."""
    service = RoomService(redis_client)
    assert service.get_frame_count("test") == 0


@pytest.mark.parametrize(
    "redis_value,expected",
    [
        (0, 0),
        (5, 5),
        (None, 0),
        (-1, 0),
        ("invalid", 0),
    ],
)
def test_get_current_frame_handles_invalid_values(redis_client, redis_value, expected):
    """Test get_current_frame handles various invalid values."""
    service = RoomService(redis_client)

    if redis_value is not None:
        redis_client.set("room:test:current_frame", redis_value)

    assert service.get_current_frame("test") == expected
```

**File:** `tests/services/test_client_service.py`

```python
from zndraw.services.client_service import ClientService


def test_update_client_room(redis_client):
    """Test updating client's current room."""
    service = ClientService(redis_client)
    service.update_client_room("client1", "room1")

    assert redis_client.hget("client:client1", "currentRoom") == "room1"


def test_update_client_room_overwrites_previous(redis_client):
    """Test updating client room overwrites previous value."""
    service = ClientService(redis_client)
    service.update_client_room("client1", "room1")
    service.update_client_room("client1", "room2")

    assert redis_client.hget("client:client1", "currentRoom") == "room2"


def test_add_client_to_room(redis_client):
    """Test adding client to room's client set."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")

    assert redis_client.sismember("room:room1:clients", "client1")


def test_add_client_to_room_idempotent(redis_client):
    """Test adding same client twice is idempotent."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")
    service.add_client_to_room("room1", "client1")

    assert redis_client.scard("room:room1:clients") == 1


def test_remove_client_from_room(redis_client):
    """Test removing client from room's client set."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")
    service.remove_client_from_room("room1", "client1")

    assert not redis_client.sismember("room:room1:clients", "client1")


def test_remove_client_from_room_idempotent(redis_client):
    """Test removing nonexistent client is idempotent."""
    service = ClientService(redis_client)
    service.remove_client_from_room("room1", "client1")  # Should not raise

    assert not redis_client.sismember("room:room1:clients", "client1")


def test_get_room_clients(redis_client):
    """Test getting all clients in a room."""
    service = ClientService(redis_client)
    service.add_client_to_room("room1", "client1")
    service.add_client_to_room("room1", "client2")

    clients = service.get_room_clients("room1")
    assert "client1" in clients
    assert "client2" in clients
    assert len(clients) == 2


def test_get_room_clients_empty_room(redis_client):
    """Test getting clients from room with no clients."""
    service = ClientService(redis_client)

    clients = service.get_room_clients("room1")
    assert len(clients) == 0
```

## Migration Strategy

**Note:** Per CLAUDE.md guidelines, this is a new application - no backwards compatibility required.

### Implementation Phases

#### Phase 1: Service Layer (Week 1)
- Implement `RoomService` and `ClientService`
- Write comprehensive unit tests
- Target: >90% test coverage
- Ensure all tests pass before proceeding

#### Phase 2: Backend Refactoring (Week 2)
- **Replace** existing /join endpoint (no compatibility layer)
- Add lazy loading endpoints
- Update all error responses to use consistent format
- Integration tests with JWT authentication

#### Phase 3: Frontend Migration (Week 3)
- Update `useRestManager` to use new minimal join response
- Implement lazy loading hooks
- Remove dead code that handled old response format
- Update error handling to parse new error format
