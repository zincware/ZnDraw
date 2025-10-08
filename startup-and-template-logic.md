# Current State (Before Changes)
One or multiple files are uploaded to a named room each.
Once they are uploaded, the rooms are locked and defined as templates.
The first room finished uploading is defined as "default", until then, default is "empty".
Rooms can be modified if not locked, e.g. frames added. But the underlying data is not modified, just new data appended and the mapping is updated from room to data.

**Issues with Current Design:**
- Template concept conflates two concerns: reusability and immutability
- Opening the app too early might show an empty room or only the default room
- Room locking mechanism is permanent and uses trajectory locks
- No explicit room duplication mechanism
- Auto-promotion to template happens during file upload
- First file automatically becomes default template

# New Architecture Design

## Core Principles
1. **Separate Concerns**: Rooms can be locked (immutable) OR used as templates (duplicatable), but these are independent properties
2. **Explicit Actions**: Users explicitly set default rooms and duplicate existing rooms
3. **No Auto-Magic**: File uploads don't automatically promote rooms to templates or set defaults
4. **Persistent UI**: Always show room list unless explicitly navigating to a specific room
5. **Room Metadata**: Track room properties (locked, hidden, frame count) for better UX

## Data Model Changes

### Room Properties (Redis Keys)
- `room:{room_id}:description` - Optional text description for the room
- `room:{room_id}:locked` - Boolean flag (0 or 1) indicating if room is locked (immutable)
- `room:{room_id}:hidden` - Boolean flag (0 or 1) indicating if room should be hidden from list
- `default_room` - Redis key storing the default room_id (replaces `default_template`)
- Remove: `room:{room_id}:template` field (we don't track template source)
- Remove: Permanent trajectory locks for templates (use the new locked flag instead)

## Startup Logic (Root Path `/`)

```
START -> Check room count
         |
         v
    No rooms? -> Navigate to `/rooms/{new-uuid}/{user-uuid}?template=empty`
         |
         v
    One room? -> Navigate to `/rooms/{room_id}/{user-uuid}`
         |
         v
    Multiple rooms?
         |
         +-> Has default_room set? -> Navigate to `/rooms/{default_room}/{user-uuid}`
         |
         +-> No default? -> Navigate to `/rooms` (show list)
```

## Existing API Structure Analysis

### Current Room-Related Endpoints

**Existing (Keep & Enhance):**
- `GET /api/rooms` - Returns list of rooms with id and template
- `GET /api/rooms/{room_id}` - Returns room details (id, template)
- `POST /api/rooms/{room_id}/join` - Join a room, create if doesn't exist
- `GET /api/rooms/{room_id}/frames` - Get frames data
- `POST /api/rooms/{room_id}/frames` - Append/extend/replace/insert frames
- `DELETE /api/rooms/{room_id}/frames` - Delete frames
- `PATCH /api/rooms/{room_id}/frames/bulk` - Bulk replace frames
- `GET /api/rooms/{room_id}/frames/{frame_id}/metadata` - Get frame metadata
- `POST /api/rooms/{room_id}/geometries` - Create/update geometry
- `GET /api/rooms/{room_id}/geometries` - List all geometries
- `GET /api/rooms/{room_id}/geometries/{key}` - Get specific geometry
- `DELETE /api/rooms/{room_id}/geometries/{key}` - Delete geometry
- `POST /api/rooms/{room_id}/figures` - Create/update figure
- `GET /api/rooms/{room_id}/figures` - List all figures
- `GET /api/rooms/{room_id}/figures/{key}` - Get specific figure
- `DELETE /api/rooms/{room_id}/figures/{key}` - Delete figure
- `GET /api/rooms/{room_id}/chat/messages` - Get chat messages
- `GET /api/rooms/{room_id}/schema/{category}` - Get extension schema
- `POST /api/rooms/{room_id}/extensions/register` - Register extension
- `GET /api/rooms/{room_id}/jobs` - List jobs
- `POST /api/rooms/{room_id}/jobs/next` - Get next job for worker
- `POST /api/rooms/{room_id}/renormalize` - Renormalize frame indices

**Template-Related (To Modify/Remove):**
- `GET /api/templates` - List templates
- `GET /api/templates/default` - Get default template
- `PUT /api/templates/default` - Set default template
- `POST /api/rooms/{room_id}/promote` - Promote room to template (locks it)

### Current Redis Key Patterns

**Per-Room Keys (from code analysis):**
```
room:{room_id}:template             - Template ID (string)
room:{room_id}:trajectory:indices   - Sorted set of frame mapping
room:{room_id}:current_frame        - Current frame number (string/int)
room:{room_id}:selection:default    - JSON array of selected atoms
room:{room_id}:frame_selection:default - JSON array of selected frames
room:{room_id}:presenter_lock       - Presenter SID (string, TTL)
room:{room_id}:lock:{target}        - Distributed lock (string, TTL)
room:{room_id}:bookmarks            - Hash: physical_key -> label
room:{room_id}:geometries           - Hash: key -> JSON geometry data
room:{room_id}:figures              - Hash: key -> JSON figure data
room:{room_id}:clients              - Set of client IDs in room
room:{room_id}:user:{user}:settings - Hash: category -> JSON settings
room:{room_id}:chat:counter         - Message counter (int)
room:{room_id}:chat:index           - Sorted set: message_id -> timestamp
room:{room_id}:chat:data            - Hash: message_id -> JSON message
room:{room_id}:extensions:{category} - Hash: extension_name -> JSON schema
room:{room_id}:extensions:{category}:{extension}:idle_workers - Set of worker IDs
room:{room_id}:extensions:{category}:{extension}:progressing_workers - Set of worker IDs
room:{room_id}:extensions:{category}:{extension}:queue - List of job IDs
```

**Global Keys:**
```
default_template        - Default template ID (string)
template:{template_id}  - Hash with id, name, description
client:{client_id}      - Hash with userName, currentRoom, currentSid
sid:{sid}               - Client ID for socket ID (string)
join_token:{token}      - Token data JSON (string, TTL 60s)
job:{job_id}            - Hash with job data
```

## API Changes

### Modified Endpoints

```json
[
  {
    "id": "room1",
    "description": "Initial dataset from file.xyz",
    "frameCount": 42,
    "locked": false,
    "hidden": false,
    "isDefault": false
  },
  {
    "id": "room2",
    "description": "Modified copy for testing",
    "frameCount": 42,
    "locked": true,
    "hidden": false,
    "isDefault": true
  }
]
```

**GET /api/rooms/{room_id}** (Enhanced)
**Current:** Returns `{"id": "room1", "template": "empty"}`
**New:** Add description, frameCount, locked, hidden flags
```json
{
  "id": "room1",
  "description": "My room",
  "frameCount": 42,
  "locked": false,
  "hidden": false
}
```

**PATCH /api/rooms/{room_id}** (New)
Update room metadata (description, locked, hidden).
```json
{
  "description": "My custom description",  // Optional
  "locked": true,                          // Optional
  "hidden": false                          // Optional
}
```
Response: `{"status": "ok"}`

**GET /api/rooms/default** (New)
Get the default room ID.
Response: `{"roomId": "room1"}` or `{"roomId": null}`

**PUT /api/rooms/default** (New - replaces /api/templates/default)
Set the default room.
```json
{
  "roomId": "room1"  // or null to unset
}
```
Response: `{"status": "ok"}`

**POST /api/rooms/{room_id}/duplicate** (New)
Duplicate a room by copying all frame mappings.
```json
{
  "newRoomId": "new-room-uuid",  // Optional, auto-generated if not provided
  "description": "Copy of room1"  // Optional, description for new room
}
```
Response:
```json
{
  "status": "ok",
  "roomId": "new-room-uuid",
  "frameCount": 42
}
```

**Endpoints to Deprecate/Remove:**
- `GET /api/templates` - No longer needed (use GET /api/rooms)
- `GET /api/templates/default` - Replaced by GET /api/rooms/default
- `PUT /api/templates/default` - Replaced by PUT /api/rooms/default
- `POST /api/rooms/{room_id}/promote` - Locking is now via PATCH /api/rooms/{room_id}
**Current:** Returns `[{"id": "room1", "template": "empty"}]`
**New:** Add frameCount, locked, hidden, isDefault flags
```json
[
  {
    "id": "room1",
    "description": "Initial dataset from file.xyz",
    "frameCount": 42,
    "locked": false,
    "hidden": false,
    "isDefault": false
  },
  {
    "id": "room2",
    "description": "Modified copy for testing",
    "frameCount": 42,
    "locked": true,
    "hidden": false,
    "isDefault": true
  }
]
```

**PUT /api/rooms/{room_id}/description** (New)
Set or update the room description.
```json
{
  "description": "My custom description"  // or null to clear
}
```

**GET /api/rooms/default** (New)
Returns the default room_id or null if not set.

**POST /api/rooms/{room_id}/duplicate** (New)
Creates a new room by copying all frames from the source room.
```json
{
  "newRoomId": "new-room-uuid",  // Optional, auto-generated if not provided
  "description": "Copy of room1"  // Optional, description for new room
}
```
Response:
```json
{
  "status": "ok",
  "newRoomId": "new-room-uuid",
  "frameCount": 42
}
```

**PUT /api/rooms/{room_id}/lock** (New)
Lock/unlock a room to make it immutable.
```json
{
  "locked": true
}
```

**PUT /api/rooms/{room_id}/hide** (New)
Hide/show a room in the room list.
```json
{
  "hidden": true
}
```

**PUT /api/rooms/default** (Modified)
Set the default room (replaces /api/templates/default).
```json
{
  "room_id": "room1"  // or null to unset
}
```

**POST /api/rooms/{room_id}/promote** (Deprecated)
This endpoint should be removed.

## Frontend Changes

### Route Structure
- `/` - Startup logic page (implements logic above)
- `/rooms` - Room list page (MUI DataGrid)
- `/rooms/{roomId}/{userId}` - Room view (3D viewer + controls)

### Room List Page (`/rooms`)

**Components:**
- MUI-X DataGrid with columns:
  - Room ID/Name
  - Description (editable inline)
  - Frame Count
  - Locked Status (icon/badge)
  - Hidden Status (toggle)
  - Default Indicator (star icon)
  - Actions (buttons)

**Action Buttons:**
- Enter Room - Navigate to room view
- Duplicate - Create new room from this one
- Lock/Unlock - Toggle immutability
- Set as Default - Mark as default room
- Hide/Show - Toggle visibility in list  (use mui-x feature)

**Behavior:**
- Always show list, even if only one room exists
- Filter out hidden rooms by default (with toggle to show all) (use mui-x feature)
- Highlight default room visually

### Room View Changes

**Lower Right Corner Panel:**
- Lock Status Indicator (icon + text)
- Lock/Unlock Button
- Duplicate Room Button
- Hide Room Button
- Back to Room List Button

## Backend Implementation Details

### Room Duplication Logic
```python
def duplicate_room(source_room_id, new_room_id, description=None):
    """Duplicate a room by copying all mappings and metadata.
    
    Frame data is shared - only the mappings are copied.
    This is efficient as frames are referenced by room_id:physical_index.
    """
    r = current_app.extensions["redis"]
    
    # 1. Check source room exists
    if not r.exists(f"room:{source_room_id}:trajectory:indices"):
        return {"error": "Source room not found"}, 404
    
    # 2. Copy trajectory indices (sorted set) - this shares frame data
    source_indices = r.zrange(
        f"room:{source_room_id}:trajectory:indices", 
        0, -1, 
        withscores=True
    )
    if source_indices:
        r.zadd(
            f"room:{new_room_id}:trajectory:indices",
            {member: score for member, score in source_indices}
        )
    
    # 3. Copy geometries hash
    geometries = r.hgetall(f"room:{source_room_id}:geometries")
    if geometries:
        r.hset(f"room:{new_room_id}:geometries", mapping=geometries)
    
    # 4. Copy bookmarks hash (physical keys remain valid)
    bookmarks = r.hgetall(f"room:{source_room_id}:bookmarks")
    if bookmarks:
        r.hset(f"room:{new_room_id}:bookmarks", mapping=bookmarks)
    
    # 5. Initialize new room metadata
    r.set(f"room:{new_room_id}:current_frame", 0)
    r.set(f"room:{new_room_id}:locked", 0)
    r.set(f"room:{new_room_id}:hidden", 0)
    if description:
        r.set(f"room:{new_room_id}:description", description)
    
    # 6. DO NOT copy:
    #    - template field (we don't track lineage)
    #    - selection/frame_selection (user-specific)
    #    - presenter_lock (session-specific)
    #    - user settings (user-specific)
    #    - chat messages (history is room-specific)
    #    - figures (visualization-specific, can be recreated)
    #    - extensions/workers/jobs (session-specific)
    
    return {
        "status": "ok",
        "roomId": new_room_id,
        "frameCount": len(source_indices)
    }
```

**What Gets Duplicated:**
- ✅ Frame mapping (trajectory:indices) - shares actual data
- ✅ Geometries (rendering configuration)
- ✅ Bookmarks (frame annotations)
- ✅ Description (if provided)

**What Doesn't Get Duplicated:**
- ❌ Template field - we don't track lineage
- ❌ User selections - these are user-specific
- ❌ User settings - these are user-specific
- ❌ Chat history - separate conversation per room
- ❌ Figures - can be regenerated
- ❌ Extensions/workers/jobs - session state
- ❌ Current frame - always starts at 0
- ❌ Locks - always starts unlocked

### Lock Mechanism - Understanding Current System

**Existing Lock Types:**

1. **Trajectory Lock (`trajectory:meta`)** - Used by client context manager
   - Key: `room:{room_id}:lock:trajectory:meta`
   - Purpose: Prevents concurrent frame modifications during bulk operations
   - Acquired via: `with vis.lock:` in Python client
   - TTL: 60 seconds with automatic renewal
   - Implementation: Redis SETNX with expiry
   - Used in: File uploads, bulk frame operations

2. **Presenter Lock** - Used for animation/scrubbing control
   - Key: `room:{room_id}:presenter_lock`
   - Purpose: Grants exclusive control for continuous frame updates (scrubbing)
   - Acquired via: `request_presenter_token` socket event
   - TTL: 5 seconds with heartbeat renewal every 3 seconds
   - Behavior: Blocks `set_frame_atomic` calls from other clients
   - Auto-cleanup: On disconnect or explicit release

3. **Template Lock (Current, TO BE REMOVED)** - Permanent room lock
   - Key: `room:{room_id}:lock:trajectory:meta`
   - Purpose: Makes promoted template rooms permanently immutable
   - Implementation: Redis SET without TTL (permanent)
   - Value: `"template-lock"` (special marker)
   - Problem: Conflates reusability with immutability

**New Lock Strategy:**

Replace the permanent template lock with a simple flag:
- Add `room:{room_id}:locked` boolean flag (0 or 1)
- Check this flag before ANY mutating operation:
  - Frame append/extend/replace/insert/delete
  - Geometry modifications
- Return 403 Forbidden if `locked == 1`
- The existing trajectory lock (`trajectory:meta`) continues to work for:
  - Preventing concurrent modifications during bulk operations
  - File uploads with `with vis.lock:`
- The presenter lock continues to work for animation control

**Migration of Existing Locks:**
```python
# On server startup or migration script:
# 1. Find all permanent trajectory locks (value == "template-lock")
# 2. Set room:{room_id}:locked = 1 for those rooms
# 3. Delete the permanent trajectory lock key
# 4. The room is now properly locked via the flag, not the lock mechanism
```

**Implementation Changes:**
```python
# Before any mutating operation (in routes.py):
def check_room_locked(room_id):
    r = current_app.extensions["redis"]
    if r.get(f"room:{room_id}:locked") == "1":
        return {"error": "Room is locked and cannot be modified"}, 403
    return None

# In append_frame, delete_frames, etc:
@main.route("/api/rooms/<string:room_id>/frames", methods=["POST"])
def append_frame(room_id):
    # Check if room is locked
    lock_error = check_room_locked(room_id)
    if lock_error:
        return lock_error
    
    # Proceed with normal operation...
    # The trajectory lock (via vis.lock) still works here
```

**Benefits of This Approach:**
- Separates concerns: locking vs. synchronization
- Trajectory locks remain for their original purpose (preventing concurrent access)
- Presenter locks remain for animation control
- Room lock is a simple, clear flag for immutability
- No confusion between "template" and "locked"
- Easy to toggle on/off via API

### File Upload Changes (CLI)
```python
# Remove auto-promotion and default setting
# Use explicit lock via vis.lock context manager
with vis.lock:
    for batch in batches:
        vis.extend(batch)
```

## Migration Strategy

### Phase 1: Add New Endpoints (Non-Breaking)
- Implement new endpoints alongside existing ones
- Add room metadata fields (locked, hidden)
- Keep existing template system working

### Phase 2: Frontend Updates
- Create new room list page at `/rooms`
- Implement startup logic at `/`
- Update routing in App.tsx
- Add room management UI components

### Phase 3: Backend Refactoring
- Remove permanent lock mechanism from promote endpoint
- Update file upload task to use context manager
- Add room duplication logic

### Phase 4: Testing & Migration
- Update all tests for new behavior
- Create migration script for existing rooms
- Deprecation notices for old endpoints

### Phase 5: Cleanup
- Remove deprecated endpoints
- Remove old template locking logic
- Update documentation

## Implementation Checklist

### Backend (Python)
- [x] Design architecture
- [x] Analyze existing lock mechanisms
- [x] Add `room:{room_id}:description` field support
- [x] Add `room:{room_id}:locked` flag support
- [x] Add `room:{room_id}:hidden` flag support
- [x] Add `default_room` key handling
- [ ] Remove `room:{room_id}:template` field from new rooms (kept for backward compatibility)
- [x] Enhance GET /api/rooms endpoint (include description, frameCount, locked, hidden, isDefault)
- [x] Enhance GET /api/rooms/{room_id} endpoint (include description, frameCount, locked, hidden)
- [x] Add PATCH /api/rooms/{room_id} endpoint (update description, locked, hidden)
- [x] Add GET /api/rooms/default endpoint
- [x] Add PUT /api/rooms/default endpoint
- [x] Add POST /api/rooms/{room_id}/duplicate endpoint
- [x] Create `check_room_locked()` helper function
- [x] Add lock checks to all mutating operations (POST /api/rooms/{room_id}/frames, DELETE /api/rooms/{room_id}/frames)
- [x] Remove auto-promotion logic from tasks.py
- [x] Update read_file task to use new default room API
- [ ] Remove permanent trajectory lock from promote endpoint (can be deprecated later)
- [ ] Create migration script to convert permanent locks to locked flags (can be done later)

### Frontend (TypeScript/React)
- [ ] Create RoomListPage component with MUI-X DataGrid
- [ ] Add editable description column in DataGrid
- [ ] Add room list API calls (listRooms, duplicateRoom, setDescription, etc.)
- [ ] Create StartupLogicPage component for `/`
- [ ] Add room management buttons to room view
- [ ] Add lock status indicator to room view
- [ ] Update routing in App.tsx
- [ ] Add room duplication dialog with description input
- [ ] Add default room indicator (star icon) in list
- [ ] Add locked room indicator (lock icon) in list

### Testing
- [ ] Test room description setting/clearing
- [ ] Test room duplication with and without description
- [ ] Test lock/unlock functionality
- [ ] Test that locked rooms reject mutations (append, delete, etc.)
- [ ] Test that locked rooms still allow trajectory lock acquisition (for reads)
- [ ] Test hide/show functionality
- [ ] Test default room setting/unsetting
- [ ] Test startup logic with different room counts
- [ ] Update existing template tests to use new locked flag
- [ ] Test CLI file upload without auto-promotion
- [ ] Test migration script for converting permanent locks
- [ ] Test presenter lock doesn't interfere with room lock

### Documentation
- [ ] Update REST API documentation
- [ ] Update architecture documentation
- [ ] Add migration guide
- [ ] Update CLI documentation