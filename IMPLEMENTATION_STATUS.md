# Implementation Status - Room Management Architecture

## Overview
This document tracks the implementation of the new room management architecture that separates room locking (immutability) from room reusability (duplication).

## ✅ Completed Backend Implementation

### 1. Data Model Changes (routes.py, tasks.py)
All new Redis keys are now supported:
- `room:{room_id}:description` - Optional room description
- `room:{room_id}:locked` - Boolean flag (0/1) for immutability
- `room:{room_id}:hidden` - Boolean flag (0/1) for visibility
- `default_room` - Global key for default room ID

**Note**: `room:{room_id}:template` field is kept for backward compatibility.

### 2. Enhanced Existing Endpoints

#### GET /api/rooms
**Location**: `src/zndraw/app/routes.py` lines 889-938
**Returns**:
```json
[
  {
    "id": "room1",
    "description": "My custom description",  // null if not set
    "frameCount": 42,
    "locked": false,
    "hidden": false,
    "isDefault": true
  }
]
```

#### GET /api/rooms/{room_id}
**Location**: `src/zndraw/app/routes.py` lines 941-982
**Returns**:
```json
{
  "id": "room1",
  "description": "My custom description",  // null if not set
  "frameCount": 42,
  "locked": false,
  "hidden": false
}
```

### 3. New API Endpoints

#### PATCH /api/rooms/{room_id}
**Location**: `src/zndraw/app/routes.py` lines 985-1026
**Purpose**: Update room metadata (description, locked, hidden)
**Request**:
```json
{
  "description": "New description",  // Optional, set to null to clear
  "locked": true,                    // Optional
  "hidden": false                    // Optional
}
```
**Response**: `{"status": "ok"}` (200)

#### GET /api/rooms/default
**Location**: `src/zndraw/app/routes.py` lines 1029-1038
**Purpose**: Get the default room ID
**Response**: `{"roomId": "room1"}` or `{"roomId": null}`

#### PUT /api/rooms/default
**Location**: `src/zndraw/app/routes.py` lines 1041-1070
**Purpose**: Set or unset the default room
**Request**: `{"roomId": "room1"}` or `{"roomId": null}`
**Response**: `{"status": "ok"}` (200)

#### POST /api/rooms/{room_id}/duplicate
**Location**: `src/zndraw/app/routes.py` lines 1073-1156
**Purpose**: Duplicate a room by copying frame mappings and metadata
**Request**:
```json
{
  "newRoomId": "new-room-uuid",      // Optional, auto-generated if not provided
  "description": "Copy of room1"     // Optional
}
```
**Response**:
```json
{
  "status": "ok",
  "roomId": "new-room-uuid",
  "frameCount": 42
}
```

**What gets copied**:
- ✅ Trajectory indices (sorted set) - shares frame data efficiently
- ✅ Geometries hash
- ✅ Bookmarks hash
- ✅ Current frame set to 0
- ✅ Locked/hidden flags initialized to 0
- ✅ Description (from request or empty)
- ❌ Chat messages (not copied)
- ❌ Room settings (not copied)
- ❌ Selections (not copied)

### 4. Lock Enforcement

#### check_room_locked() Helper
**Location**: `src/zndraw/app/routes.py` lines 60-65
**Purpose**: Check if room is locked and return error response if so
```python
def check_room_locked(room_id: str) -> tuple[dict[str, str], int] | None:
    """Check if a room is locked. Returns error response tuple if locked, None if not locked."""
    redis_client = current_app.extensions["redis"]
    locked = redis_client.get(f"room:{room_id}:locked")
    if locked == "1":
        return {"error": "Room is locked and cannot be modified"}, 403
    return None
```

#### Endpoints with Lock Checks
- ✅ `POST /api/rooms/{room_id}/frames` - All actions (append, extend, replace, insert)
- ✅ `DELETE /api/rooms/{room_id}/frames` - Frame deletion

**Response when locked**: `{"error": "Room is locked and cannot be modified"}` (403)

### 5. File Upload Task Changes

#### read_file Task Update
**Location**: `src/zndraw/app/tasks.py` lines 28-118
**Changes**:
- ❌ Removed auto-promotion logic (`POST /api/rooms/{room_id}/promote`)
- ✅ Uses new `PUT /api/rooms/default` API when `make_default=True`
- ✅ No longer sets permanent trajectory locks
- ✅ Room remains unlocked and editable after upload

## 🚧 Remaining Backend Work

### Low Priority (Can Be Done Later)
1. **Deprecate promote endpoint**: The `POST /api/rooms/{room_id}/promote` endpoint can still be used but should eventually be removed
2. **Migration script**: Convert existing permanent trajectory locks to `room:{room_id}:locked` flags
3. **Remove template field**: Eventually phase out `room:{room_id}:template` field completely

### Additional Lock Checks (Optional)
Consider adding lock checks to:
- `PATCH /api/rooms/{room_id}/frames/bulk` - Bulk replace frames
- `POST /api/rooms/{room_id}/geometries` - Geometry updates (or allow for visualization changes)
- `POST /api/rooms/{room_id}/figures` - Figure updates (or allow for visualization changes)
- `DELETE /api/rooms/{room_id}/geometries/{key}` - Geometry deletion
- `DELETE /api/rooms/{room_id}/figures/{key}` - Figure deletion

**Design Decision**: Should locked rooms allow visualization changes (geometries/figures)?
- **Option A**: Lock only trajectory data, allow geometry/figure changes
- **Option B**: Full immutability - lock everything
- **Current Implementation**: Only trajectory data is locked

## 📋 Pending Frontend Implementation

### 1. TypeScript API Client
**File**: `app/src/myapi/rooms.ts` (or similar)
**Needed**:
```typescript
export interface Room {
  id: string;
  description?: string;
  frameCount: number;
  locked: boolean;
  hidden: boolean;
  isDefault?: boolean;
}

export async function listRooms(): Promise<Room[]>
export async function getRoom(roomId: string): Promise<Room>
export async function updateRoom(roomId: string, updates: Partial<Pick<Room, 'description' | 'locked' | 'hidden'>>): Promise<void>
export async function duplicateRoom(roomId: string, newRoomId?: string, description?: string): Promise<{roomId: string, frameCount: number}>
export async function getDefaultRoom(): Promise<string | null>
export async function setDefaultRoom(roomId: string | null): Promise<void>
```

### 2. Room List Page
**Component**: `RoomListPage.tsx`
**Features**:
- MUI-X DataGrid showing all rooms
- Columns: Name/ID, Description (editable), Frame Count, Actions (Lock, Hide, Set Default, Duplicate, Open)
- Icons: Lock icon for locked rooms, Star icon for default room, Eye-slash for hidden
- Inline editing for descriptions
- Action buttons/menu per row

### 3. Startup Logic
**Component**: Update routing in `App.tsx`
**Logic**:
```typescript
// Root path "/" logic
async function handleRootPath() {
  const rooms = await listRooms();
  
  if (rooms.length === 0) {
    // No rooms - create empty template
    const newRoomId = generateUUID();
    navigate(`/rooms/${newRoomId}/${userUUID}?template=empty`);
  } else if (rooms.length === 1) {
    // One room - navigate to it
    navigate(`/rooms/${rooms[0].id}/${userUUID}`);
  } else {
    // Multiple rooms
    const defaultRoomId = await getDefaultRoom();
    if (defaultRoomId) {
      navigate(`/rooms/${defaultRoomId}/${userUUID}`);
    } else {
      navigate('/rooms'); // Show room list
    }
  }
}
```

### 4. Room Management UI
**Add to room view**:
- Lock/Unlock button (only for owners/admins)
- "Set as Default" button
- "Duplicate Room" button with dialog
- Description editor (text field in header)
- Lock indicator badge when viewing locked room

## 🧪 Testing Strategy

### Backend Tests (pytest)
Create `tests/test_room_management.py`:
```python
def test_list_rooms_includes_metadata()
def test_get_room_details()
def test_update_room_description()
def test_update_room_locked_flag()
def test_update_room_hidden_flag()
def test_locked_room_rejects_append()
def test_locked_room_rejects_delete()
def test_locked_room_allows_reads()
def test_set_default_room()
def test_get_default_room()
def test_unset_default_room()
def test_duplicate_room_copies_frames()
def test_duplicate_room_copies_geometries()
def test_duplicate_room_copies_bookmarks()
def test_duplicate_room_does_not_copy_chat()
def test_duplicate_room_with_description()
def test_duplicate_room_with_auto_id()
def test_duplicate_nonexistent_room_fails()
def test_duplicate_to_existing_room_fails()
```

### Manual Testing Checklist
- [ ] Verify GET /api/rooms returns all metadata fields
- [ ] Verify PATCH updates work individually and in combination
- [ ] Verify locked rooms reject POST /api/rooms/{id}/frames
- [ ] Verify locked rooms reject DELETE /api/rooms/{id}/frames
- [ ] Verify default room can be set and retrieved
- [ ] Verify room duplication creates proper copies
- [ ] Verify duplicated room shares frame data (check Redis keys)
- [ ] Upload file with CLI and verify it's not auto-promoted
- [ ] Upload file with make_default=True and verify default_room is set

### Frontend Tests
- [ ] Test room list renders correctly
- [ ] Test description editing works
- [ ] Test lock/unlock buttons work
- [ ] Test duplicate dialog works
- [ ] Test default room indicator updates
- [ ] Test startup logic with 0, 1, and N rooms
- [ ] Test navigation to default room works

## 📝 Code Quality Notes

### Lint Warnings
The following pre-existing lint warnings appear in routes.py (not introduced by these changes):
- Line 1831: Unused variable `user_extensions_key`
- Line 1504: Type error in celery job worker call
- Line 1932, 1977: Type errors in job ID manipulation
- Line 1857: Route decorator type mismatch

These should be addressed in a separate cleanup PR.

### Design Decisions Made

1. **Backward Compatibility**: Kept `room:{room_id}:template` field to avoid breaking existing code
2. **Lock Scope**: Only trajectory mutations are checked for locks (geometries/figures still allowed)
3. **Default Geometries**: Duplicated rooms get default geometries if source had none
4. **Frame Data Sharing**: Duplication copies references (sorted set), not actual frame data
5. **Description Optional**: All rooms can have descriptions, but they're optional
6. **HTTP Status Codes**: 403 for locked rooms, 404 for not found, 409 for conflicts

## 🚀 Deployment Notes

### Redis Key Migration (Optional)
If you want to convert existing permanent trajectory locks to locked flags:
```python
# Pseudocode for migration script
for room_id in redis.scan_iter("room:*:template"):
    lock_key = f"room:{room_id}:lock:trajectory"
    if redis.ttl(lock_key) == -1:  # Permanent lock
        redis.set(f"room:{room_id}:locked", "1")
        redis.delete(lock_key)
```

### API Documentation Updates
Update `rest-endpoints.md` with:
- New endpoint specifications
- Updated request/response examples
- Lock error response documentation
- Room metadata field descriptions

### Rollout Strategy
1. Deploy backend changes (backward compatible)
2. Test with existing clients (should continue working)
3. Deploy frontend changes
4. Run migration script (if needed)
5. Update documentation
6. Deprecate old template endpoints

## 📚 Related Documentation
- `startup-and-template-logic.md` - Detailed architecture design
- `rest-endpoints.md` - API endpoint documentation
- `index-mapping.md` - Frame data storage explanation
- `communication.md` - WebSocket communication patterns
