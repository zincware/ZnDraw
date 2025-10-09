# Room Metadata Implementation Plan

## Overview
Add comprehensive metadata tracking for rooms when files are uploaded. Metadata includes file information, upload method, and original frame count. This enables:
- Tracking which files are loaded in which rooms
- Creating new rooms from the same physical storage without re-uploading
- Searching rooms by their source files
- Enhanced file browser UI showing which files are already in use


## Requirements Summary

### Core Features
1. **Metadata Storage**: Store metadata in Redis as `room:{roomId}:metadata` hash
2. **Python Client API**: Add `vis.metadata` property as MutableMapping for easy access
3. **REST Endpoints**: Create GET/POST/DELETE endpoints for metadata management
4. **Lock Protection**: Metadata writes must respect the room lock
5. **File Browser Integration**: Show file usage status in file browser
6. **Room Duplication**: Allow creating new rooms from existing file metadata
7. **Search Capability**: Enable regex search across file metadata
8. **Search Capability**: Enable regex search across file browser / filenames in filesystem browser

### Metadata Schema
```json
{
  "relative_file_path": "data/file.xyz",
  "root_path": "/path/to/root",
  "upload_method": "file_browser",
  "frame_count_original": 100,
  "file_size": 12345678,
  "file_modified_timestamp": "2025-10-08T12:34:56.789Z",
  "upload_timestamp": "2025-10-08T12:35:00.123Z"
}
```

**New Fields for File Change Detection**:
- `file_size`: Size of the original file in bytes
- `file_modified_timestamp`: Last modification time of the file (ISO 8601 format)
- `upload_timestamp`: When the file was uploaded to the room

**Purpose**: Enable detection of file changes since upload. If current file size or modification time differs from stored metadata, the file has been modified and should be re-uploaded instead of reused.

## Lock Handling Philosophy

**IMPORTANT**: Metadata mutation protection follows the same pattern as other room mutations:

### Backend Lock Checking (Server-Side)
- All metadata mutation endpoints (`POST`, `DELETE`) MUST call `check_room_locked(room_id)` FIRST
- If the room is locked (either permanently or via metadata lock), the endpoint returns an error
- Lock checking happens at the HTTP endpoint level, not in the manager class

### Client-Side Behavior
- Clients DO NOT need to wrap metadata writes in `with vis.lock:`
- The backend endpoint will reject the request if the room is locked
- This is consistent with how other mutation operations work (append, delete, etc.)

### Example Flow
```python
# Client code (NO lock needed):
vis.metadata["key"] = "value"  # Sends POST request to backend

# Backend endpoint (routes.py):
@main.route("/api/rooms/<room_id>/metadata", methods=["POST"])
def update_room_metadata(room_id: str):
    # FIRST: Check if room is locked
    lock_error = check_room_locked(room_id)
    if lock_error:
        return lock_error  # Returns 423 Locked or similar
    
    # THEN: Perform mutation
    manager = RoomMetadataManager(redis_client, room_id)
    manager.update(data)
```

**Why**: 
- Consistent with existing mutation endpoints (delete_frames, append_frame, etc.)
- Centralized lock checking in one place (backend)
- Simpler client code
- Server enforces invariants, not client

## File Change Detection Strategy

### Purpose
**ONLY**: When loading a file through the file browser, check if that exact file (same size + modification time) is already loaded in an existing room. If yes, offer to open a new room loading from zarr file storage instead of re-uploading.

### Simple Logic
When user clicks "Load File" in file browser:
1. Check all existing rooms for metadata matching this file path
2. Compare current file's size and mtime with stored metadata
3. If **exact match found**: Offer to create a new room loading from zarr file storage
4. If **no match OR file changed**: Upload to new room as normal

### Storage on Upload
```python
file_stat = Path(file).stat()
vis.metadata["file_size"] = str(file_stat.st_size)
vis.metadata["file_modified_timestamp"] = datetime.fromtimestamp(
    file_stat.st_mtime, tz=timezone.utc
).isoformat()
```

### File Browser UI Flow

#### Case 1: First Time Loading File
```
User clicks "Load my_very_long_file_name.xyz"
  ↓
Backend checks: Does any room have this exact file?
  ↓
NO → Generate short room name: "my_very_long_file_n"
   → If name exists, add hash: "my_very_long_f_a3c4d2"
   → Proceed with upload
```

#### Case 2: File Already Loaded
```
User clicks "Load data.xyz" (already loaded in "experiment_1")
  ↓
Backend returns: {status: "file_already_loaded", existingRoom: "experiment_1"}
  ↓
Frontend shows dialog:
  "This file is already loaded in room 'experiment_1'"
  
  [Open Existing Room] → Navigate to experiment_1
  
  [Create New Room]    → POST to /create-room-from-file
                         → Creates "data_xyz_abc123" (short, unique name)
                         → Reuses physical storage (fast!)
                         → Opens new room
  
  [Upload Anyway]      → POST to /load with force_upload: true
                         → Ignores existing match
                         → Creates new room with full upload
```

#### Case 3: File Changed Since Upload
```
User clicks "Load data.xyz" (exists but file modified)
  ↓
Backend: No exact match (size/mtime different)
  ↓
Generate new room name and proceed with upload
```

**Key Features**:
1. **Short room names**: "my_very_long_file_name.xyz" → "my_very_long_file_n"
2. **Collision handling**: Automatically adds hash if name taken
3. **Fast room creation**: Reuse storage with identity mapping (no re-upload!)
4. **Force upload option**: User can override duplicate detection
5. **Independent rooms**: Each room has its own trajectory but shares storage

### Implementation: Two Simple Endpoints

**1. Check if file exists** (in `load_file`):
```python
existing_room = _find_room_with_exact_file(target_path, redis_client)
if existing_room:
    return {"status": "file_already_loaded", "existingRoom": existing_room}
```

**2. Create new room from existing** (`create-room-from-file`):
```python
# Copy Zarr storage reference
# Create identity mapping: {0: 0, 1: 1, 2: 2, ...}
# New room shares physical storage, no re-upload!
```

### That's It!
- **Check on upload**: Is this exact file already here?
- **If yes**: Offer to create new room reusing storage (fast) OR open existing
- **If no/changed**: Normal upload
- **No UI indicators** about changes in directory listings
- **No automatic updates** to existing rooms
- Simple, focused, solves the actual problem: fast room creation without re-upload

## Implementation Plan

### Phase 1: Backend Storage & API (Core Infrastructure)

#### 1.1 Create Metadata Manager Class
**File**: `src/zndraw/app/metadata_manager.py` (NEW)

**Purpose**: Encapsulate metadata operations following separation of concerns

**Implementation**:
```python
class RoomMetadataManager:
    """Manages room metadata in Redis.
    
    Provides CRUD operations for room metadata with proper validation
    and atomic operations.
    """
    
    def __init__(self, redis_client, room_id: str):
        self.redis = redis_client
        self.room_id = room_id
        self.key = f"room:{room_id}:metadata"
    
    def get_all(self) -> dict[str, str]:
        """Get all metadata for the room."""
        
    def get(self, field: str) -> str | None:
        """Get specific metadata field."""
        
    def set(self, field: str, value: str) -> None:
        """Set a metadata field."""
        
    def update(self, data: dict[str, str]) -> None:
        """Update multiple metadata fields atomically."""
        
    def delete(self, field: str) -> None:
        """Delete a metadata field."""
        
    def clear(self) -> None:
        """Clear all metadata for the room."""
        
    def exists(self) -> bool:
        """Check if metadata exists for the room."""
```

**Why**: 
- Single responsibility principle
- Reusable across different parts of the application
- Easier to test in isolation
- Consistent Redis key naming

#### 1.2 Add REST Endpoints
**File**: `src/zndraw/app/routes.py` (UPDATE)

**Endpoints to Add**:

1. **GET `/api/rooms/<room_id>/metadata`**
   - Purpose: Retrieve room metadata
   - Response: `{"metadata": {"relative_file_path": "...", ...}}`
   - Error handling: 404 if room doesn't exist, 200 with empty dict if no metadata

2. **POST `/api/rooms/<room_id>/metadata`**
   - Purpose: Update room metadata (partial or full)
   - Body: `{"relative_file_path": "data/file.xyz", ...}`
   - **Lock check**: Call `check_room_locked(room_id)` before writing - returns error if locked
   - Validation: Ensure all values are strings
   - Response: `{"success": true, "metadata": {...}}`
   - **Important**: Lock check happens in endpoint, NOT in client code

3. **DELETE `/api/rooms/<room_id>/metadata/<field>`**
   - Purpose: Delete specific metadata field
   - **Lock check**: Call `check_room_locked(room_id)` before deletion
   - Response: `{"success": true}`

**Implementation Details**:
```python
@main.route("/api/rooms/<string:room_id>/metadata", methods=["GET"])
def get_room_metadata(room_id: str):
    """Get all metadata for a room."""
    # Use RoomMetadataManager
    # Return empty dict if no metadata exists
    
@main.route("/api/rooms/<string:room_id>/metadata", methods=["POST"])
def update_room_metadata(room_id: str):
    """Update room metadata. Respects room lock."""
    # CRITICAL: Check lock FIRST using check_room_locked(room_id)
    # If locked, return error response immediately
    lock_error = check_room_locked(room_id)
    if lock_error:
        return lock_error
    
    # Validate input (all values must be strings)
    # Use RoomMetadataManager.update()
    # Return updated metadata

@main.route("/api/rooms/<string:room_id>/metadata/<string:field>", methods=["DELETE"])
def delete_room_metadata_field(room_id: str, field: str):
    """Delete specific metadata field. Respects room lock."""
    # CRITICAL: Check lock FIRST
    lock_error = check_room_locked(room_id)
    if lock_error:
        return lock_error
    
    # Use RoomMetadataManager.delete(field)
    # Return success response
```

**Why**:
- Follows existing REST API patterns in routes.py
- Lock checking consistent with other mutation endpoints
- Proper error handling and validation

#### 1.3 Update Existing Endpoints
**File**: `src/zndraw/app/routes.py` (UPDATE)

**Changes**:

1. **`GET /api/rooms` (list_rooms function)**
   - Add metadata to room objects in response
   - Include metadata in the returned room data structure
   ```python
   for room_id in sorted(room_ids):
       # ... existing code ...
       metadata_manager = RoomMetadataManager(redis_client, room_id)
       metadata = metadata_manager.get_all()
       # Add to room object: "metadata": metadata
   ```

2. **`GET /api/rooms/<room_id>` (get_room function)**
   - Include metadata in the response
   ```python
   metadata_manager = RoomMetadataManager(redis_client, room_id)
   metadata = metadata_manager.get_all()
   # Add to response: "metadata": metadata
   ```

**Why**: 
- Makes metadata available to frontend without extra requests
- Consistent with RESTful design
- Minimal overhead (Redis HGETALL is O(N) where N is number of fields)

### Phase 2: File Upload Integration

#### 2.1 Update Celery Task
**File**: `src/zndraw/app/tasks.py` (UPDATE)

**Changes in `read_file` task**:
```python
@shared_task
def read_file(
    file: str,
    room: str,
    server_url: str = "http://localhost:5000",
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    make_default: bool = False,
    batch_size: int = 10,
    root_path: str | None = None,  # NEW parameter
) -> None:
    from datetime import datetime, timezone
    
    file_path = Path(file)
    vis = ZnDraw(room=room, url=server_url, user="uploader", description=f"{file}")
    
    # ... existing validation and upload code ...
    
    # AFTER successful upload, set metadata via API
    # Note: We use vis.metadata directly - the backend endpoint will check the lock
    file_stat = file_path.stat()
    
    vis.metadata["relative_file_path"] = str(file_path.relative_to(root_path)) if root_path else str(file_path)
    vis.metadata["root_path"] = root_path or str(file_path.parent)
    vis.metadata["upload_method"] = "file_browser"
    vis.metadata["frame_count_original"] = str(len(vis))
    vis.metadata["file_size"] = str(file_stat.st_size)
    vis.metadata["file_modified_timestamp"] = datetime.fromtimestamp(file_stat.st_mtime, tz=timezone.utc).isoformat()
    vis.metadata["upload_timestamp"] = datetime.now(timezone.utc).isoformat()
```

**Why**:
- Metadata set AFTER successful upload (as required)
- **No `with vis.lock:` needed** - lock checking happens in the backend REST endpoint
- Captures file size and modification time for change detection
- All values converted to strings (Redis hash requirement)
- Uses timezone-aware UTC timestamps

#### 2.2 Update File Browser Load Endpoint
**File**: `src/zndraw/app/file_browser.py` (UPDATE)

**Changes in `load_file` function**:
```python
@file_browser.route("/load", methods=["POST"])
def load_file():
    # ... existing validation code ...
    
    # Get file browser root for relative path calculation
    root = current_app.config.get("FILE_BROWSER_ROOT", ".")
    
    task = read_file.delay(
        file=str(target_path),
        room=room,
        server_url=server_url,
        start=start,
        stop=stop,
        step=step,
        make_default=make_default,
        root_path=root,  # NEW parameter
    )
```

**Why**:
- Passes root_path to task for proper relative path calculation
- Minimal change to existing code

### Phase 3: Python Client API

#### 3.1 Create Metadata Manager for Client
**File**: `src/zndraw/metadata_manager.py` (NEW)

**Purpose**: Provide MutableMapping interface for client-side metadata access

**Implementation**:
```python
from collections.abc import MutableMapping
import typing as t

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class RoomMetadata(MutableMapping):
    """Client-side interface for room metadata.
    
    Provides dict-like access to room metadata with automatic
    synchronization to the server.
    
    Examples:
        vis.metadata["file"] = "data.xyz"
        file = vis.metadata["file"]
        del vis.metadata["file"]
        for key in vis.metadata:
            print(key, vis.metadata[key])
    """
    
    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance
        self._cache: dict[str, str] = {}
        self._loaded = False
    
    def _ensure_loaded(self) -> None:
        """Lazy load metadata from server."""
        if not self._loaded:
            self._cache = self.vis.api.get_metadata()
            self._loaded = True
    
    def __getitem__(self, key: str) -> str:
        self._ensure_loaded()
        return self._cache[key]
    
    def __setitem__(self, key: str, value: str) -> None:
        """Set metadata field. Requires connection."""
        if not self.vis.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        self.vis.api.set_metadata({key: value})
        self._cache[key] = value
    
    def __delitem__(self, key: str) -> None:
        """Delete metadata field. Requires connection."""
        if not self.vis.socket.sio.connected:
            raise RuntimeError("Client is not connected.")
        self.vis.api.delete_metadata_field(key)
        if key in self._cache:
            del self._cache[key]
    
    def __iter__(self):
        self._ensure_loaded()
        return iter(self._cache)
    
    def __len__(self) -> int:
        self._ensure_loaded()
        return len(self._cache)
    
    def __repr__(self) -> str:
        self._ensure_loaded()
        return f"RoomMetadata({self._cache!r})"
```

**Why**:
- Follows same pattern as Geometries and Figures managers
- Lazy loading for performance
- Dict-like interface is intuitive for users
- Local cache reduces network calls

#### 3.2 Update APIManager
**File**: `src/zndraw/api_manager.py` (UPDATE)

**Add methods**:
```python
def get_metadata(self) -> dict[str, str]:
    """Get all metadata for the room."""
    response = requests.get(f"{self.url}/api/rooms/{self.room}/metadata")
    response.raise_for_status()
    return response.json().get("metadata", {})

def set_metadata(self, data: dict[str, str]) -> None:
    """Update room metadata."""
    response = requests.post(
        f"{self.url}/api/rooms/{self.room}/metadata",
        json=data
    )
    response.raise_for_status()

def delete_metadata_field(self, field: str) -> None:
    """Delete a metadata field."""
    response = requests.delete(
        f"{self.url}/api/rooms/{self.room}/metadata/{field}"
    )
    response.raise_for_status()
```

**Why**:
- Centralizes HTTP logic in APIManager (existing pattern)
- Provides clean interface for RoomMetadata to use

#### 3.3 Update ZnDraw Class
**File**: `src/zndraw/zndraw.py` (UPDATE)

**Changes**:
1. Add import: `from zndraw.metadata_manager import RoomMetadata`
2. Add field: `_metadata: RoomMetadata | None = dataclasses.field(default=None, init=False)`
3. Add property:
```python
@property
def metadata(self) -> RoomMetadata:
    """Access room metadata as a dict-like object."""
    if self._metadata is None:
        self._metadata = RoomMetadata(self)
    return self._metadata
```

**Why**:
- Lazy initialization (only created when accessed)
- Consistent with geometries and figures properties
- Type hints ensure proper usage

### Phase 4: File Browser Enhancement

#### 4.0 Add Room Name Generation Utility
**File**: `src/zndraw/utils.py` (UPDATE)

**Move and enhance `path_to_room` function**:
```python
import hashlib
import re

def generate_room_name(
    file_path: str,
    redis_client=None,
    max_length: int = 20,
    hash_length: int = 6
) -> str:
    """Generate a short, unique room name from a file path.
    
    Rules:
    1. Extract filename without extension
    2. Replace non-alphanumeric with underscore
    3. Truncate to max_length characters
    4. If room name exists in Redis, append hash of full path
    
    Parameters
    ----------
    file_path : str
        Full or relative path to file
    redis_client : Redis, optional
        Redis client to check for existing rooms
    max_length : int
        Maximum length of base name before hash
    hash_length : int
        Length of hash suffix if name collision occurs
    
    Returns
    -------
    str
        Unique room name
    
    Examples
    --------
    >>> generate_room_name("data/my_very_long_file_name.xyz")
    'my_very_long_file_n'
    
    >>> generate_room_name("data/my_very_long_file_name.xyz", redis_client)
    'my_very_long_f_a3c4d2'  # if collision detected
    """
    from pathlib import Path
    
    # Extract filename without extension
    filename = Path(file_path).stem
    
    # Replace non-alphanumeric with underscore
    clean_name = re.sub(r"[^a-zA-Z0-9\-]", "_", filename)
    
    # Truncate to max_length
    base_name = clean_name[:max_length]
    
    # Check if room exists (if redis_client provided)
    if redis_client is None:
        return base_name
    
    # Check for collision
    room_exists = redis_client.exists(f"room:{base_name}:trajectory:indices")
    
    if not room_exists:
        return base_name
    
    # Collision detected - add hash suffix
    # Hash the full path for uniqueness
    path_hash = hashlib.md5(file_path.encode()).hexdigest()[:hash_length]
    
    # Make room for hash by truncating base name
    truncated_base = clean_name[:max_length - hash_length - 1]
    
    return f"{truncated_base}_{path_hash}"


def path_to_room(path: str) -> str:
    """Legacy function for backward compatibility.
    
    Converts path to room name without truncation or collision checking.
    Use generate_room_name() for new code.
    """
    return re.sub(r"[^a-zA-Z0-9\-]", "_", path)
```

**Update cli.py to use new function**:
```python
from zndraw.utils import generate_room_name, path_to_room

# In CLI code where rooms are created:
# OLD: room = path_to_room(p)
# NEW: room = generate_room_name(p, redis_client=redis_client)
```

**Why**:
- **Short names**: Truncates long filenames to manageable length
- **Unique**: Adds hash suffix if name collision occurs
- **Consistent**: Same logic used in CLI and file browser
- **Centralized**: Single source of truth in utils module

#### 4.1 Update File Load Endpoint with Duplicate Detection
**File**: `src/zndraw/app/file_browser.py` (UPDATE)

**Update load_file function**:
```python
@file_browser.route("/load", methods=["POST"])
def load_file():
    """Load a file into ZnDraw, checking for existing room with same file first.
    
    Request Body:
        path: str - Relative path to file
        room: str, optional - Custom room name
        start/stop/step: int, optional - Frame range
        make_default: bool, optional
        force_upload: bool, optional - If true, skip duplicate check and upload anyway
    """
    # ... existing validation code ...
    
    data = request.get_json()
    force_upload = data.get("force_upload", False)
    
    redis_client = current_app.extensions["redis"]
    
    # Check if this exact file is already loaded (unless force_upload is true)
    if not force_upload:
        existing_room = _find_room_with_exact_file(target_path, redis_client)
        
        if existing_room:
            return jsonify({
                "status": "file_already_loaded",
                "existingRoom": existing_room,
                "message": f"This file is already loaded in room '{existing_room}'",
                "filePath": str(target_path),
                "options": {
                    "openExisting": f"Open room '{existing_room}'",
                    "createNew": "Create new room (reuse storage)",
                    "forceUpload": "Upload anyway (ignore existing)"
                }
            }), 200
    
    # Get or generate room name using new utility
    room = data.get("room")
    if not room:
        from zndraw.utils import generate_room_name
        room = generate_room_name(str(target_path), redis_client=redis_client)
    
    # File not loaded, has changed, or force_upload=true - proceed with upload
    # ... existing upload code ...
```

**Helper function**:
```python
def _find_room_with_exact_file(file_path: Path, redis_client) -> str | None:
    """Find a room that has this exact file (matching path, size, and mtime).
    
    Returns:
        Room ID if exact match found, None otherwise.
    """
    from datetime import datetime, timezone
    from zndraw.app.metadata_manager import RoomMetadataManager
    
    # Get current file stats
    current_stat = file_path.stat()
    current_size = str(current_stat.st_size)
    current_mtime = datetime.fromtimestamp(
        current_stat.st_mtime, tz=timezone.utc
    ).isoformat()
    
    # Scan all rooms for matching metadata
    for key in redis_client.scan_iter(match="room:*:metadata"):
        room_id = key.split(":")[1]
        manager = RoomMetadataManager(redis_client, room_id)
        metadata = manager.get_all()
        
        # Check exact match: same path, size, and mtime
        if (metadata.get("relative_file_path") == str(file_path) and
            metadata.get("file_size") == current_size and
            metadata.get("file_modified_timestamp") == current_mtime):
            return room_id
    
    return None
```

**Why**:
- Check before upload: "Is this exact file already here?"
- Three options if match found:
  1. **Open existing room**
  2. **Create new room** (reuse storage via create-room-from-file)
  3. **Force upload** (set `force_upload: true` to ignore existing)
- Uses new `generate_room_name()` for short, unique room names
- If no match or file changed, proceed with normal upload

#### 4.2 Add Room Creation from Existing File Endpoint
**File**: `src/zndraw/app/file_browser.py` (UPDATE)

**New endpoint**:
```python
@file_browser.route("/create-room-from-file", methods=["POST"])
def create_room_from_existing_file():
    """Create a new room reusing physical storage from an existing room.
    
    Body:
        {
            "sourceRoom": "room1",
            "newRoom": "room2",  # optional, will generate if not provided
            "description": "Copy of room1"
        }
    
    Creates a new room with identity mapping {i: i for i in range(frame_count_original)}
    WITHOUT re-uploading the file. Just reuses the same physical Zarr storage.
    """
    data = request.get_json()
    source_room = data.get("sourceRoom")
    new_room = data.get("newRoom")
    description = data.get("description")
    
    if not source_room:
        return jsonify({"error": "sourceRoom is required"}), 400
    
    redis_client = current_app.extensions["redis"]
    
    # Get source room metadata
    from zndraw.app.metadata_manager import RoomMetadataManager
    source_metadata_manager = RoomMetadataManager(redis_client, source_room)
    source_metadata = source_metadata_manager.get_all()
    
    if not source_metadata:
        return jsonify({"error": "Source room has no metadata"}), 400
    
    frame_count = source_metadata.get("frame_count_original")
    if not frame_count:
        return jsonify({"error": "Source room missing frame_count_original"}), 400
    
    # Generate new room name if not provided
    if not new_room:
        from zndraw.utils import generate_room_name
        # Use the original file path from metadata to generate consistent name
        file_path = source_metadata.get("relative_file_path", source_room)
        new_room = generate_room_name(file_path, redis_client=redis_client)
    
    # Use existing duplicate_room logic to copy the room with identity mapping
    # This creates: room:new_room:trajectory:indices as sorted set with {i: i} mapping
    from zndraw.app.routes import get_storage
    
    source_storage = get_storage(source_room)
    new_storage = get_storage(new_room)
    
    # Copy all frames with identity mapping
    frame_count_int = int(frame_count)
    indices_key = f"room:{new_room}:trajectory:indices"
    
    # Create identity mapping: frame 0 -> 0, frame 1 -> 1, etc.
    for i in range(frame_count_int):
        physical_key = f"frame_{i}"
        redis_client.zadd(indices_key, {physical_key: i})
        
        # Copy physical frame data from source to new storage
        if i in source_storage:
            new_storage[i] = source_storage[i]
    
    # Copy metadata to new room
    new_metadata_manager = RoomMetadataManager(redis_client, new_room)
    new_metadata_manager.update(source_metadata)
    
    # Set room description and properties
    if description:
        redis_client.set(f"room:{new_room}:description", description)
    redis_client.set(f"room:{new_room}:locked", "0")
    redis_client.set(f"room:{new_room}:hidden", "0")
    
    return jsonify({
        "status": "success",
        "roomId": new_room,
        "sourceRoom": source_room,
        "frameCount": frame_count_int,
        "message": f"Room '{new_room}' created from '{source_room}' without re-uploading"
    }), 201
```

**Why**:
- **Fast room creation** - no file re-upload needed
- **Reuses physical storage** - just creates new logical mapping
- **Identity mapping** - {0: 0, 1: 1, 2: 2, ...} points to original frames
- **Preserves metadata** - new room inherits file info from source

### Phase 5: Search Functionality

#### 5.1 Add Search to Room List
**File**: `src/zndraw/app/routes.py` (UPDATE)

**Update list_rooms endpoint**:
```python
@main.route("/api/rooms", methods=["GET"])
def list_rooms():
    """List all active rooms with metadata.
    
    Query params:
        search: Optional regex pattern to search in metadata
    """
    redis_client = current_app.extensions["redis"]
    search_pattern = request.args.get("search")
    
    # ... existing room collection code ...
    
    rooms = []
    for room_id in sorted(room_ids):
        # ... existing code to get room info ...
        
        metadata_manager = RoomMetadataManager(redis_client, room_id)
        metadata = metadata_manager.get_all()
        
        # Filter by search pattern if provided
        if search_pattern:
            import re
            try:
                pattern = re.compile(search_pattern, re.IGNORECASE)
                # Search in metadata values
                if not any(pattern.search(str(v)) for v in metadata.values()):
                    continue
            except re.error:
                # Invalid regex, skip filtering
                pass
        
        room_data = {
            # ... existing fields ...
            "metadata": metadata
        }
        rooms.append(room_data)
    
    return jsonify(rooms)
```

**Why**:
- Enables finding rooms by source file
- Regex provides flexible search
- Server-side filtering reduces client load

#### 5.2 Add Search to File Browser
**File**: `src/zndraw/app/file_browser.py` (UPDATE)

**Update list_directory endpoint**:
```python
@file_browser.route("/list", methods=["GET"])
def list_directory():
    """List contents of a directory.
    
    Query params:
        path: Relative path from browser root
        search: Optional regex pattern for file names
    """
    # ... existing code ...
    
    search_pattern = request.args.get("search")
    
    items = []
    for item in sorted(target_path.iterdir(), ...):
        # ... existing filtering ...
        
        # Apply search filter
        if search_pattern:
            import re
            try:
                pattern = re.compile(search_pattern, re.IGNORECASE)
                if not pattern.search(item.name):
                    continue
            except re.error:
                pass  # Invalid regex, include all
        
        # ... rest of existing code ...
```

**Why**:
- Helps users find files in large directories
- Consistent with room list search
- Regex allows complex patterns

### Phase 6: Testing

#### 6.1 Create Metadata Manager Tests
**File**: `tests/test_metadata.py` (NEW)

**Test cases**:
```python
def test_metadata_manager_crud(redis_client):
    """Test basic CRUD operations."""
    
def test_metadata_get_nonexistent(redis_client):
    """Test getting metadata for room that doesn't exist."""
    
def test_metadata_update_atomic(redis_client):
    """Test atomic updates."""
    
def test_metadata_clear(redis_client):
    """Test clearing all metadata."""
```

#### 6.2 Create REST Endpoint Tests
**File**: `tests/test_metadata.py` (UPDATE)

**Test cases**:
```python
def test_get_room_metadata_empty(server):
    """Test GET when no metadata exists."""
    
def test_post_room_metadata(server):
    """Test creating/updating metadata."""
    
def test_metadata_respects_lock(server):
    """Test that metadata writes check room lock."""
    
def test_delete_metadata_field(server):
    """Test deleting specific field."""
```

#### 6.3 Create Client API Tests
**File**: `tests/test_metadata.py` (UPDATE)

**Test cases**:
```python
def test_client_metadata_get_set(server):
    """Test vis.metadata dict-like interface."""
    
def test_client_metadata_requires_connection(server):
    """Test that writes require active connection."""
    
def test_client_metadata_lazy_load(server):
    """Test lazy loading behavior."""
```

#### 6.4 Create Integration Tests
**File**: `tests/test_metadata.py` (UPDATE)

**Test cases**:
```python
def test_read_file_sets_metadata(server, tmp_path):
    """Test that read_file task sets metadata after upload."""
    
def test_room_list_includes_metadata(server):
    """Test that room list includes metadata."""

def test_file_already_loaded_detection(server, tmp_path):
    """Test that loading same file twice detects existing room."""
    # Upload file once, try to upload again, should get existing room response
    
def test_file_changed_allows_new_upload(server, tmp_path):
    """Test that modified file can be uploaded to new room."""
    # Upload file, modify it (change size or mtime), upload again should succeed

def test_metadata_includes_file_stats(server, tmp_path):
    """Test that metadata includes file_size and file_modified_timestamp."""

def test_force_upload_bypasses_duplicate_check(server, tmp_path):
    """Test that force_upload=true allows uploading duplicate file."""
    # Upload file, try to upload again with force_upload=true, should succeed

def test_create_room_from_file_reuses_storage(server, tmp_path):
    """Test creating new room from existing file shares physical storage."""
```

**File**: `tests/test_room_naming.py` (NEW)

**Test cases**:
```python
def test_generate_room_name_short_filename():
    """Test that short filenames are used as-is."""
    assert generate_room_name("data.xyz") == "data"

def test_generate_room_name_truncates_long_filename():
    """Test that long filenames are truncated to max_length."""
    long_name = "my_very_long_file_name_that_exceeds_twenty_chars.xyz"
    result = generate_room_name(long_name, max_length=20)
    assert len(result) <= 20
    assert result == "my_very_long_file_n"

def test_generate_room_name_collision_adds_hash(redis_client):
    """Test that collision detection adds hash suffix."""
    # Create existing room
    redis_client.set("room:test_file:trajectory:indices", "")
    
    # Generate name for file that would create same room name
    result = generate_room_name("test/file.xyz", redis_client, max_length=20)
    assert result.startswith("test_file_")
    assert len(result.split("_")[-1]) == 6  # hash length

def test_generate_room_name_replaces_special_chars():
    """Test that special characters are replaced with underscores."""
    assert generate_room_name("data/file-name.xyz") == "data_file_name"

@pytest.mark.parametrize("filename,expected_base", [
    ("a.xyz", "a"),
    ("my_file.pdb", "my_file"),
    ("data/structure.h5", "data_structure"),
])
def test_generate_room_name_various_inputs(filename, expected_base):
    """Test room name generation with various inputs."""
    result = generate_room_name(filename)
    assert result.startswith(expected_base[:20])
```

**Why**:
- Comprehensive tests for room naming logic
- Tests truncation, collision detection, and special character handling
- Uses parametrize to avoid code duplication (AGENTS.md)
- Focused tests for actual use cases
- Tests follow existing patterns in conftest.py
- Each test is specific and focused (AGENTS.md requirement)

### Phase 7: Documentation

#### 7.1 Update API Documentation
**Files**: 
- `rest-endpoints.md` - Document new metadata endpoints
- `README.md` - Add metadata usage examples

#### 7.2 Add Docstrings
- All new functions must have numpy-style docstrings (AGENTS.md requirement)
- Include type hints throughout
- Keep docstrings concise

## Summary of Files to Create/Update

### New Files
1. `src/zndraw/app/metadata_manager.py` - Backend metadata manager with Redis operations
2. `src/zndraw/metadata_manager.py` - Client-side metadata interface (MutableMapping)
3. `tests/test_metadata.py` - Comprehensive test suite including duplicate detection tests
4. `tests/test_room_naming.py` - Tests for room name generation and collision handling

### Files to Update
1. `src/zndraw/utils.py` - Add `generate_room_name()` function with truncation and collision detection
2. `src/zndraw/cli.py` - Update to use `generate_room_name()` instead of `path_to_room()`
3. `src/zndraw/app/routes.py` - Add 3 REST endpoints (GET/POST/DELETE metadata), update room list/get to include metadata, add lock checking
4. `src/zndraw/app/tasks.py` - Update read_file to capture and store file size/mtime metadata after upload
5. `src/zndraw/app/file_browser.py` - Add duplicate detection, force_upload option, create-room-from-file endpoint, use generate_room_name()
6. `src/zndraw/api_manager.py` - Add 3 metadata API methods (get/set/delete)
7. `src/zndraw/zndraw.py` - Add metadata property returning RoomMetadata instance
8. `rest-endpoints.md` - Document new metadata and create-room-from-file endpoints
9. `README.md` - Add metadata usage examples and room naming explanation

### Key Features Added
- **Lock Protection**: Backend endpoints check `check_room_locked()` before mutations
- **Smart File Detection**: When loading a file, check if exact same file is already in a room
- **Fast Room Creation**: Create new rooms reusing physical storage with identity mapping (no re-upload!)
- **Flexible Options**: Open existing room OR create new room from same file OR force upload anyway
- **Smart Room Naming**: Auto-truncate long filenames with collision detection and hash suffixes
- **Collision Handling**: Automatically appends hash when room name already exists
- **Search Capability**: Regex search in both file browser and room list

## Design Principles Applied

### KISS (Keep It Simple, Stupid)
- Simple Redis hash for storage
- Standard MutableMapping interface
- Straightforward REST endpoints

### DRY (Don't Repeat Yourself)
- RoomMetadataManager encapsulates Redis operations
- Reusable helper functions for file usage tracking
- Shared validation logic

### SOLID
- **Single Responsibility**: Each manager class has one job
- **Open/Closed**: Easy to extend with new metadata fields
- **Liskov Substitution**: MutableMapping implementation is standard
- **Interface Segregation**: Clean separation of client/server concerns
- **Dependency Inversion**: Depends on abstractions (Redis, APIManager)

### YAGNI (You Aren't Gonna Need It)
- No complex metadata versioning (not required)
- No caching beyond simple lazy load (not needed yet)
- No metadata types beyond strings (sufficient for requirements)

## Implementation Order

1. **Phase 1**: Core backend infrastructure (manager + endpoints)
2. **Phase 2**: File upload integration (store file stats in metadata)
3. **Phase 3**: Python client API (`vis.metadata` property)
4. **Phase 4**: File browser duplicate detection (check before upload)
5. **Phase 5**: Search functionality (regex in file browser and room list)
6. **Phase 6**: Comprehensive testing
7. **Phase 7**: Documentation

This order ensures each phase builds on previous work and can be tested incrementally.

## Key Simplification
This plan focuses **only** on smart room creation without re-uploading files. We:
- ✅ Store file size and mtime when uploading
- ✅ Check if exact same file exists before uploading again
- ✅ Offer to **create new room reusing physical storage** (identity mapping)
- ✅ OR open existing room
- ✅ OR force upload anyway (with `force_upload: true` flag)
- ✅ **Auto-generate short, unique room names** from filenames
- ✅ Handle name collisions with hash suffixes
- ❌ Do NOT track file changes in existing rooms
- ❌ Do NOT update or warn about stale rooms
- ❌ Do NOT show change indicators in directory listings

**The Wins**: 
1. Create multiple rooms from the same file without re-uploading! Each room has independent trajectory (can delete frames, etc.) but shares the original physical storage.
2. Clean room names: "my_very_long_file_name.xyz" → "my_very_long_file_n" (or "my_very_long_f_a3c4d2" if collision)

