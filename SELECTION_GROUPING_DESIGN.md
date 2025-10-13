We have

```python
vis.geometries = {
    "particles": {
        "type": "Sphere",
        "position": "arrays.positions",
        "color": "arrays.colors"
    },
    "forces": {
        "type": "Arrow",
        "position": "arrays.positions",
        "direction": "calc.forces",
        "color": "#FF0000"
    },
    "bonds": {
        "type": "Cylinder",
        "indices": "arrays.bonds",  # Bond connectivity pairs
        "color": "#888888"
    },
    "boxes": {
        "type": "Box",
        "corners": [[0,0,0], [10,10,10]],
        "color": "#0000FF",
        "selectable": False  # Decorations usually not selectable
    }
}
```


```python
# Basic
vis.selection = [1, 2, 3]  # particles only

# Set / read other geometries
print(vis.selections["particles"])  # [1, 2, 3]
vis.selections["forces"] = [2, 3]
print(vis.selections["forces"])     # [2, 3]
print(vis.selection)                # [1, 2, 3] (particles unchanged)
print(vis.active_selection_group)   # None

# Named groups
vis.selection_groups = {
    "strong-forces": {"particles": [1, 3], "forces": [1, 3]}
}

# Load a group
vis.load_selection_group("strong-forces")

# Query
print(vis.selection)                # [1, 3]
print(vis.selections["forces"])     # [1, 3]
print(vis.active_selection_group)   # "strong-forces"
```

---

# Implementation Plan

## Architecture Overview

Following the existing patterns in the codebase (geometries, figures), selections will use:
- **REST API** for CRUD operations (easier to test, follows existing patterns)
- **Socket.IO invalidate events** for real-time updates across clients
- **Redis** for storage

This replaces the current `selection:set` / `selection:update` socket-only approach.

## Current State Analysis

**Existing Selection API:**
- `src/zndraw/app/events.py:268`: `selection:set` handler - stores in `room:{room}:selection:default`
- `src/zndraw/app/events.py:300`: `frame_selection:set` handler - stores in `room:{room}:frame_selection:default`
- Both emit `selection:update` / `frame_selection:update` to notify clients
- Frontend `app/src/hooks/useSocketManager.ts:348`: listens to `selection:update`

**Similar Patterns to Follow:**
- **Geometries**: REST API + `invalidate:geometry` (lines 2390-2513 in `routes.py`)
- **Figures**: REST API + `invalidate:figure` (lines 2532-2588 in `routes.py`)
- Both use targeted invalidation with optional `key` parameter

## Backend Implementation (Python)

### 1. Redis Schema (`src/zndraw/app/redis_keys.py` - new functions)

```python
def get_selections_key(room: str) -> str:
    """Hash storing current selections per geometry."""
    return f"room:{room}:selections"

def get_selection_groups_key(room: str) -> str:
    """Hash storing named selection groups."""
    return f"room:{room}:selection_groups"

def get_active_selection_group_key(room: str) -> str:
    """String storing the active selection group name."""
    return f"room:{room}:active_selection_group"
```

**Data Format:**
- `room:{room}:selections` - Hash: `{"particles": "[1,2,3]", "forces": "[2,3]"}`
- `room:{room}:selection_groups` - Hash: `{"group1": "{\"particles\": [1,3], \"forces\": [1,3]}"}`
- `room:{room}:active_selection_group` - String: `"group1"` or empty

### 2. REST API Routes (`src/zndraw/app/routes.py`)

Add after line 2588 (after figures routes):

```python
@main.route("/api/rooms/<string:room_id>/selections", methods=["GET"])
def get_all_selections(room_id: str):
    """Get all current selections and groups.

    Returns:
        {
            "selections": {"particles": [1,2,3], "forces": [2,3]},
            "groups": {"group1": {"particles": [1,3], "forces": [1,3]}},
            "activeGroup": "group1" | null
        }
    """
    r = current_app.extensions["redis"]

    # Get current selections
    selections_raw = r.hgetall(f"room:{room_id}:selections")
    selections = {k: json.loads(v) for k, v in selections_raw.items()}

    # Get selection groups
    groups_raw = r.hgetall(f"room:{room_id}:selection_groups")
    groups = {k: json.loads(v) for k, v in groups_raw.items()}

    # Get active group
    active_group = r.get(f"room:{room_id}:active_selection_group")

    return jsonify({
        "selections": selections,
        "groups": groups,
        "activeGroup": active_group
    })


@main.route("/api/rooms/<string:room_id>/selections/<string:geometry>", methods=["GET"])
def get_selection(room_id: str, geometry: str):
    """Get selection for a specific geometry."""
    r = current_app.extensions["redis"]
    selection = r.hget(f"room:{room_id}:selections", geometry)

    if selection is None:
        return jsonify({"selection": []}), 200

    return jsonify({"selection": json.loads(selection)})


@main.route("/api/rooms/<string:room_id>/selections/<string:geometry>", methods=["PUT"])
def update_selection(room_id: str, geometry: str):
    """Update selection for a specific geometry.

    Body: {"indices": [1, 2, 3]}
    """
    r = current_app.extensions["redis"]
    data = request.get_json()

    indices = data.get("indices", [])
    if not isinstance(indices, list):
        return jsonify({"error": "indices must be a list"}), 400

    if any(not isinstance(idx, int) or idx < 0 for idx in indices):
        return jsonify({"error": "All indices must be non-negative integers"}), 400

    # Store selection
    r.hset(f"room:{room_id}:selections", geometry, json.dumps(indices))

    # Clear active group (manual edit breaks group association)
    r.delete(f"room:{room_id}:active_selection_group")

    # Emit invalidation
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION,
        {"geometry": geometry},
        to=f"room:{room_id}"
    )

    return jsonify({"success": True})


@main.route("/api/rooms/<string:room_id>/selections/groups/<string:group_name>", methods=["GET"])
def get_selection_group(room_id: str, group_name: str):
    """Get a specific selection group."""
    r = current_app.extensions["redis"]
    group = r.hget(f"room:{room_id}:selection_groups", group_name)

    if group is None:
        return jsonify({"error": "Group not found"}), 404

    return jsonify({"group": json.loads(group)})


@main.route("/api/rooms/<string:room_id>/selections/groups/<string:group_name>", methods=["PUT"])
def create_update_selection_group(room_id: str, group_name: str):
    """Create or update a selection group.

    Body: {"particles": [1, 3], "forces": [1, 3]}
    """
    r = current_app.extensions["redis"]
    data = request.get_json()

    # Validate data is a dict of geometry -> indices
    if not isinstance(data, dict):
        return jsonify({"error": "Group must be a dictionary"}), 400

    for geometry, indices in data.items():
        if not isinstance(indices, list):
            return jsonify({"error": f"Indices for '{geometry}' must be a list"}), 400
        if any(not isinstance(idx, int) or idx < 0 for idx in indices):
            return jsonify({"error": f"Invalid indices for '{geometry}'"}), 400

    # Store group
    r.hset(f"room:{room_id}:selection_groups", group_name, json.dumps(data))

    # Emit invalidation (groups list changed)
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION,
        {"operation": "groups_changed"},
        to=f"room:{room_id}"
    )

    return jsonify({"success": True})


@main.route("/api/rooms/<string:room_id>/selections/groups/<string:group_name>", methods=["DELETE"])
def delete_selection_group(room_id: str, group_name: str):
    """Delete a selection group."""
    r = current_app.extensions["redis"]

    # Check if group exists
    if not r.hexists(f"room:{room_id}:selection_groups", group_name):
        return jsonify({"error": "Group not found"}), 404

    # Delete group
    r.hdel(f"room:{room_id}:selection_groups", group_name)

    # Clear active group if it was this one
    active_group = r.get(f"room:{room_id}:active_selection_group")
    if active_group == group_name:
        r.delete(f"room:{room_id}:active_selection_group")

    # Emit invalidation
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION,
        {"operation": "groups_changed"},
        to=f"room:{room_id}"
    )

    return jsonify({"success": True})


@main.route("/api/rooms/<string:room_id>/selections/groups/<string:group_name>/load", methods=["POST"])
def load_selection_group(room_id: str, group_name: str):
    """Load a selection group (apply it to current selections).

    This sets the active group and updates all selections to match the group.
    """
    r = current_app.extensions["redis"]

    # Get group
    group_data = r.hget(f"room:{room_id}:selection_groups", group_name)
    if group_data is None:
        return jsonify({"error": "Group not found"}), 404

    group = json.loads(group_data)

    # Apply group to current selections
    for geometry, indices in group.items():
        r.hset(f"room:{room_id}:selections", geometry, json.dumps(indices))

    # Set as active group
    r.set(f"room:{room_id}:active_selection_group", group_name)

    # Emit invalidation (all selections changed)
    socketio.emit(
        SocketEvents.INVALIDATE_SELECTION,
        {"operation": "group_loaded", "group": group_name},
        to=f"room:{room_id}"
    )

    return jsonify({"success": True})
```

### 3. Socket Events Updates

**`src/zndraw/app/constants.py`** - Add new constant:
```python
class SocketEvents:
    # ... existing events ...
    INVALIDATE_SELECTION = "invalidate:selection"
```

**`src/zndraw/app/events.py`** - Remove old handlers (lines 268-329):
- remove `selection:set` and `frame_selection:set`

### 4. Python Client (`src/zndraw/zndraw.py`)

Create new selection manager classes:

```python
class Selections(MutableMapping):
    """Accessor for per-geometry selections."""

    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def __getitem__(self, geometry: str) -> frozenset[int]:
        response = self.vis.api.get_selection(geometry)
        return frozenset(response.get("selection", []))

    def __setitem__(self, geometry: str, indices: t.Iterable[int]) -> None:
        self.vis.api.update_selection(geometry, list(indices))

    def __delitem__(self, geometry: str) -> None:
        self.vis.api.update_selection(geometry, [])

    def __iter__(self):
        data = self.vis.api.get_all_selections()
        return iter(data["selections"].keys())

    def __len__(self) -> int:
        data = self.vis.api.get_all_selections()
        return len(data["selections"])


class SelectionGroups(MutableMapping):
    """Accessor for named selection groups."""

    def __init__(self, zndraw_instance: "ZnDraw") -> None:
        self.vis = zndraw_instance

    def __getitem__(self, group_name: str) -> dict[str, list[int]]:
        response = self.vis.api.get_selection_group(group_name)
        return response["group"]

    def __setitem__(self, group_name: str, group_data: dict[str, t.Iterable[int]]) -> None:
        # Convert iterables to lists
        data = {k: list(v) for k, v in group_data.items()}
        self.vis.api.create_update_selection_group(group_name, data)

    def __delitem__(self, group_name: str) -> None:
        self.vis.api.delete_selection_group(group_name)

    def __iter__(self):
        data = self.vis.api.get_all_selections()
        return iter(data["groups"].keys())

    def __len__(self) -> int:
        data = self.vis.api.get_all_selections()
        return len(data["groups"])


class ZnDraw(MutableSequence):
    # ... existing code ...

    @property
    def selections(self) -> Selections:
        """Access selections by geometry name."""
        if not hasattr(self, "_selections_accessor"):
            self._selections_accessor = Selections(self)
        return self._selections_accessor

    @property
    def selection_groups(self) -> SelectionGroups:
        """Access named selection groups."""
        if not hasattr(self, "_selection_groups_accessor"):
            self._selection_groups_accessor = SelectionGroups(self)
        return self._selection_groups_accessor

    @property
    def active_selection_group(self) -> str | None:
        """Get the currently active selection group name."""
        data = self.api.get_all_selections()
        return data.get("activeGroup")

    def load_selection_group(self, group_name: str) -> None:
        """Load a selection group (apply it to current selections)."""
        self.api.load_selection_group(group_name)

    @property
    def selection(self) -> frozenset[int]:
        """Get selection for 'particles' geometry (backward compatible)."""
        return self.selections["particles"]

    @selection.setter
    def selection(self, value: t.Iterable[int] | None):
        """Set selection for 'particles' geometry (backward compatible)."""
        if value is None:
            value = []
        self.selections["particles"] = value
```

### 5. API Client Updates (`src/zndraw/api.py` or similar)

Add new API methods to the HTTP client class.

## Frontend Implementation (TypeScript/React)

### 1. API Client (`app/src/myapi/client.ts`)

Add new functions:

```typescript
export async function getAllSelections(roomId: string) {
  const response = await fetch(`/api/rooms/${roomId}/selections`);
  return await response.json();
}

export async function getSelection(roomId: string, geometry: string) {
  const response = await fetch(`/api/rooms/${roomId}/selections/${geometry}`);
  return await response.json();
}

export async function updateSelection(roomId: string, geometry: string, indices: number[]) {
  const response = await fetch(`/api/rooms/${roomId}/selections/${geometry}`, {
    method: "PUT",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ indices }),
  });
  return await response.json();
}

export async function getSelectionGroup(roomId: string, groupName: string) {
  const response = await fetch(`/api/rooms/${roomId}/selections/groups/${groupName}`);
  return await response.json();
}

export async function createUpdateSelectionGroup(
  roomId: string,
  groupName: string,
  groupData: Record<string, number[]>
) {
  const response = await fetch(`/api/rooms/${roomId}/selections/groups/${groupName}`, {
    method: "PUT",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(groupData),
  });
  return await response.json();
}

export async function deleteSelectionGroup(roomId: string, groupName: string) {
  const response = await fetch(`/api/rooms/${roomId}/selections/groups/${groupName}`, {
    method: "DELETE",
  });
  return await response.json();
}

export async function loadSelectionGroup(roomId: string, groupName: string) {
  const response = await fetch(`/api/rooms/${roomId}/selections/groups/${groupName}/load`, {
    method: "POST",
  });
  return await response.json();
}
```

### 2. Store Updates (`app/src/store.tsx`)

Update selection state:

```typescript
interface AppStore {
  // Old (keep for backward compat during transition)
  selection: number[] | null;
  setSelection: (selection: number[] | null) => void;

  // New
  selections: Record<string, number[]>;  // {"particles": [1,2,3], "forces": [2,3]}
  selectionGroups: Record<string, Record<string, number[]>>;
  activeSelectionGroup: string | null;
  setSelections: (selections: Record<string, number[]>) => void;
  setSelectionGroups: (groups: Record<string, Record<string, number[]>>) => void;
  setActiveSelectionGroup: (group: string | null) => void;
  // ... rest
}
```

### 3. Socket Manager Updates (`app/src/hooks/useSocketManager.ts`)

Replace lines 348-349:

```typescript
// Remove old handler:
// socket.on("selection:update", onSelectionUpdate);

// Add new handler:
async function onSelectionInvalidate(data: any) {
  if (!roomId) return;

  try {
    if (data.operation === "group_loaded") {
      // Full refresh - all selections changed
      const response = await getAllSelections(roomId);
      setSelections(response.selections);
      setSelectionGroups(response.groups);
      setActiveSelectionGroup(response.activeGroup);
    } else if (data.operation === "groups_changed") {
      // Only groups list changed
      const response = await getAllSelections(roomId);
      setSelectionGroups(response.groups);
    } else if (data.geometry) {
      // Specific geometry selection changed
      const response = await getSelection(roomId, data.geometry);
      setSelections((prev) => ({
        ...prev,
        [data.geometry]: response.selection,
      }));
      // Clear active group if not from group load
      setActiveSelectionGroup(null);
    }
  } catch (error) {
    console.error("Error handling selection invalidation:", error);
  }
}

socket.on("invalidate:selection", onSelectionInvalidate);
```

### 4. UI Component Updates

Update selection-related components to:
- Use `selections` instead of `selection` for multi-geometry support
- Add UI for managing selection groups (optional, can be deferred)
- Ensure keeping `vis.selection` for selecting on the "particles" geometry

## Migration Strategy

### Phase 1: Add New API (Backward Compatible)
1. Add new REST endpoints alongside existing socket handlers
2. Add new socket event `invalidate:selection`
3. Remove old `selection:set` / `selection:update`
4. Update backend ZnDraw client to use new API, update `selection` property to work with new structure

### Phase 2: Frontend Transition
1. Add new frontend API client functions
2. Update store to support new structure
3. Add `invalidate:selection` handler
4. remove old `selection:update` handler
5. Update UI components gradually

### Phase 3:
1. remove old socket events, if still there

## Testing Strategy

### Backend Tests (`tests/`)

```python
def test_get_all_selections(client, room):
    """Test GET /api/rooms/{room}/selections"""
    response = client.get(f"/api/rooms/{room}/selections")
    assert response.status_code == 200
    data = response.json
    assert "selections" in data
    assert "groups" in data
    assert "activeGroup" in data

def test_update_selection(client, room):
    """Test PUT /api/rooms/{room}/selections/{geometry}"""
    response = client.put(
        f"/api/rooms/{room}/selections/particles",
        json={"indices": [1, 2, 3]}
    )
    assert response.status_code == 200

    # Verify stored correctly
    response = client.get(f"/api/rooms/{room}/selections/particles")
    assert response.json["selection"] == [1, 2, 3]

def test_selection_groups_crud(client, room):
    """Test selection groups create/read/update/delete"""
    # Create group
    group_data = {"particles": [1, 3], "forces": [1, 3]}
    response = client.put(
        f"/api/rooms/{room}/selections/groups/test-group",
        json=group_data
    )
    assert response.status_code == 200

    # Read group
    response = client.get(f"/api/rooms/{room}/selections/groups/test-group")
    assert response.json["group"] == group_data

    # Delete group
    response = client.delete(f"/api/rooms/{room}/selections/groups/test-group")
    assert response.status_code == 200

def test_load_selection_group(client, room):
    """Test loading a selection group"""
    # Create group
    group_data = {"particles": [1, 3], "forces": [1, 3]}
    client.put(f"/api/rooms/{room}/selections/groups/test-group", json=group_data)

    # Load group
    response = client.post(f"/api/rooms/{room}/selections/groups/test-group/load")
    assert response.status_code == 200

    # Verify selections match group
    response = client.get(f"/api/rooms/{room}/selections/particles")
    assert response.json["selection"] == [1, 3]

    # Verify active group is set
    response = client.get(f"/api/rooms/{room}/selections")
    assert response.json["activeGroup"] == "test-group"

def test_vis_selections_accessor(s22, server):
    """Test ZnDraw.selections dict-like interface"""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Set selections for different geometries
    vis.selections["particles"] = [1, 2, 3]
    vis.selections["forces"] = [2, 3]

    assert vis.selections["particles"] == frozenset([1, 2, 3])
    assert vis.selections["forces"] == frozenset([2, 3])

def test_vis_selection_groups(s22, server):
    """Test ZnDraw selection groups"""
    vis = ZnDraw(url=server, room="testroom", user="testuser")
    vis.extend(s22)

    # Create group
    vis.selection_groups["group1"] = {
        "particles": [1, 3],
        "forces": [1, 3]
    }

    # Load group
    vis.load_selection_group("group1")

    # Verify active
    assert vis.active_selection_group == "group1"
    assert vis.selections["particles"] == frozenset([1, 3])
    assert vis.selections["forces"] == frozenset([1, 3])
```

### Frontend Tests

Use existing test patterns with REST API mocking.

## Files to Create/Modify

### Backend (Python)
- [ ] `src/zndraw/app/routes.py` - Add selection REST API endpoints
- [ ] `src/zndraw/app/constants.py` - Add `INVALIDATE_SELECTION` constant
- [ ] `src/zndraw/app/events.py` - Mark old handlers as deprecated, keep for compatibility
- [ ] `src/zndraw/zndraw.py` - Add `selections`, `selection_groups`, `active_selection_group` properties
- [ ] `src/zndraw/api.py` (or client API file) - Add REST client methods
- [ ] `tests/test_vis_selections.py` - Create comprehensive tests

### Frontend (TypeScript/React)
- [ ] `app/src/myapi/client.ts` - Add selection API functions
- [ ] `app/src/store.tsx` - Update selection state structure
- [ ] `app/src/hooks/useSocketManager.ts` - Replace `selection:update` with `invalidate:selection`
- [ ] Selection UI components - Update to use new API (as needed)
