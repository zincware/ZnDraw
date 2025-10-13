# Selection Groups UI and Multi-Geometry Selection Visualization Plan

## Overview
Extend the ZnDraw selection system with:
1. A new `INVALIDATE_SELECTION_GROUPS` socket event for efficient group-only updates
2. A sidebar UI for managing selection groups
3. Visual selection indicators for Bonds and Arrows (like Particles already have)
4. Update Zustand store to properly manage selection state for different geometries

## Part 1: Backend Changes

### 1.1 Add INVALIDATE_SELECTION_GROUPS Constant
**File**: `src/zndraw/app/constants.py`

Add new constant:
```python
INVALIDATE_SELECTION_GROUPS = "invalidate:selection_groups"
```

**Rationale**: Separate event for group changes allows clients to update only the groups list without refetching all selections. This is more efficient when only group metadata changes (create/delete/rename).

### 1.2 Update Backend Routes
**File**: `src/zndraw/app/routes.py`

Update the following endpoints to emit `INVALIDATE_SELECTION_GROUPS`:
- **PUT** `/api/rooms/{room_id}/selections/groups/{group_name}` - Emit after creating/updating group
- **DELETE** `/api/rooms/{room_id}/selections/groups/{group_name}` - Emit after deleting group

Keep `INVALIDATE_SELECTION` for:
- **PUT** `/api/rooms/{room_id}/selections/{geometry}` - Selection changes
- **POST** `/api/rooms/{room_id}/selections/groups/{group_name}/load` - Loading groups (changes selections)

**Rationale**: Different events for different operations allows frontend to handle them differently. Groups-only changes don't need to refetch selections.

## Part 2: Frontend Store Updates

### 2.1 Zustand Store Strategy

**Current Situation**:
- `selection: number[] | null` - Backward compatibility, maps to `selections["particles"]`
- `selections: Record<string, number[]>` - Per-geometry selections
- `selectionGroups: Record<string, Record<string, number[]>>` - Named groups
- `activeSelectionGroup: string | null` - Currently active group

**Components Currently Using Selection**:
1. **Particles.tsx** (line 57): `const { selection } = useAppStore();`
2. **SingleBonds.tsx** (line 63): `const { selection } = useAppStore();`
3. **Arrow.tsx** (line 64): `const { selection } = useAppStore();`

**Problem**: All components currently read from the single `selection` state which maps to "particles" geometry. We need components to read from their own geometry's selection.

**Solution**:
- **Option A (Minimal Change)**: Add a helper function in store: `getSelectionForGeometry(geometryKey: string)`
- **Option B (Better Design)**: Each component reads `selections[geometryKey]` directly

**Recommendation**: Use Option B for cleaner code:

```typescript
// In Particles.tsx
const selections = useAppStore(state => state.selections);
const particleSelection = selections["particles"] || [];

// In SingleBonds.tsx
const selections = useAppStore(state => state.selections);
const bondSelection = selections["bonds"] || [];

// In Arrow.tsx
const selections = useAppStore(state => state.selections);
const arrowSelection = selections["arrows"] || [];
```

### 2.2 Socket Manager Updates
**File**: `app/src/hooks/useSocketManager.ts`

Add new handler for `invalidate:selection_groups`:

```typescript
async function onSelectionGroupsInvalidate(data: any) {
  if (!roomId) return;

  try {
    console.log("Received invalidate:selection_groups event:", data);

    // Fetch only the groups (more efficient than fetching everything)
    const response = await getAllSelections(roomId);

    // Update only groups and active group (not selections)
    setSelectionGroups(response.groups);
    setActiveSelectionGroup(response.activeGroup);

    console.log("Updated selection groups from server:", response.groups);
  } catch (error) {
    console.error("Error fetching selection groups after invalidation:", error);
  }
}
```

Register the handler:
```typescript
socket.on("invalidate:selection_groups", onSelectionGroupsInvalidate);
```

## Part 3: Selection Visualization for Geometries

### 3.1 Particles (Already Implemented ✓)
**Current Implementation**: Lines 174-186 in Particles.tsx
- Separate `selectionMeshRef` for selection visualization
- Uses `SELECTION_SCALE` constant for sizing
- Renders with semi-transparent overlay

**No changes needed** - this is the reference implementation.

### 3.2 Bonds Selection Visualization
**File**: `app/src/components/three/SingleBonds.tsx`

**Current State**: Lines 276-286 have commented-out selection/hover meshes

**Implementation Plan**:

1. **Track Selected Bonds** (add near line 63):
```typescript
const selections = useAppStore(state => state.selections);
const bondSelection = selections["bonds"] || [];
const bondSelectionSet = useMemo(() => new Set(bondSelection), [bondSelection]);
```

2. **Calculate Selected Bond Instances**:
```typescript
// Bonds are rendered as 2 instances per bond (half-bonds)
// Need to map particle selection to bond instances
const selectedBondIndices = useMemo(() => {
  const indices: number[] = [];
  // Logic to determine which bond instances are selected
  // based on whether either endpoint particle is selected
  return indices;
}, [bondSelectionSet, /* bond data */]);
```

3. **Enable Selection Mesh** (uncomment and update lines 276-286):
```typescript
{selecting.enabled && selectedBondIndices.length > 0 && (
  <instancedMesh
    key={`selection-${selectedBondIndices.length}`}
    ref={selectionMeshRef}
    args={[bondGeometry, undefined, selectedBondIndices.length]}
  >
    <meshBasicMaterial
      transparent
      opacity={selecting.opacity}
      color={selecting.color || "#FFFF00"}
      side={THREE.BackSide}
    />
  </instancedMesh>
)}
```

4. **Update Selection Mesh in useEffect** (add to mesh update logic):
```typescript
// After main mesh update
if (selecting.enabled && selectionMeshRef.current && selectedBondIndices.length > 0) {
  const selectionMesh = selectionMeshRef.current;
  selectedBondIndices.forEach((instanceIdx, arrayIdx) => {
    // Copy matrix from main mesh and scale up
    mainMesh.getMatrixAt(instanceIdx, _matrix);
    _matrix.scale(_vec3.set(SELECTION_SCALE, SELECTION_SCALE, SELECTION_SCALE));
    selectionMesh.setMatrixAt(arrayIdx, _matrix);
  });
  selectionMesh.instanceMatrix.needsUpdate = true;
}
```

**Key Decision**: How to determine bond selection?
- **Option A**: Bond selected if EITHER endpoint is selected
- **Option B**: Bond selected if BOTH endpoints are selected
- **Recommendation**: Option A (more intuitive - selecting particle shows connected bonds)

### 3.3 Arrows Selection Visualization
**File**: `app/src/components/three/Arrow.tsx`

**Current State**: Lines 213-217 change color inline when selected

**Two Approaches**:

**Approach A: Keep Inline Color Change (Simplest)**
- Update to read from `selections["arrows"]` instead of generic `selection`
- No mesh changes needed

```typescript
const selections = useAppStore(state => state.selections);
const arrowSelection = selections["arrows"] || [];
const arrowSelectionSet = useMemo(() => new Set(arrowSelection), [arrowSelection]);

// In update loop (line 213)
if (arrowSelectionSet.has(i)) {
  _color.setRGB(SELECTION_COLOR[0], SELECTION_COLOR[1], SELECTION_COLOR[2]);
} else {
  _color.setRGB(finalColors[i3], finalColors[i3 + 1], finalColors[i3 + 2]);
}
```

**Approach B: Add Selection Mesh (Like Particles)**
- Create separate selection mesh with scaled arrows
- More consistent with Particles
- Slightly more complex

**Recommendation**: Start with Approach A for consistency with current implementation. Can upgrade to Approach B later if needed.

## Part 4: Selection Groups Sidebar UI

### 4.1 Component Structure

Create new component: `app/src/components/selection/SelectionGroupsPanel.tsx`

**Features**:
1. **Group List** - Display all available selection groups
2. **Active Group Indicator** - Show which group is currently loaded
3. **Group Actions**:
   - Load/Apply group (updates all geometry selections)
   - Rename group
   - Delete group
4. **Create New Group** - Save current selections as a new named group
5. **Group Details** - Show which geometries have selections in each group

### 4.2 UI Layout

```
┌─────────────────────────────────────┐
│ Selection Groups                    │
├─────────────────────────────────────┤
│ ┌─────────────────────────────────┐ │
│ │ Create New Group                │ │
│ │ [Current Selection] [+ Save]    │ │
│ └─────────────────────────────────┘ │
│                                     │
│ Active Groups:                      │
│ ┌─────────────────────────────────┐ │
│ │ ● Group 1          [Load][Del]  │ │ ← Active
│ │   particles: 5, bonds: 3        │ │
│ ├─────────────────────────────────┤ │
│ │ ○ Group 2          [Load][Del]  │ │
│ │   particles: 2, arrows: 4       │ │
│ ├─────────────────────────────────┤ │
│ │ ○ Group 3          [Load][Del]  │ │
│ │   particles: 10                 │ │
│ └─────────────────────────────────┘ │
└─────────────────────────────────────┘
```

### 4.3 Component Implementation

```typescript
import { useState } from "react";
import {
  Box,
  Button,
  TextField,
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  IconButton,
  Typography,
  Chip,
  Stack,
} from "@mui/material";
import DeleteIcon from "@mui/icons-material/Delete";
import PlayArrowIcon from "@mui/icons-material/PlayArrow";
import RadioButtonCheckedIcon from "@mui/icons-material/RadioButtonChecked";
import RadioButtonUncheckedIcon from "@mui/icons-material/RadioButtonUnchecked";
import { useAppStore } from "../../store";
import { createUpdateSelectionGroup, deleteSelectionGroup } from "../../myapi/client";

export default function SelectionGroupsPanel() {
  const {
    roomId,
    selections,
    selectionGroups,
    activeSelectionGroup,
    loadSelectionGroup,
    setSelectionGroups,
  } = useAppStore();

  const [newGroupName, setNewGroupName] = useState("");

  const handleSaveCurrentSelection = async () => {
    if (!roomId || !newGroupName.trim()) return;

    try {
      // Save current selections as a new group
      await createUpdateSelectionGroup(roomId, newGroupName, selections);
      setNewGroupName("");
      // Groups will update via socket invalidation
    } catch (error) {
      console.error("Failed to save selection group:", error);
    }
  };

  const handleLoadGroup = (groupName: string) => {
    loadSelectionGroup(groupName);
  };

  const handleDeleteGroup = async (groupName: string) => {
    if (!roomId) return;

    try {
      await deleteSelectionGroup(roomId, groupName);
      // Groups will update via socket invalidation
    } catch (error) {
      console.error("Failed to delete selection group:", error);
    }
  };

  // Count geometries with selections in each group
  const getGroupSummary = (group: Record<string, number[]>) => {
    return Object.entries(group)
      .filter(([_, indices]) => indices.length > 0)
      .map(([geom, indices]) => `${geom}: ${indices.length}`)
      .join(", ");
  };

  const hasCurrentSelection = Object.values(selections).some(
    (indices) => indices.length > 0
  );

  return (
    <Box sx={{ p: 2, height: "100%", overflow: "auto" }}>
      <Typography variant="h6" gutterBottom>
        Selection Groups
      </Typography>

      {/* Create New Group */}
      <Box sx={{ mb: 3, p: 2, bgcolor: "background.default", borderRadius: 1 }}>
        <Typography variant="subtitle2" gutterBottom>
          Save Current Selection
        </Typography>
        <Stack direction="row" spacing={1}>
          <TextField
            size="small"
            fullWidth
            placeholder="Group name"
            value={newGroupName}
            onChange={(e) => setNewGroupName(e.target.value)}
            disabled={!hasCurrentSelection}
          />
          <Button
            variant="contained"
            onClick={handleSaveCurrentSelection}
            disabled={!newGroupName.trim() || !hasCurrentSelection}
          >
            Save
          </Button>
        </Stack>
        {!hasCurrentSelection && (
          <Typography variant="caption" color="text.secondary" sx={{ mt: 1 }}>
            Select some elements first to save a group
          </Typography>
        )}
      </Box>

      {/* Groups List */}
      <Typography variant="subtitle2" gutterBottom>
        Saved Groups
      </Typography>
      <List>
        {Object.entries(selectionGroups).map(([groupName, groupData]) => (
          <ListItem
            key={groupName}
            sx={{
              border: 1,
              borderColor: "divider",
              borderRadius: 1,
              mb: 1,
              bgcolor:
                activeSelectionGroup === groupName
                  ? "action.selected"
                  : "background.paper",
            }}
          >
            <ListItemIcon>
              {activeSelectionGroup === groupName ? (
                <RadioButtonCheckedIcon color="primary" />
              ) : (
                <RadioButtonUncheckedIcon />
              )}
            </ListItemIcon>
            <ListItemText
              primary={groupName}
              secondary={getGroupSummary(groupData)}
            />
            <IconButton
              size="small"
              onClick={() => handleLoadGroup(groupName)}
              title="Load this group"
            >
              <PlayArrowIcon />
            </IconButton>
            <IconButton
              size="small"
              onClick={() => handleDeleteGroup(groupName)}
              title="Delete group"
            >
              <DeleteIcon />
            </IconButton>
          </ListItem>
        ))}
        {Object.keys(selectionGroups).length === 0 && (
          <Typography variant="body2" color="text.secondary" sx={{ mt: 2, textAlign: "center" }}>
            No saved groups yet
          </Typography>
        )}
      </List>
    </Box>
  );
}
```

### 4.4 Integrate into Sidebar

**File**: `app/src/components/SideBar.tsx`

Update navItems to include selection groups:
```typescript
{
  name: "selection-groups",
  icon: <BookmarkIcon />,  // Or another appropriate icon
  schemaType: "selection-groups",
  description: "Manage selection groups"
}
```

Update render logic:
```typescript
{selectedCategory === "geometries" ? (
  <GeometryPanel />
) : selectedCategory === "selection-groups" ? (
  <SelectionGroupsPanel />
) : (
  <SecondaryPanel
    key={selectedCategory}
    panelTitle={selectedCategory}
  />
)}
```

## Part 5: Testing Strategy

### 5.1 Manual Testing Checklist

**Selection Visualization**:
- [ ] Particles show selection overlay (already working)
- [ ] Bonds show selection overlay when particles are selected
- [ ] Arrows change color when selected
- [ ] Selection persists across frame changes
- [ ] Multiple geometries can be selected simultaneously

**Selection Groups UI**:
- [ ] Can create new group from current selection
- [ ] Can load existing group
- [ ] Loading group updates all geometry selections
- [ ] Active group indicator shows correctly
- [ ] Can delete groups
- [ ] Groups persist across page reload
- [ ] Groups sync between multiple clients

**Socket Events**:
- [ ] `invalidate:selection` updates selections in other clients
- [ ] `invalidate:selection_groups` updates groups in other clients
- [ ] Loading group triggers selection invalidation

### 5.2 Automated Tests

Add to `tests/test_vis_selections.py`:
```python
def test_invalidate_selection_groups_event(s22, server):
    """Test that INVALIDATE_SELECTION_GROUPS event is emitted correctly."""
    # TODO: Implement socket event testing
    pass
```

## Part 6: Implementation Order

1. **Backend** (Quick wins first):
   - [ ] Add `INVALIDATE_SELECTION_GROUPS` constant
   - [ ] Update routes to emit new event
   - [ ] Test with existing backend tests

2. **Frontend Store** (Foundation):
   - [ ] Add `invalidate:selection_groups` handler to useSocketManager
   - [ ] Update components to read from `selections[geometryKey]`

3. **Selection Visualization**:
   - [ ] Update Bonds to show selection mesh
   - [ ] Update Arrows to use geometry-specific selection
   - [ ] Test visual feedback

4. **Selection Groups UI** (User-facing):
   - [ ] Create SelectionGroupsPanel component
   - [ ] Integrate into Sidebar
   - [ ] Polish UI/UX

5. **Testing & Refinement**:
   - [ ] Manual testing
   - [ ] Bug fixes
   - [ ] Performance optimization if needed

## Key Design Decisions Summary

1. **Separate Socket Events**: `INVALIDATE_SELECTION` for selections, `INVALIDATE_SELECTION_GROUPS` for groups - allows efficient updates
2. **Per-Geometry Selection**: Each geometry reads `selections[geometryKey]` for its own selection state
3. **Selection Visualization Approach**:
   - Particles: Separate mesh (current)
   - Bonds: Separate mesh (new - matches particles)
   - Arrows: Color change (current - simpler)
4. **Bond Selection Logic**: Bond selected if EITHER endpoint is selected (more intuitive)
5. **UI Placement**: Separate sidebar panel for selection groups management

## Performance Considerations

1. **Socket Events**: Separate events prevent unnecessary data fetching
2. **Zustand Selectors**: Use selective subscriptions to prevent re-renders
3. **Selection Mesh Updates**: Only update when selection changes, not every frame
4. **Memoization**: Use useMemo for selection sets and derived data
