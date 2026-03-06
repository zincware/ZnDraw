# Constraint Visualization Design

## Problem

The old ZnDraw `room_service.py` created a `constraints-fixed-atoms` default geometry
using `InArrayTransform` to overlay red wireframe spheres on constrained atoms.
The new FastAPI codebase has the full `InArrayTransform` infrastructure (backend model,
frontend rendering, transform evaluation) but no default constraint geometry and no
user-friendly constraint selector in the TransformEditor.

## Decision

Enhance the existing system with two changes:

1. Add `constraints-fixed-atoms` to default room geometries
2. Enhance `TransformEditor` with a constraint selector (no free-solo)

## Design

### 1. Default Geometry: `constraints-fixed-atoms`

Added to `_initialize_default_geometries()` in `rooms.py`:

```python
"constraints-fixed-atoms": (
    "Sphere",
    {
        "active": True,
        "position": InArrayTransform(
            source="constraints",
            path="0.kwargs.indices",
            filter="arrays.positions",
        ).model_dump(),
        "radius": InArrayTransform(
            source="constraints",
            path="0.kwargs.indices",
            filter="arrays.radii",
        ).model_dump(),
        "color": ["#FF0000"],
        "material": {"wireframe": True},
        "scale": [[0.71, 0.71, 0.71]],
        "selecting": {"enabled": False},
        "hovering": {"enabled": False},
    },
),
```

- Red wireframe spheres at 71% scale on top of constrained atoms
- Targets first constraint (`0.kwargs.indices`) by default
- Non-interactive (selecting/hovering disabled)

### 2. TransformEditor Constraint Selector

When `source === "constraints"`, the TransformEditor:

1. Fetches constraint data from the current frame via `getFrames(roomId, currentFrame, ["constraints"])`
2. Displays a selectable list of constraints, e.g.:
   - `#0: FixAtoms - indices: [0, 2]`
   - `#1: FixedLine - indices: [1]`
3. On selection, sets `path` to `{index}.kwargs.indices`
4. No free-solo / manual editing - only selectable options
5. Graceful fallback when no constraints exist: "No constraints in current frame"
6. Selector only active for `source === "constraints"` - other sources keep plain text fields

### 3. Data Flow

```
Constraint list (asebytes encoded in frame):
[
  {"name": "FixAtoms", "kwargs": {"indices": [0, 2]}},
  {"name": "FixedLine", "kwargs": {"indices": [1], "direction": [1, 0, 0]}}
]

TransformEditor fetches this, shows:
  #0: FixAtoms — indices: [0, 2]
  #1: FixedLine — indices: [1]

User selects #0 → path = "0.kwargs.indices"

evaluateTransform():
  1. sourceData = frameData["constraints"] → list above
  2. extractPath(sourceData, "0.kwargs.indices") → [0, 2]
  3. filterData = frameData["arrays.positions"] → Float32Array
  4. Filter positions at indices [0, 2] → positions of fixed atoms only
```

## Testing

### Backend
- `_initialize_default_geometries` creates `constraints-fixed-atoms` with correct InArrayTransform config
- Constraint geometry config deserializes to valid Sphere model

### Playwright E2E
- Load room with FixAtoms constraints → red wireframe spheres visible
- Open geometry editor → TransformEditor shows constraint selector
- Select different constraint → path updates

## Documentation
- Python client example: build a molecule (e.g., acetic acid via molify/RDKit) where the
  carbon chain is constrained (FixAtoms) but the COOH group is free. Upload and show the
  red wireframe overlay on only the constrained atoms.
- How `constraints-fixed-atoms` geometry auto-visualizes
- Customization guide (color, scale, target constraint)
- Screenshot of the TransformEditor constraint selector with the molecule example

## Files to Change

| File | Change |
|------|--------|
| `src/zndraw/routes/rooms.py` | Add `constraints-fixed-atoms` to defaults |
| `frontend/src/components/jsonforms-renderers/TransformEditor.tsx` | Add constraint selector |
| `tests/test_constraints.py` | Add default geometry tests |
| `tests/e2e/` | Add Playwright constraint tests |
| `docs/` | Add constraint visualization docs |
