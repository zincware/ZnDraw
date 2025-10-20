# Plotly Figure Interactions - Implementation Summary

This document summarizes the implementation of bidirectional synchronization between Plotly plots and 3D geometry in ZnDraw.

## Overview

Users can now create interactive Plotly figures where:
- **Clicking/selecting points in plots** updates frame, particle selection, or other state
- **State changes in 3D view** highlight corresponding points in plots

The system uses an explicit `meta.interactions` schema in Plotly figures to control how interactions work.

## Files Modified

### Backend (Python)

#### 1. `src/zndraw/figures_manager.py`
- **Added extensive docstring** to `Figures` class documenting the new interactive features with examples
- **Added `validate_interaction_schema()` function** to validate customdata and interaction schemas
  - Checks dimension alignment between customdata and schema
  - Validates action names (must be "step" or valid geometry names)
  - Provides detailed error messages for debugging
  - Handles numpy arrays and lists

### Frontend (TypeScript/React)

#### 1. `app/src/store.tsx` (Zustand Store)
- **Added `hoveredFrame` state** to track hovered frames from plot interactions
- **Added `setHoveredFrame()` action** to update hovered frame state

#### 2. `app/src/components/FigureWindow.tsx` (Complete Rewrite)
- **Implemented interaction schema parsing** - reads `meta.interactions` from figures
- **Implemented schema validation** - checks customdata dimensions match schema
- **Implemented batch action collection** - collects actions from clicked/selected points
- **Implemented Plotly event bindings:**
  - `plotly_click` - handle point clicks
  - `plotly_selected` - handle lasso/box selections
  - `plotly_deselect` - handle deselection
  - `plotly_hover` - handle hover events
  - `plotly_unhover` - handle hover exit
- **Implemented batch processing** - atomically applies all state changes
- **Implemented state synchronization:**
  - Built customdata → point index lookup map for O(1) lookups
  - Implemented consolidated state effect that highlights points based on:
    - Current frame selection
    - Frame selection range
    - Geometry selections
  - Added debounced plot highlighting (100ms) for performance
- **Implemented error handling:**
  - Displays validation errors in alert pill
  - Prevents state updates on invalid schemas
  - Figure still renders even with errors

## Key Features

### 1. Flexible Schema System

Users completely control data structure via `customdata` and define interaction semantics in `meta.interactions`:

```python
fig.update_traces(
    customdata=np.column_stack([steps, particle_ids]),
    meta={
        "interactions": [
            {"click": "step", "select": "step"},      # Dimension 0
            {"click": "particles", "select": "particles"}  # Dimension 1
        ]
    }
)
```

### 2. Standard Action Names

- **`"step"`** - Controls frame (click jumps, select adds to selection, hover shows feedback)
- **`"<geometry_name>"`** - Maps to geometries in vis.geometries (particles, forces, bonds, etc.)

### 3. Bidirectional Sync

- **Plot → 3D:** Interactions update store state and sync with 3D view
- **3D → Plot:** State changes trigger plot highlighting of corresponding points

### 4. Graceful Degradation

- Plots without interaction schemas work normally (no interactions)
- Invalid schemas display clear error messages
- Figure continues rendering even with validation errors

### 5. Performance Optimized

- O(1) customdata lookups using Map
- Debounced plot updates (100ms) to batch multiple changes
- Efficient opacity-based highlighting

## Examples Provided

### 1. `examples/plotly_interactions_demo.py`

Demonstrates four common usage patterns:

1. **Step mapping only** - Energy prediction plot with frame selection
2. **Geometry mapping** - Force distribution histogram with geometry selection
3. **Multi-dimensional** - Particle forces with both step and geometry selection
4. **Sparse interactions** - Selective actions (click only vs select only)

## Documentation

### 1. `docs/PLOTLY_INTERACTIONS.md`

Comprehensive guide including:
- Quick start example
- Core concepts explanation
- Four detailed examples with expected behavior
- Validation guide
- Important design constraints
- Performance considerations
- Troubleshooting guide
- API reference

### 2. Inline Documentation

- Extensive docstring in Figures class with multiple examples
- Type definitions and interfaces throughout FigureWindow component
- Helper function documentation explaining algorithm complexity

## Testing

### Backend Validation Tests

All validation tests pass:
- ✅ Valid 1D schema
- ✅ Valid 2D schema
- ✅ Dimension mismatch detection
- ✅ Invalid geometry detection

### Frontend TypeScript

- ✅ All custom TypeScript errors fixed
- ✅ Only pre-existing @types/react-plotly.js issue remains (unrelated)

## Schema Validation Examples

```python
from zndraw.figures_manager import validate_interaction_schema

# Valid: Simple step mapping
customdata = np.arange(100)
schema = [{"click": "step", "select": "step"}]
is_valid, errors = validate_interaction_schema(customdata, schema)
# Returns: (True, [])

# Invalid: Dimension mismatch
customdata = np.column_stack([steps, particles])  # 2D
schema = [{"click": "step"}]  # Only 1 schema entry
is_valid, errors = validate_interaction_schema(customdata, schema)
# Returns: (False, ["Dimension mismatch: customdata has 2 dimension(s), schema has 1"])

# Invalid: Unknown geometry
customdata = np.arange(100)
schema = [{"click": "unknown_geometry"}]
geometries = {"particles": {}, "forces": {}}
is_valid, errors = validate_interaction_schema(customdata, schema, geometries)
# Returns: (False, ["'unknown_geometry' not found in vis.geometries..."])
```

## Design Decisions

### 1. Explicit Schema vs. Convention

**Decision:** Use explicit `meta.interactions` schema rather than assuming fixed structure
**Rationale:** Supports arbitrary plot types (scatter, histogram, heatmap, etc.) and data structures

### 2. Customdata Alignment

**Decision:** Customdata values should map directly to geometry selection indices
**Rationale:** Simplest approach, users understand this is how selections work

### 3. Per-Trace Interactions

**Decision:** Each trace has independent `meta.interactions` schema
**Rationale:** Multi-trace figures need different interaction semantics per trace

### 4. Debounced Highlighting

**Decision:** 100ms debounce for plot highlighting
**Rationale:** Prevents excessive DOM updates when multiple state changes happen rapidly

## Limitations & Future Work

### Current Limitations

1. Only first trace is processed for interactions (multi-trace figures use first trace schema)
2. Hover only supports single item (not multiple simultaneous hovers)
3. Geometry selection highlighting requires customdata to directly match selection indices

### Potential Enhancements

1. Support multiple traces with independent interaction schemas
2. Add reverse index (geometry ID → customdata values) for non-aligned data
3. Support multi-selection hover with visual feedback
4. Add animation/transition effects for highlighting
5. Add keyboard modifiers support (Ctrl+Click for multi-select)

## Files Created

- `/Users/fzills/tools/zndraw-communication-testing/examples/plotly_interactions_demo.py` - Example patterns
- `/Users/fzills/tools/zndraw-communication-testing/docs/PLOTLY_INTERACTIONS.md` - User guide
- `/Users/fzills/tools/zndraw-communication-testing/IMPLEMENTATION_SUMMARY.md` - This file

## Integration Checklist

- ✅ Backend validation function implemented
- ✅ Frontend event bindings implemented
- ✅ State synchronization (highlighting) implemented
- ✅ Zustand store extended with hoveredFrame
- ✅ Type definitions for all interfaces
- ✅ Error handling and validation
- ✅ Performance optimization (debouncing, O(1) lookups)
- ✅ Documentation (docstrings, guide, examples)
- ✅ Testing (validation function verified)
- ✅ TypeScript compilation successful

## Conclusion

The Plotly interactions feature is fully implemented according to the design document. Users can now create rich, interactive visualizations that synchronize bidirectionally with the 3D geometry view. The implementation prioritizes flexibility (user-controlled schemas), maintainability (clear error messages), and performance (efficient lookups and debouncing).
