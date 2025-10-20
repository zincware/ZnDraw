# Plotly Figure Interactions

ZnDraw now supports **bidirectional synchronization** between interactive Plotly plots and 3D geometry:

- **Plot → 3D**: Clicking/selecting points in plots updates frame, particle selection, or other state
- **3D → Plot**: State changes (frame, selection) highlight corresponding points in plots

This feature allows users to create rich, interactive visualizations where plots and 3D geometry are tightly coupled.

## Quick Start

### Basic Example: Frame Selection

```python
import numpy as np
import plotly.express as px
from zndraw import ZnDraw

vis = ZnDraw()

# Create a scatter plot
fig = px.scatter(x=energy_predicted, y=energy_true)

# Add interaction schema
frame_indices = np.arange(len(energy_predicted))
fig.update_traces(
    customdata=frame_indices,
    meta={
        "interactions": [
            {"click": "step", "select": "step"}
        ]
    }
)

# Add to visualization
vis.figures["energy"] = fig

# Now clicking points in the plot jumps to that frame!
```

## Core Concepts

### Schema Format: `meta.interactions`

Each element in the `interactions` array corresponds to a dimension in `customdata`:

```python
fig.update_traces(
    customdata=np.column_stack([steps, particle_ids, force_mags]),
    meta={
        "interactions": [
            {"click": "step", "select": "step"},        # dimension 0
            {"click": "particles", "select": "particles"},  # dimension 1
            None                                        # dimension 2 (no interaction)
        ]
    }
)
```

### Standard Action Names

#### Special Action: `"step"`
Controls the current frame or frame selection:
- **`click: "step"`** - Jump to a specific frame (sets `currentFrame`)
- **`select: "step"`** - Add frames to selection (sets `frame_selection`)
- **`hover: "step"`** - Visual feedback for hovered frame

#### Geometry Actions: `"<geometry_name>"`
Maps to a geometry in `vis.geometries` (e.g., `"particles"`, `"forces"`, `"bonds"`):
- **`click: "<geometry_name>"`** - Exclusively select this item
- **`select: "<geometry_name>"`** - Exclusively select items matching box/lasso selection
- **`hover: "<geometry_name>"`** - Hover feedback

**Key semantics:**
- All selections are **exclusive per-geometry** (each geometry has one selection)
- Multiple geometries can be selected simultaneously (each maintains independent selection)
- Single click on point sets `[value]`, lasso/box selection sets `[...all_values]`

## Examples

### Example 1: Energy Prediction Plot

Simple 1D mapping for frame selection:

```python
import numpy as np
import plotly.express as px

n_frames = 100
pred_energy = np.random.rand(n_frames)
true_energy = np.random.rand(n_frames)
frame_indices = np.arange(n_frames)

fig = px.scatter(x=pred_energy, y=true_energy)
fig.update_traces(
    customdata=frame_indices,
    meta={
        "interactions": [
            {"click": "step", "select": "step"}
        ]
    }
)
```

**Result:**
- Clicking a point → jumps to that frame
- Selecting a region → adds frames to frame_selection
- Frame highlighting in plot when frame changes in 3D view

### Example 2: Force Distribution

Map histogram bins to a geometry:

```python
import numpy as np
import plotly.express as px

n_forces = 100
forces = np.random.exponential(scale=0.5, size=n_forces)
force_ids = np.arange(n_forces)

fig = px.histogram(x=forces, nbins=10)
fig.update_traces(
    customdata=force_ids,
    meta={
        "interactions": [
            {"click": "forces", "select": "forces"}
        ]
    }
)
```

**Result:**
- Clicking a bin → selects all forces in that range
- Updates `vis.selections["forces"]` (highlights forces in 3D)
- Force highlighting in plot when selection changes

### Example 3: Particle Trajectories (Multi-Dimensional)

Combine frame and geometry selection:

```python
import numpy as np
import plotly.express as px

n_steps = 10
n_particles = 50
steps = np.repeat(np.arange(n_steps), n_particles)
particle_ids = np.tile(np.arange(n_particles), n_steps)
forces = np.random.rand(n_steps * n_particles)

fig = px.scatter(x=particle_ids, y=forces)
fig.update_traces(
    customdata=np.column_stack([steps, particle_ids]),
    meta={
        "interactions": [
            {"click": "step", "select": "step"},
            {"click": "particles", "select": "particles"}
        ]
    }
)
```

**Result:**
- Clicking a point → jumps to frame AND selects particle
- Selecting region → adds frames AND selects particles
- Both frame and particle highlighting in plot

### Example 4: Sparse Interactions

Use `None` to disable dimensions and mix action types:

```python
fig.update_traces(
    customdata=np.column_stack([steps, particle_ids, force_mags]),
    meta={
        "interactions": [
            {"click": "step"},               # click only
            {"select": "particles"},         # select only
            None                             # no interaction
        ]
    }
)
```

**Result:**
- Clicking: only sets frame, doesn't select particles
- Selecting: only selects particles, doesn't add to frame_selection
- Force magnitude data not used for interactions

## Validation

Use `validate_interaction_schema()` to check your configuration:

```python
from zndraw.figures_manager import validate_interaction_schema

customdata = np.column_stack([steps, particle_ids])
interactions = [{"click": "step"}, {"click": "particles"}]

is_valid, errors = validate_interaction_schema(
    customdata,
    interactions,
    vis.geometries  # optional: validates geometry names exist
)

if not is_valid:
    for error in errors:
        print(f"Error: {error}")
```

### Common Validation Errors

1. **Dimension mismatch:**
   ```
   customdata has 2 dimension(s), schema has 1
   ```
   Solution: Add schema entries for each customdata dimension

2. **Invalid geometry:**
   ```
   'unknown_geo' not found in vis.geometries. Available: ['particles', 'forces']
   ```
   Solution: Use valid geometry names or create the geometry first

3. **Schema type error:**
   ```
   Schema dimension 0: interaction must be a dict or None
   ```
   Solution: Ensure each element is either a dict or None

## Important Design Constraints

### Multi-Trace Figures

Each trace needs its own `meta.interactions` schema. The system assumes one schema per trace:

```python
# Figure with 2 traces, each with different interactions
fig.add_trace(go.Scatter(..., customdata=step_ids, name="Energy"))
fig.add_trace(go.Scatter(..., customdata=particle_ids, name="Forces"))

# Set schema for first trace only (current implementation)
fig.data[0].meta = {"interactions": [{"click": "step"}]}
```

### Customdata-to-Selection Alignment

Customdata values must be valid indices into the geometry's selection:

```python
# GOOD: customdata uses valid indices
force_indices = [0, 1, 2, ...]  # Direct indices into vis.geometries["forces"]
fig.update_traces(
    customdata=force_indices,
    meta={"interactions": [{"click": "forces"}]}
)

# BAD: customdata uses arbitrary IDs
force_ids = [999, 1005, 1234]  # Won't match geometry selection!
fig.update_traces(
    customdata=force_ids,
    meta={"interactions": [{"click": "forces"}]}
)
```

### Selection vs. Hover Exclusivity

- **Selections**: Persistent, exclusive per-geometry (one selection at a time)
- **Hover**: Ephemeral, single-value only (cannot hover multiple items)
- Both operate independently (can have selection + hover on different items)

### Frame Selection Behavior

- **`select: "step"`** on lasso → adds frames to `frame_selection` (cumulative)
- **`click: "step"`** on point → replaces `currentFrame` (immediate jump)
- **Deselect** clears geometry selections but NOT `frame_selection`

## Frontend Error Handling

When validation fails:
- Figure still renders (no crash)
- Error alert displays problem
- User interactions are ignored
- Figure can be corrected and re-submitted

```
Interaction schema mismatch with customdata
```

## Performance Considerations

- **Efficient lookup**: O(1) customdata → point mapping using Map
- **Debounced updates**: Plot highlighting batched (100ms debounce)
- **O(n) highlighting**: n = number of points (unavoidable for per-point updates)

## Troubleshooting

### Plot not responding to clicks

**Check:**
1. Is `customdata` set on the trace?
2. Does `meta.interactions` exist?
3. Are dimensions matching? (run `validate_interaction_schema()`)

### Selections not highlighting in plot

**Check:**
1. Is the geometry name correct? (matches `vis.geometries` key)
2. Are customdata values valid indices? (0 to geometry_size-1)
3. Are state changes happening? (check console logs)

### Performance issues

**Solutions:**
1. Reduce number of points
2. Increase debounce delay
3. Use simpler plot types (scatter instead of scattergl)

## API Reference

### `validate_interaction_schema(customdata, interactions_meta, vis_geometries=None)`

Validates schema configuration.

**Parameters:**
- `customdata`: array passed to `fig.update_traces(customdata=...)`
- `interactions_meta`: schema from `meta={"interactions": [...]}`
- `vis_geometries`: optional dict of available geometries

**Returns:** `(is_valid: bool, errors: list[str])`

### Store Actions

**Frame Actions:**
- `setCurrentFrame(frame)` - Jump to frame
- `setFrameSelection(frames)` - Set frame selection
- `setHoveredFrame(frame)` - Set hovered frame

**Geometry Actions:**
- `updateSelectionForGeometry(geometry_name, indices)` - Set geometry selection
- `setHoveredGeometryInstance(geometry_name, instance_id)` - Set hovered instance

## Contributing

For issues or improvements related to interactions:
1. Check the design document: `plotly-updates.md`
2. Validate schema first: `validate_interaction_schema()`
3. Report with example code and error messages
