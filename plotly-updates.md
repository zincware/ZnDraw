# Plotly Figure Interactions - Design & Implementation Plan

## Problem Statement

Currently, Plotly figures in ZnDraw are static visualizations. We need bidirectional synchronization between plots and 3D geometry:
- Clicking/selecting points in plots → updates frame, particle selection, or other state
- State changes (frame, selection) → highlight corresponding points in plots

The challenge: Different plots have different data structures and interaction semantics. A fixed `[step, particle]` assumption breaks for many common plots (energy predictions, force distributions, dimensionality reductions).

## Proposed Solution: Explicit Interaction Schema via `meta.interactions`

### Core Design Principles

1. **User controls `customdata` completely** - no assumptions or helper functions
2. **Explicit interaction contract via `meta.interactions`** - array indexed by customdata dimension
3. **Frontend as dumb parser** - reads schema, dispatches state updates
4. **Works with any plot type** - scatter, histogram, heatmap, box, contour, etc.
5. **Graceful degradation** - plots without interactions still display normally

### Data Structure

```python
# User manually sets customdata with their desired structure
fig.update_traces(
    customdata=np.column_stack([step_indices, particle_indices]),
    meta={
        "interactions": [
            {"click": "step", "select": "step"},      # index 0: customdata[0]
            {"click": "particle", "select": "particle"} # index 1: customdata[1]
        ]
    }
)
```

### Schema Format: `meta.interactions`

Each element in the `interactions` array corresponds to the same dimension in `customdata`:

```python
# Define which actions trigger for each customdata dimension
meta={
    "interactions": [
        {"click": "step", "select": "step"},              # dimension 0
        {"click": "particle", "select": None},            # dimension 1 - select disabled
        None,                                             # dimension 2 - no interaction
        {"select": "forces"}                              # dimension 3 - select only
    ]
}
```

### Standard Action Names

Frontend recognizes these action values and their semantics differ by event type:

#### Special Action: `"step"`
- **Click action (`"click": "step"`)**: Jump to a specific frame (sets `currentFrame` in zustand to the exact value)
- **Select action (`"select": "step"`)**: Add frames to the frame selection (updates `frame_selection` in zustand with all customdata values from the selected points)
- **Hover action (`"hover": "step"`)**: Highlight current frame in progressbar (include hovered frame in store.tsx)

#### Geometry Actions: `"<geometry_name>"`
Where `<geometry_name>` matches a key in `vis.geometries` (e.g., `"particles"`, `"forces"`, `"bonds"`, etc.)

- **Click action** (`"click": "<geometry_name>"`): Exclusively select this single item in the geometry (calls `setSelection(geometryName, [customdata_value])`)
- **Select action** (`"select": "<geometry_name>"`): Exclusively select all items matching the lasso/box selection (calls `setSelection(geometryName, [all_customdata_values])`)
- **Hover action** (`"hover": "<geometry_name>"`): Highlight the hovered item (updates `hoveredGeometryInstance` in zustand; will be extended to support multi-geometry)

**Key semantics:**
- **All selections are exclusive for that geometry** – clicking/selecting a geometry replaces any previous selection in that geometry
- **Multiple geometries can be selected simultaneously** – each geometry maintains its own independent selection via `selections[<geometry_name>]`
- **Single point vs. region selection** – click on a point sets `[value]`, lasso/box selection sets `[...all_values]`
- **Deselect** – only affects the specified geometry, others remain unchanged

### Examples

#### Example 1: Energy vs Accuracy (step mapping only)
```python
import plotly.express as px
import numpy as np

pred_energy = np.random.rand(10)
true_energy = np.random.rand(10)
frame_indices = np.arange(10)

fig = px.scatter(x=pred_energy, y=true_energy)
fig.update_traces(
    customdata=frame_indices,
    meta={
        "interactions": [
            {"click": "step", "select": "step"}
        ]
    }
)
# Clicking a point sets current frame
# Selecting region sets frame_selection
```

#### Example 2: Force Distribution (forces geometry mapping)
```python
forces = np.random.rand(100)
force_ids = np.arange(100)  # IDs of force vectors in vis.geometries["forces"]

fig = px.histogram(x=forces)
fig.update_traces(
    customdata=force_ids,  # 1D array
    meta={
        "interactions": [
            {"click": "forces", "select": "forces"}  # maps to vis.geometries["forces"] selection
        ]
    }
)
# Clicking a bin selects all forces in that range
# This works because we explicitly map customdata[0] to the "forces" geometry
```

#### Example 3: Multi-dimensional (step + particles geometry)
```python
pred_forces = np.random.rand(50)
particle_ids = np.tile(np.arange(10), 5)  # 10 particles, 5 steps (indices into vis.geometries["particles"])
steps = np.repeat(np.arange(5), 10)       # repeated for each particle

fig = px.scatter(x=particle_ids, y=pred_forces)
fig.update_traces(
    customdata=np.column_stack([steps, particle_ids]),
    meta={
        "interactions": [
            {"click": "step", "select": "step"},
            {"click": "particles", "select": "particles"}  # maps to vis.geometries["particles"] selection
        ]
    }
)
# Clicking point selects both the current step (jumps to it) and the particle
# The selection updates vis.selections["particles"]
```

#### Example 4: Sparse interactions (nulls and selective actions)
```python
# Only step selection on click, particles on select
fig.update_traces(
    customdata=np.column_stack([steps, particle_ids, force_magnitudes]),
    meta={
        "interactions": [
            {"click": "step"},                    # click only
            {"select": "particles"},              # select only, maps to vis.geometries["particles"]
            None                                  # force_magnitudes dimension: no interaction
        ]
    }
)
```

#### Example 5: Multiple geometries
```python
# Map to different geometry types (e.g., vis.geometries["bonds"] and vis.geometries["atoms"])
fig.update_traces(
    customdata=np.column_stack([bond_ids, atom_ids]),
    meta={
        "interactions": [
            {"click": "bonds", "select": "bonds"},      # maps to vis.geometries["bonds"]
            {"click": "atoms", "select": "atoms"}       # maps to vis.geometries["atoms"]
        ]
    }
)
```

## Plotly.js Event System (Reference)

Plotly.js provides these events:
- `plotly_click` - user clicked a data point
- `plotly_selected` - user completed a selection (box/lasso)
- `plotly_selecting` - user is actively selecting (fires during drag)
- `plotly_deselect` - user deselected
- `plotly_hover` / `plotly_unhover` - mouse over/out

Event data structure for click:
```typescript
{
  points: [
    {
      curveNumber: 0,        // trace index
      pointNumber: 5,        // point index in trace
      x: 3.5,                // x value
      y: 12.4,               // y value
      customdata: [0, 5],    // User's customdata (what we use!)
      data: {...},           // Reference to trace data
      // ... other fields
    }
  ]
}
```

Binding events:
```typescript
plotDiv.on('plotly_click', (data) => {
  // Handle click
});

plotDiv.on('plotly_selected', (data) => {
  // Handle selection
});
```

## Implementation Plan

### Phase 1: Backend Documentation & Examples
**Goal**: Enable users to add interactions to figures

Tasks:
- [ ] Add docstring to `Figures` class explaining `meta.interactions` format
- [ ] Create example script showing 3 common patterns (step, particle, step+particle)
- [ ] Add validation helper function:
  ```python
  def validate_interaction_schema(customdata, interactions_meta):
      """Verify customdata dimensions match interactions array length and validate the meta schema, keys are valid and the geometries exist in vis.geometries."""
  ```

### Phase 2: Frontend Event Binding
**Goal**: Wire up Plotly events to state updates

Files to modify: `app/src/components/FigureWindow.tsx`

Tasks:
- [ ] Parse `meta.interactions` schema from figure
  ```typescript
  interface InteractionSchema {
    click?: string | null;
    select?: string | null;
    hover?: string | null;
  }

  function parseInteractions(figure: PlotlyData): InteractionSchema[] {
    return figure.meta?.interactions || [];
  }

  function validateSchema(customdata: any, schema: InteractionSchema[]): boolean {
    // Customdata can be 1D array or 2D array
    // For 2D, first element's length determines dimensionality
    if (!customdata) return true; // No customdata is fine
    const expectedDims = Array.isArray(customdata[0]) ? customdata[0].length : 1;
    if (schema.length > 0 && schema.length !== expectedDims) {
      console.error(
        `Schema mismatch: customdata has ${expectedDims} dimension(s), schema has ${schema.length}`
      );
      return false;
    }
    return true;
  }
  ```

- [ ] Implement batch action processor
  ```typescript
  interface ActionBatch {
    step?: number; // Click: single frame value
    frameSelection?: number[]; // Select: multiple frame values
    geometries: Record<string, Set<number>>; // geometry_name -> Set of customdata values (for O(1) lookup)
    hoveredFrame?: number | null;
    hoveredGeometry?: [string, number] | null; // [geometry_name, index]
  }

  function processBatch(batch: ActionBatch) {
    // Handle step click (replace current frame)
    if (batch.step !== undefined) {
      setFrameAtomic(batch.step);
    }

    // Handle frame selection (add frames)
    if (batch.frameSelection?.length) {
      batch.frameSelection.forEach(frameIdx => addFrameSelection(frameIdx));
    }

    // Handle geometry selections (exclusive per-geometry, replaces previous)
    Object.entries(batch.geometries).forEach(([geometryName, valuesSet]) => {
      if (geometries[geometryName]) {
        // Convert Set to array for setSelection
        setSelection(geometryName, Array.from(valuesSet));
      }
    });

    // Handle hover (visual feedback only)
    if (batch.hoveredFrame !== undefined) {
      setHoveredFrame(batch.hoveredFrame);
    }
    if (batch.hoveredGeometry !== undefined) {
      setHoveredGeometryInstance(batch.hoveredGeometry);
    }
  }

  function collectActions(
    points: PlotlyPoint[],
    schema: InteractionSchema[],
    eventType: 'click' | 'select' | 'hover'
  ): ActionBatch {
    const batch: ActionBatch = { geometries: {} };

    points.forEach(point => {
      // Handle both 1D and 2D customdata
      const customdata = Array.isArray(point.customdata?.[0])
        ? point.customdata
        : [point.customdata];

      customdata?.forEach((value, dimensionIndex) => {
        const interaction = schema[dimensionIndex];
        if (!interaction) return;

        if (eventType === 'click' && interaction.click) {
          if (interaction.click === 'step') {
            batch.step = value as number; // Last value (single point click)
          } else {
            // Geometry action: collect to Set
            if (!batch.geometries[interaction.click]) {
              batch.geometries[interaction.click] = new Set<number>();
            }
            batch.geometries[interaction.click].add(value as number);
          }
        } else if (eventType === 'select' && interaction.select) {
          if (interaction.select === 'step') {
            if (!batch.frameSelection) batch.frameSelection = [];
            batch.frameSelection.push(value as number);
          } else {
            // Geometry action: collect to Set
            if (!batch.geometries[interaction.select]) {
              batch.geometries[interaction.select] = new Set<number>();
            }
            batch.geometries[interaction.select].add(value as number);
          }
        } else if (eventType === 'hover' && interaction.hover) {
          if (interaction.hover === 'step') {
            batch.hoveredFrame = value as number;
          } else {
            // Geometry hover: only first value (single hover)
            batch.hoveredGeometry = [interaction.hover, value as number];
          }
        }
      });
    });

    return batch;
  }
  ```

- [ ] Bind `plotly_click` event
  ```typescript
  plotDiv.on('plotly_click', (data) => {
    const schema = parseInteractions(currentFigure);
    if (!validateSchema(data.points[0]?.customdata, schema)) {
      // Show error pill: "Interaction schema mismatch"
      return;
    }

    const batch = collectActions(data.points, schema, 'click');
    processBatch(batch);
  });
  ```

- [ ] Bind `plotly_selected` event
  ```typescript
  plotDiv.on('plotly_selected', (data) => {
    const schema = parseInteractions(currentFigure);
    if (!validateSchema(data.points[0]?.customdata, schema)) {
      // Show error pill: "Interaction schema mismatch"
      return;
    }

    const batch = collectActions(data.points, schema, 'select');
    processBatch(batch);
  });
  ```

- [ ] Bind `plotly_hover` event
  ```typescript
  plotDiv.on('plotly_hover', (data) => {
    const schema = parseInteractions(currentFigure);
    const batch = collectActions([data.points[0]], schema, 'hover');
    processBatch(batch);
  });
  ```

- [ ] Bind `plotly_unhover` and `plotly_deselect` events
  ```typescript
  plotDiv.on('plotly_unhover', () => {
    setHoveredFrame(null);
    setHoveredGeometryInstance(null);
  });

  plotDiv.on('plotly_deselect', () => {
    // Clear all geometry-specific selections for this figure
    const schema = parseInteractions(currentFigure);
    const geometryNames = new Set<string>();

    schema.forEach(interaction => {
      if (interaction?.select && interaction.select !== 'step') {
        geometryNames.add(interaction.select);
      }
    });

    geometryNames.forEach(name => setSelection(name, []));
    // Note: frame_selection is NOT cleared on deselect
  });
  ```

### Phase 3: State Synchronization (Plot → 3D)
**Goal**: Update plot highlighting when state changes

Files to modify: `app/src/components/FigureWindow.tsx`

**Key implementation notes:**
- **Consolidated effect**: All state changes trigger a single `useEffect` that recomputes all highlights at once (prevents race conditions)
- **Efficient lookup**: Build a `Map<customdata_value, point_index>` once at plot creation, enabling O(1) lookups during highlighting
- **Debounced updates**: Multiple state changes within a short time are batched into a single DOM update

Tasks:
- [ ] Build customdata ↔ point index lookup on plot creation
  ```typescript
  function buildCustomdataIndex(figure: PlotlyFigure, traceIndex: number = 0): Map<number, number[]> {
    // Maps customdata value -> array of point indices
    // Handles both 1D and 2D customdata
    const index = new Map<number, number[]>();

    const data = figure.data[traceIndex];
    if (!data?.customdata) return index;

    const customdata = data.customdata;
    customdata.forEach((value, pointIndex) => {
      // For 2D customdata, use first dimension only (step)
      const key = Array.isArray(value) ? value[0] : value;
      if (!index.has(key)) {
        index.set(key, []);
      }
      index.get(key)!.push(pointIndex);
    });

    return index;
  }

  // Create once and memoize
  const customdataIndex = useMemo(
    () => buildCustomdataIndex(currentFigure),
    [currentFigure]
  );
  ```

- [ ] Consolidate all state synchronization into single effect
  ```typescript
  useEffect(() => {
    if (!customdataIndex.size) return;

    const schema = parseInteractions(currentFigure);
    const subscribedActions = new Set<string>();

    schema.forEach(interaction => {
      if (interaction?.click) subscribedActions.add(interaction.click);
      if (interaction?.select) subscribedActions.add(interaction.select);
      if (interaction?.hover) subscribedActions.add(interaction.hover);
    });

    // Collect all points that should be highlighted
    const highlightedPoints = new Set<number>();

    // Check for "step" action subscribers
    if (subscribedActions.has('step')) {
      // Highlight current frame
      customdataIndex.get(currentFrame)?.forEach(idx => highlightedPoints.add(idx));
    }

    // Check for "frame_selection" action subscribers
    if (subscribedActions.has('step')) { // frame_selection uses "step" dimension
      frame_selection.forEach(frameIdx => {
        customdataIndex.get(frameIdx)?.forEach(idx => highlightedPoints.add(idx));
      });
    }

    // Check for geometry action subscribers
    Object.entries(selections).forEach(([geometryName, selectedIds]) => {
      if (subscribedActions.has(geometryName)) {
        // For each selected ID, find all points with that customdata value
        selectedIds.forEach(id => {
          // Find which dimension holds this geometry
          schema.forEach((interaction, dimIndex) => {
            if (interaction?.select === geometryName || interaction?.click === geometryName) {
              // This dimension maps to the geometry
              // Need to map selected ID back to customdata values...
              // (This requires reverse index from geometry ID -> customdata values)
            }
          });
        });
      }
    });

    updatePlotHighlightWithDebounce(highlightedPoints);
  }, [currentFrame, frame_selection, selections, currentFigure, customdataIndex]);
  ```

- [ ] Implement efficient highlighting with debouncing
  ```typescript
  const updatePlotHighlightWithDebounce = useMemo(
    () => debounce((highlightedPoints: Set<number>) => {
      if (!plotDiv.current) return;

      // Calculate opacities: O(n) where n = total points
      const opacities = Array.from({ length: data.length }, (_, i) =>
        highlightedPoints.has(i) ? 1.0 : 0.2
      );

      // Update plot efficiently
      Plotly.restyle(plotDiv.current, { 'marker.opacity': [opacities] }, [0]);
    }, 100), // 100ms debounce
    [plotDiv]
  );
  ```

**Critical notes on geometry selection highlighting:**
- Current design assumes geometry IDs in `selections[geometryName]` correspond directly to customdata values
- This requires users to ensure their customdata indices match the IDs they use in selections!

---

## Important Design Constraints

### Multi-Trace Figures
Each trace in a Plotly figure **must have its own independent `meta.interactions` schema**. The system assumes:
- One schema applies to one trace at a time
- When binding events, events are **trace-specific** (Plotly provides `curveNumber` for trace identification)
- If a figure has multiple traces, handle each trace's interactions independently

Example:
```python
# Figure with 2 traces, each with different interactions
fig.add_trace(go.Scatter(..., customdata=step_ids, name="Energy"))
fig.add_trace(go.Scatter(..., customdata=particle_ids, name="Forces"))

# Each trace needs its own schema in different figures
fig1_energy.data[0].meta = {"interactions": [{"click": "step"}]}
fig2_forces.data[0].meta = {"interactions": [{"click": "forces"}]}
```

### Customdata-to-Selection Alignment
When using geometry actions (e.g., `"click": "forces"`), the **customdata values must be valid indices into the geometry's selection array**:

```python
# BAD: customdata uses arbitrary IDs that don't match selection indices
force_ids = [999, 1005, 1234]  # Arbitrary IDs
fig.update_traces(customdata=force_ids, meta={"interactions": [{"click": "forces"}]})
# When user clicks, it tries to set selections["forces"] = [999]
# But vis.geometries["forces"] might only have 10 items!
# THIS SHOULD BE VALIDATED when setting the figure but also we need good error handling in frontend

# GOOD: customdata uses valid indices
force_indices = [0, 1, 2]  # Direct indices into geometry
fig.update_traces(customdata=force_indices, meta={"interactions": [{"click": "forces"}]})
# When user clicks, it sets selections["forces"] = [0]
```

### Selection vs. Hover Exclusivity
- **Selections** are persistent and exclusive per-geometry (one selection state at a time)
- **Hover** is ephemeral and single-value only (cannot hover multiple items simultaneously in current design)
- Selections and hovers operate independently (can have selection + hover on different items)

### Frame Selection Behavior
- **`select: "step"` on lasso** adds frames to `frame_selection` (cumulative)
- **`click: "step"` on point** replaces `currentFrame` (immediate jump)
- Deselecting points does NOT clear `frame_selection` – user must explicitly clear via UI or deselect all

---

## Validation & Error Handling

### Backend Validation (Python)
The backend validation function should enforce:

```python
def validate_interaction_schema(customdata, interactions_meta, vis_geometries):
    """
    Validate that meta.interactions aligns with data and available geometries.

    Parameters:
    - customdata: The array passed to fig.update_traces()
    - interactions_meta: The schema from meta={"interactions": [...]}
    - vis_geometries: Available geometries from vis.geometries

    Returns:
    - (is_valid, errors_list)

    Checks:
    1. customdata and schema have matching dimensions
    2. All action names are either "step" or valid geometry names
    3. No invalid geometry references
    4. Schema is not empty or None
    """
    errors = []

    # Check 1: Dimension alignment
    if not customdata:
        errors.append("customdata is empty")
    else:
        expected_dims = len(customdata[0]) if isinstance(customdata[0], (list, tuple)) else 1
        if len(interactions_meta) != expected_dims:
            errors.append(
                f"Dimension mismatch: customdata has {expected_dims} dimension(s), "
                f"schema has {len(interactions_meta)}"
            )

    # Check 2 & 3: Valid action names
    for idx, interaction in enumerate(interactions_meta):
        if interaction is None:
            continue
        for action_type in ['click', 'select', 'hover']:
            action_name = interaction.get(action_type)
            if action_name is None:
                continue
            if action_name == 'step':
                continue  # Valid
            if action_name not in vis_geometries:
                errors.append(
                    f"Schema dimension {idx}, action '{action_type}': "
                    f"'{action_name}' not found in vis.geometries. "
                    f"Available: {list(vis_geometries.keys())}"
                )

    # Check 4: Schema not empty
    if not interactions_meta or all(x is None for x in interactions_meta):
        errors.append("Interaction schema is completely empty (all None)")

    return len(errors) == 0, errors
```

### Frontend Error Display
When validation fails, display an error pill in the figure:

```typescript
interface FigureState {
  validationErrors: string[];
  hasInteractionErrors: boolean;
}

// In FigureWindow component
const [figureState, setFigureState] = useState<FigureState>({
  validationErrors: [],
  hasInteractionErrors: false
});

// After parsing figure
useEffect(() => {
  const schema = parseInteractions(currentFigure);
  const isValid = validateSchema(currentFigure.data[0]?.customdata, schema);

  if (!isValid) {
    setFigureState({
      validationErrors: ['Interaction schema mismatch with customdata'],
      hasInteractionErrors: true
    });
  }
}, [currentFigure]);

// Render error pill at top of figure
{figureState.hasInteractionErrors && (
  <Alert severity="error" sx={{ mb: 1 }}>
    {figureState.validationErrors.map((err, i) => (
      <div key={i}>{err}</div>
    ))}
  </Alert>
)}
```

### Graceful Degradation
- If validation fails, the figure still renders (no crash)
- Error pill clearly indicates the problem
- User interactions are ignored (no state updates)
- Figure can be corrected by user and re-submitted

claude mcp add context7 -- bunx -y @upstash/context7-mcp --api-key ctx7sk-a067061d-9ce5-4063-8aae-b4367e5d9d6a