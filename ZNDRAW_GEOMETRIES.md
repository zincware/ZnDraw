# ZnDraw Geometries System

## Overview

ZnDraw is a 3D visualization tool for molecular dynamics and scientific data. The `vis.geometries` object provides a flexible way to render various geometric primitives (particles, arrows, etc.) with both static and dynamic data.

## Basic Usage

```python
from zndraw import ZnDraw

vis = ZnDraw()
vis.geometries = {
    "particles": {
        "position": "arrays.positions",  # Dynamic: fetch from frame data
        "color": "#FF0000",              # Static: hex color
        "radius": 1.0,                   # Static: single value
        "material": "standard"
    },
    "forces": {
        "position": "arrays.positions",
        "direction": "calc.forces",
        "color": "arrays.colors",        # Dynamic: fetch per-particle colors
        "radius": [1.0, 1.5, 2.0],      # Static: array of values
        "scale": 1.0,
        "material": "standard"
    }
}
```

## Data Types: Dynamic vs Static

### Dynamic Data (String Keys)
When a property is a **string** (except hex colors), it's treated as a **data key** that fetches values from the server for each frame:

```python
"position": "arrays.positions"     # Fetches position data from server
"direction": "calc.forces"         # Fetches force vectors from server
"color": "arrays.colors"           # Fetches per-particle RGB colors
```

**How it works:**
1. Frontend detects the string value
2. Creates a query to fetch data from `/api/frames/{frame_id}` with the key
3. Server returns msgpack-encoded numpy arrays
4. Frontend decodes and applies to geometry instances

### Static Data (Values)
When a property is **not a string** (or is a hex color), it's used directly:

```python
# Single value (replicated for all instances)
"radius": 1.0
"scale": 2.5
"color": "#FF0000"                 # Hex colors are NOT fetched

# Array of values (one per instance)
"radius": [1.0, 1.5, 2.0, 1.2]
"color": [[1, 0, 0], [0, 1, 0]]   # RGB values per particle

# Nested arrays
"position": [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
```

## Hex Colors (Special Case)

Hex color strings are **treated as static data** and NOT fetched from the server:

```python
"color": "#FF0000"      # ✅ Converted to RGB locally
"color": "#abc"         # ✅ Short form also works
"color": "arrays.colors" # ✅ Fetched from server (not a hex color)
```

**Implementation:** The frontend uses `shouldFetchAsFrameData()` to distinguish hex colors from data keys using regex pattern matching.

## Supported Geometries

### Particles (Spheres)
```python
{
    "position": "arrays.positions" | [[x,y,z], ...] | [x,y,z],
    "color": "data.key" | "#RRGGBB" | [[r,g,b], ...] | [r,g,b],
    "radius": "data.key" | [r1, r2, ...] | r,
    "material": "standard" | "basic" | "phong",
    "resolution": 8,  # sphere segments
    "scale": 1.0,
    "opacity": 1.0,
    "selecting": {"enabled": true, "color": "#FFC0CB", "opacity": 0.3},
    "hovering": {"enabled": true, "color": "#FFD700", "opacity": 0.2}
}
```

### Arrows
```python
{
    "position": "arrays.positions" | [[x,y,z], ...],
    "direction": "calc.forces" | [[dx,dy,dz], ...],
    "color": "data.key" | "#RRGGBB" | [[r,g,b], ...],
    "radius": "data.key" | [r1, r2, ...] | r,
    "scale": "data.key" | [s1, s2, ...] | s,
    "material": "standard"
}
```

## Data Flow

### Dynamic Data Flow
```
Python Server                Frontend
─────────────               ────────
vis.geometries
  └─ "arrays.positions"     → useQuery()
                              └─ GET /api/frames/{id}?keys=arrays.positions
                                  ← msgpack encoded numpy array

                              decode()
                              └─ TypedArray conversion

                              processedData
                              └─ Update THREE.js instances
```

### Static Data Flow
```
Python Server                Frontend
─────────────               ────────
vis.geometries
  └─ "#FF0000"              → hexToRgb()
                              └─ [1.0, 0.0, 0.0]

                              processedData
                              └─ Replicate for all instances
                              └─ Update THREE.js instances
```

## Performance Considerations

### Query Optimization
- Each dynamic property creates a separate React Query
- Only enabled when property is a string (not hex color)
- Cached by frame number and data key
- Fetching waits until ALL required queries complete

### Avoiding Re-fetches
- Hex colors are identified via regex to prevent unnecessary API calls
- `shouldFetchAsFrameData()` returns `false` for hex patterns
- Stale `fetchedData` is cleaned when props change from dynamic to static

### Update Strategy
1. Separate queries per property (position, direction, color, radius, scale)
2. Decode in individual `useEffect` hooks (avoids re-renders)
3. Combine in `useMemo` after all fetches complete
4. Update THREE.js instances in final `useEffect`

## Example: Switching Data Sources

```python
# Frame 0: Use server data
vis.geometries = {
    "arrows": {
        "position": "arrays.positions",
        "direction": "calc.forces",
        "color": "arrays.colors",
        "radius": 0.1,
        "scale": 1.0
    }
}

# Frame 1: Switch to static hex color
vis.geometries = {
    "arrows": {
        "position": "arrays.positions",
        "direction": "calc.forces",
        "color": "#FF0000",  # Changed to hex
        "radius": 0.1,
        "scale": 1.0
    }
}
```

**What happens:**
1. `keysToFetch` updates to exclude `color`
2. `fetchedData.color` is cleaned up (removed)
3. `processedData` uses `hexToRgb("#FF0000")` instead
4. Arrows render in red without fetching from server

## Type System

```typescript
type StaticValue = number | number[] | number[][];
type DataProp = string | StaticValue;

interface ParticleProps {
    position: DataProp;
    color: DataProp;     // Can be "#RRGGBB", "data.key", or values
    radius: DataProp;
    // ...
}
```

The type system allows seamless mixing of dynamic and static data sources.
