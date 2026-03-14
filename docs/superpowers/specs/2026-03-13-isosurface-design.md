# Isosurface Geometry

Server-side marching cubes extraction for volumetric data (molecular orbitals, electron density, electrostatic potentials) rendered as translucent 3D meshes.

## Architecture

Three components:

1. **Pydantic model** (`Isosurface`) — inherits `BaseModel` directly (not `BaseGeometry`), registered in the `geometries` dict. Stored via existing geometry CRUD endpoints.
2. **REST endpoint** — `GET /v1/rooms/{room_id}/frames/{index}/isosurface` reads volumetric data from a frame, runs marching cubes, returns mesh as msgpack.
3. **React component** (`Isosurface.tsx`) — fetches mesh from endpoint, renders `BufferGeometry` with `MeshPhysicalMaterial`.

## Data Flow

```
ase.Atoms.info["orbital_homo"] = {
    "grid": np.ndarray (nx, ny, nz),    # scalar field
    "origin": np.ndarray (3,),           # world-space origin
    "cell": np.ndarray (3, 3),           # axis vectors
}
     │
     │  asebytes.encode() → b"info.orbital_homo" (single msgpack key)
     │
     ▼
GET /v1/rooms/{room_id}/frames/{index}/isosurface?cube_key=info.orbital_homo&isovalue=0.02
     │
     │  msgpack_numpy.unpackb(frame[b"info.orbital_homo"])
     │  → {"grid": ndarray, "origin": ndarray, "cell": ndarray}
     │
     │  marching_cubes(grid, level=isovalue, step_size=step_size)
     │  transform verts: origin + (verts / (shape - 1)) @ cell
     │
     ▼
msgpack { "vertices": float32 (N, 3), "faces": uint32 (M, 3) }
     │
     ▼
Isosurface.tsx → BufferGeometry → MeshPhysicalMaterial
```

## Design Decisions

### No BaseGeometry inheritance

Isosurface has no per-instance positions, no instancing, no material system, no selection/hover interactions. It inherits from `BaseModel` with `frozen=True` to match the existing geometry convention.

Fields inherited from `BaseGeometry` that do NOT apply: `position`, `material`, `selection`, `selecting`, `hovering`. Including them would be dead weight.

The `owner` and `active` fields ARE included since they're needed for the geometry CRUD system (ownership checks, visibility toggle). `cube_key` has `default=""` so that `_get_geometry_types_info()` can instantiate the model with zero arguments to produce defaults.

### Single-surface model

Each `Isosurface` instance represents ONE surface at ONE isovalue. The isovalue can be negative. For molecular orbitals (positive + negative lobes), create two geometries:

```python
vis.geometries["homo+"] = Isosurface(
    cube_key="info.orbital_homo", isovalue=0.02, color="#2244CC",
)
vis.geometries["homo-"] = Isosurface(
    cube_key="info.orbital_homo", isovalue=-0.02, color="#CC4422",
)
```

This avoids baking orbital-specific concepts (positive/negative lobes, dual colors) into the model. Works equally well for electron density, electrostatic potential, or any scalar field.

### Server-side marching cubes

Keeps the frontend thin. scikit-image's `marching_cubes` is fast (~5-50ms depending on grid size) and well-maintained.

### Frontend must NOT load cube data through the frame batch system

This is a critical constraint. Volumetric grids can be tens of megabytes — loading them into the browser through the normal frame prefetch pipeline would be catastrophic for performance.

The Isosurface component enforces this by:
1. **Not calling `useRegisterFrameKeys()`** — the batch prefetcher only fetches keys that geometry components explicitly register. Isosurface registers none.
2. **Not calling `getFrameBatched()`** — it fetches exclusively from the dedicated `/isosurface` endpoint, which returns only the extracted mesh (kilobytes, not megabytes).
3. **The `cube_key` field uses `x-custom-type: "dynamic-enum"` with `"dynamic-atom-props"`** — the dropdown is populated from frame metadata (via `useAvailableProperties`), which includes both per-atom keys (`arrays.*`) and global keys (`info.*`). The `info.orbital_homo` key will appear in the dropdown. Selecting it only stores the string value — the `DynamicEnumRenderer` never auto-fetches frame data (`enabled: false`).

The raw volumetric data stays server-side at all times. Only the extracted mesh geometry crosses the wire.

### No server-side caching

React Query handles client-side deduplication. Same `(roomId, frame, cube_key, isovalue)` tuple won't re-fetch. Optimize server-side later if profiling shows need.

### Cube data in `atoms.info` (not `atoms.arrays`)

`atoms.arrays` is per-atom (shape must be `(n_atoms, ...)`). Volumetric grids are per-frame data. `atoms.info` supports arbitrary dicts, and `asebytes.encode()` serializes the entire dict as a single msgpack key `b"info.<name>"`.

Verified: `atoms.info["orbital_homo"] = {"grid": ..., "origin": ..., "cell": ...}` encodes to `b"info.orbital_homo"` containing all three arrays.

## Pydantic Model

```python
# src/zndraw/geometries/isosurface.py
"""Isosurface geometry for volumetric data visualization."""

from pydantic import BaseModel, ConfigDict, Field


class Isosurface(BaseModel):
    """An isosurface extracted from a 3D volumetric grid.

    References a key in ``atoms.info`` that contains a dict with:
    - ``grid``: 3D float array (Nx, Ny, Nz) of scalar values
    - ``origin``: 3-vector, world-space origin of the grid
    - ``cell``: (3, 3) matrix, axis vectors spanning the grid

    The frontend fetches mesh data from a dedicated endpoint that runs
    marching cubes server-side.

    Example
    -------
    >>> from zndraw.geometries import Isosurface
    >>> iso = Isosurface(cube_key="info.orbital_homo", isovalue=0.02)
    """

    model_config = ConfigDict(frozen=True)

    owner: str | None = Field(
        default=None,
        description="User ID of the geometry owner. None means unowned.",
    )

    cube_key: str = Field(
        default="",
        description=(
            "Frame info key for the volumetric data dict. "
            "Must contain grid, origin, and cell entries."
        ),
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props"],
        },
    )

    isovalue: float = Field(
        default=0.02,
        description="Scalar threshold for surface extraction. Can be negative.",
        json_schema_extra={"format": "range", "step": 0.001, "min": -1.0, "max": 1.0},
    )

    color: str = Field(
        default="#2244CC",
        description="Surface color.",
        json_schema_extra={
            "x-custom-type": "color-picker",
            "x-features": ["color-picker"],
        },
    )

    resolution: int = Field(
        default=1,
        ge=1,
        le=8,
        description=(
            "Marching cubes step size. 1 = full resolution, "
            "2 = half, etc. Higher values give faster extraction."
        ),
        json_schema_extra={"format": "range", "step": 1},
    )

    opacity: float = Field(
        default=0.6,
        ge=0.0,
        le=1.0,
        description="Surface opacity.",
        json_schema_extra={"format": "range", "step": 0.01},
    )

    active: bool = Field(default=True, description="Show or hide this isosurface.")
```

### Registration

In `src/zndraw/geometries/__init__.py`:
- Import `Isosurface`
- Add `"Isosurface": Isosurface` to `geometries` dict
- Call `Isosurface.model_rebuild()` (for consistency, even though it has no forward refs)
- Add to `__all__`

## REST Endpoint

### `GET /v1/rooms/{room_id}/frames/{index}/isosurface`

**Query parameters:**
- `cube_key` (str, required) — frame key for the volumetric data dict
- `isovalue` (float, required) — scalar threshold for surface extraction
- `step_size` (int, default 1) — marching cubes step size (1 = full resolution, 2 = half, etc.)

**Response:** `application/x-msgpack` containing:
- `vertices`: float32 array (N, 3) — world-space vertex positions
- `faces`: uint32 array (M, 3) — triangle face indices

**Errors (RFC 9457):**
- `404 FrameNotFound` — frame index out of range
- `404 RoomNotFound` — room does not exist
- `422 UnprocessableContent` — cube_key missing from frame, or grid not 3D, or dict missing required keys

**Dependencies:** `SessionMakerDep`, `StorageDep`, `CurrentUserFactoryDep` (for auth).

**Implementation:**

```python
# src/zndraw/routes/isosurface.py

router = APIRouter(
    prefix="/v1/rooms/{room_id}/frames/{index}/isosurface",
    tags=["isosurface"],
)

@router.get("", response_class=MessagePackResponse, ...)
async def get_isosurface(
    session_maker: SessionMakerDep,
    storage: StorageDep,
    current_user: CurrentUserFactoryDep,
    room_id: str,
    index: int,
    cube_key: Annotated[str, Query(...)],
    isovalue: Annotated[float, Query(...)],
    step_size: Annotated[int, Query(ge=1, le=8)] = 1,
) -> Response:
    # 1. Verify room exists
    # 2. Check frame index bounds
    # 3. Read frame, get cube_key from frame
    # 4. Unpack dict → {"grid": ndarray, "origin": ndarray, "cell": ndarray}
    # 5. Validate dict has required keys and grid is 3D
    # 6. Run marching_cubes(grid, level=isovalue, step_size=step_size)
    # 7. Transform vertices: origin + (verts / (shape - 1)) @ cell
    #    (note: step_size changes the vertex scaling — verts are in
    #     step-adjusted index space, so divide by (shape - 1) still works
    #     because marching_cubes already accounts for step_size internally)
    # 8. Return msgpack {vertices, faces}
```

### Vertex transformation

```python
shape = np.array(grid.shape, dtype=np.float64)
fractional = verts / (shape - 1)
world = origin + fractional @ cell
```

`marching_cubes` returns vertices in grid-index space `[0, Nx-1]`. Dividing by `(shape - 1)` gives fractional coordinates `[0, 1]`. Multiplying by the cell matrix and adding the origin gives world coordinates.

### Empty surface handling

If `marching_cubes` raises `ValueError` (no surface at this isovalue) or returns zero vertices, return msgpack `{"vertices": [], "faces": []}` with status 200. Uniform response shape means the frontend doesn't need to distinguish empty from populated.

## Frontend Component

```tsx
// frontend/src/components/three/Isosurface.tsx

export default function Isosurface({
  data,
  geometryKey,
}: {
  data: IsosurfaceData;
  geometryKey: string;
}) {
  // 1. Merge with defaults via getGeometryWithDefaults()
  // 2. Read currentFrame, roomId from store
  // 3. useQuery to fetch from isosurface endpoint
  //    queryKey: ["isosurface", roomId, currentFrame, cube_key, isovalue, resolution]
  // 4. Build BufferGeometry from vertices + faces
  // 5. computeVertexNormals()
  // 6. Render <mesh> with <meshPhysicalMaterial>
}
```

### Component registration

In `Canvas.tsx`, add `Isosurface` to `SIMPLE_GEOMETRY_COMPONENTS` (receives `geometryKey` + `data`, no `pathtracingEnabled`):

```typescript
const SIMPLE_GEOMETRY_COMPONENTS = {
  Camera, DirectionalLight, AmbientLight, HemisphereLight, Fog, Isosurface,
} as const;
```

Isosurface does not use instanced meshes or the pathtracing renderer, making `SIMPLE_GEOMETRY_COMPONENTS` the correct category.

### Material properties

- `MeshPhysicalMaterial` with `side={THREE.DoubleSide}` (surfaces visible from both sides)
- `transparent={opacity < 1.0}`, `depthWrite={opacity >= 1.0}`
- `roughness={0.4}` for a slightly glossy scientific-visualization look

## Dependencies

- **Add:** `scikit-image` (for `skimage.measure.marching_cubes`)
- **Add (dev):** `pyscf` (for integration tests with real orbital data on small molecules)
- **Already present:** `msgpack-numpy`, `msgpack`, `numpy`

## Testing

**Implementation must use the `test-driven-development` skill** — write tests first, then implement.

### Unit tests (no server needed)

1. `test_extract_mesh_sphere` — create a sphere SDF grid, run `_extract_mesh`, verify vertices are within grid bounds and faces are valid triangles
2. `test_extract_mesh_no_surface` — uniform grid with isovalue outside range, verify returns `None`
3. `test_extract_mesh_negative_isovalue` — grid with negative values, extract at negative isovalue

### Integration tests (async, full server)

Using existing conftest fixtures (`auth_client`, `room_id`, frame storage):

1. `test_isosurface_basic` — store cube data in a frame, GET isosurface, verify 200 with vertices/faces
2. `test_isosurface_missing_cube_key` — GET with nonexistent key, verify 422
3. `test_isosurface_frame_not_found` — GET with out-of-range frame index, verify 404
4. `test_isosurface_empty_surface` — isovalue outside data range, verify 200 with empty vertices/faces
5. `test_isosurface_invalid_grid` — cube data with non-3D grid, verify 422
6. `test_isosurface_missing_dict_keys` — cube data dict missing "grid" key, verify 422
7. `test_isosurface_step_size` — verify coarser resolution produces fewer vertices
8. `test_isosurface_pyscf_h2` — use PySCF to generate a real H2 HOMO orbital grid, store in frame, extract isosurface, verify mesh is non-empty and vertices are in plausible coordinate range

## Files

| File | Action | Description |
|------|--------|-------------|
| `pyproject.toml` | Modify | Add `scikit-image` dependency |
| `src/zndraw/geometries/isosurface.py` | Create | Pydantic model |
| `src/zndraw/geometries/__init__.py` | Modify | Register in geometries dict |
| `src/zndraw/routes/isosurface.py` | Create | REST endpoint with marching cubes |
| `src/zndraw/app.py` | Modify | Include isosurface router |
| `tests/test_isosurface.py` | Create | Unit + integration tests |
| `frontend/src/components/three/Isosurface.tsx` | Create | React/three.js component |
| `frontend/src/components/Canvas.tsx` | Modify | Register in geometry dispatch |

## Python Client Usage

```python
import ase
import numpy as np
from pyscf.tools import cubegen

# PySCF calculation
mol = ...
mf = ...
homo_idx = ...

# Generate orbital grid
orb_on_grid = cubegen.orbital(mol, "homo.cube", mf.mo_coeff[:, homo_idx])
cc = cubegen.Cube(mol)

# Build atoms with volumetric data
atoms = ase.Atoms(...)
atoms.info["orbital_homo"] = {
    "grid": orb_on_grid,           # (nx, ny, nz) float
    "origin": cc.boxorig * 0.529,  # Bohr → Angstrom
    "cell": cc.box * 0.529,        # Bohr → Angstrom
}

# Visualize
from zndraw import ZnDraw
from zndraw.geometries import Isosurface

vis = ZnDraw(url="http://localhost:1234")
vis.append(atoms)

vis.geometries["homo+"] = Isosurface(
    cube_key="info.orbital_homo",
    isovalue=0.02,
    color="#0044FF",
    opacity=0.6,
)
vis.geometries["homo-"] = Isosurface(
    cube_key="info.orbital_homo",
    isovalue=-0.02,
    color="#FF4400",
    opacity=0.6,
)
```
