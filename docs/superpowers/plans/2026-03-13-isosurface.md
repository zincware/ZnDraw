# Isosurface Geometry Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add server-side marching cubes isosurface extraction for volumetric data, with a Pydantic model, REST endpoint, and React/Three.js rendering component.

**Architecture:** Isosurface is a geometry type that doesn't inherit `BaseGeometry` (no positions/instancing/materials). The frontend sends `(roomId, frame, cube_key, isovalue, resolution)` to a dedicated REST endpoint, which reads volumetric data from the frame, runs `skimage.measure.marching_cubes`, transforms vertices to world space, and returns a msgpack mesh. The React component renders the mesh as a `BufferGeometry` with `MeshPhysicalMaterial`.

**Tech Stack:** Python (FastAPI, scikit-image, msgpack-numpy, Pydantic), TypeScript (React, Three.js, @tanstack/react-query, @msgpack/msgpack)

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `pyproject.toml` | Modify | Add `scikit-image` dependency |
| `src/zndraw/geometries/isosurface.py` | Create | Pydantic model (frozen, no BaseGeometry) |
| `src/zndraw/geometries/__init__.py` | Modify | Register Isosurface in `geometries` dict |
| `src/zndraw/routes/isosurface.py` | Create | REST endpoint: marching cubes extraction |
| `src/zndraw/app.py` | Modify | Include isosurface router |
| `tests/test_isosurface.py` | Create | Unit + integration tests |
| `frontend/src/components/three/Isosurface.tsx` | Create | React/Three.js mesh renderer |
| `frontend/src/components/Canvas.tsx` | Modify | Register in `SIMPLE_GEOMETRY_COMPONENTS` |

---

## Chunk 1: Backend — Pydantic Model + Registration

### Task 1: Add scikit-image dependency

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Add scikit-image to dependencies**

In `pyproject.toml`, add `"scikit-image>=0.24.0"` to the `dependencies` list (after `"redis>=7.1.0"`):

```toml
    "scikit-image>=0.24.0",
```

- [ ] **Step 2: Sync dependencies**

Run: `uv sync`
Expected: Resolves and installs scikit-image successfully.

- [ ] **Step 3: Verify import**

Run: `uv run python -c "from skimage.measure import marching_cubes; print('OK')"`
Expected: `OK`

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml uv.lock
git commit -m "feat(isosurface): add scikit-image dependency"
```

---

### Task 2: Create Isosurface Pydantic model

**Files:**
- Create: `src/zndraw/geometries/isosurface.py`
- Test: `tests/test_isosurface.py` (model validation tests)

- [ ] **Step 1: Write model validation tests**

Create `tests/test_isosurface.py`:

```python
"""Tests for Isosurface geometry model and endpoint."""

import pytest


# =============================================================================
# Model Validation Tests
# =============================================================================


def test_isosurface_defaults():
    """Default construction produces expected field values."""
    from zndraw.geometries.isosurface import Isosurface

    iso = Isosurface()
    assert iso.cube_key == ""
    assert iso.isovalue == 0.02
    assert iso.color == "#2244CC"
    assert iso.resolution == 1
    assert iso.opacity == 0.6
    assert iso.active is True
    assert iso.owner is None


def test_isosurface_custom_values():
    """Construction with custom values works."""
    from zndraw.geometries.isosurface import Isosurface

    iso = Isosurface(
        cube_key="info.orbital_homo",
        isovalue=-0.05,
        color="#FF0000",
        resolution=2,
        opacity=0.8,
        owner="user123",
    )
    assert iso.cube_key == "info.orbital_homo"
    assert iso.isovalue == -0.05
    assert iso.color == "#FF0000"
    assert iso.resolution == 2
    assert iso.opacity == 0.8
    assert iso.owner == "user123"


def test_isosurface_frozen():
    """Model is immutable (frozen=True)."""
    from zndraw.geometries.isosurface import Isosurface

    iso = Isosurface()
    with pytest.raises(Exception):
        iso.isovalue = 0.5  # type: ignore[misc]


def test_isosurface_resolution_bounds():
    """Resolution must be between 1 and 8."""
    from pydantic import ValidationError

    from zndraw.geometries.isosurface import Isosurface

    with pytest.raises(ValidationError):
        Isosurface(resolution=0)
    with pytest.raises(ValidationError):
        Isosurface(resolution=9)


def test_isosurface_opacity_bounds():
    """Opacity must be between 0.0 and 1.0."""
    from pydantic import ValidationError

    from zndraw.geometries.isosurface import Isosurface

    with pytest.raises(ValidationError):
        Isosurface(opacity=-0.1)
    with pytest.raises(ValidationError):
        Isosurface(opacity=1.1)


def test_isosurface_schema_has_dynamic_enum():
    """cube_key field schema includes dynamic-enum metadata."""
    from zndraw.geometries.isosurface import Isosurface

    schema = Isosurface.model_json_schema()
    cube_key_props = schema["properties"]["cube_key"]
    assert cube_key_props["x-custom-type"] == "dynamic-enum"
    assert "dynamic-atom-props" in cube_key_props["x-features"]


def test_isosurface_in_geometries_dict():
    """Isosurface is registered in the geometries dict."""
    from zndraw.geometries import geometries

    assert "Isosurface" in geometries
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_isosurface.py -v`
Expected: All tests FAIL (module not found).

- [ ] **Step 3: Create the Isosurface model**

Create `src/zndraw/geometries/isosurface.py`:

```python
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

- [ ] **Step 4: Register in geometries dict**

Modify `src/zndraw/geometries/__init__.py`:

Add import at top (after the `from zndraw.geometries.fog import Fog` line):

```python
from zndraw.geometries.isosurface import Isosurface
```

Add `Isosurface.model_rebuild()` after the existing `model_rebuild()` calls (after `Shape.model_rebuild()`):

```python
Isosurface.model_rebuild()
```

Add to `geometries` dict (after the `"PropertyInspector": PropertyInspector,` line):

```python
    "Isosurface": Isosurface,
```

Add `"Isosurface"` to the `__all__` list (alphabetically, after `"InteractionSettings"`):

```python
    "Isosurface",
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `uv run pytest tests/test_isosurface.py -v`
Expected: All 7 tests PASS.

- [ ] **Step 6: Run linting**

Run: `uv run ruff check --select I --fix . && uv run ruff format .`
Expected: Clean.

- [ ] **Step 7: Commit**

```bash
git add src/zndraw/geometries/isosurface.py src/zndraw/geometries/__init__.py tests/test_isosurface.py
git commit -m "feat(isosurface): add Pydantic model and register in geometries dict"
```

---

## Chunk 2: Backend — Mesh Extraction Helper + Unit Tests

### Task 3: Write mesh extraction helper with unit tests

The extraction logic is a pure function that takes grid data and returns vertices/faces. This keeps the route handler thin and makes unit testing easy (no server needed).

**Files:**
- Modify: `src/zndraw/routes/isosurface.py` (create file)
- Modify: `tests/test_isosurface.py` (add unit tests)

- [ ] **Step 1: Write unit tests for mesh extraction**

Append to `tests/test_isosurface.py`:

```python
import numpy as np


# =============================================================================
# Mesh Extraction Unit Tests (no server needed)
# =============================================================================


def test_extract_mesh_sphere():
    """Extract mesh from a sphere SDF grid, verify vertices and faces."""
    from zndraw.routes.isosurface import _extract_mesh

    # Create a sphere SDF: negative inside, positive outside
    n = 30
    lin = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(lin, lin, lin, indexing="ij")
    grid = np.sqrt(x**2 + y**2 + z**2) - 0.5  # sphere radius 0.5

    origin = np.array([-1.0, -1.0, -1.0])
    cell = np.array([[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]])

    result = _extract_mesh(grid, origin, cell, isovalue=0.0, step_size=1)
    assert result is not None
    vertices, faces = result

    assert vertices.ndim == 2
    assert vertices.shape[1] == 3
    assert faces.ndim == 2
    assert faces.shape[1] == 3

    # Vertices should be in world space, roughly within the sphere
    assert np.all(vertices >= -1.5)
    assert np.all(vertices <= 1.5)

    # All face indices should be valid
    assert np.all(faces >= 0)
    assert np.all(faces < len(vertices))


def test_extract_mesh_no_surface():
    """Uniform grid with isovalue outside range returns None."""
    from zndraw.routes.isosurface import _extract_mesh

    grid = np.ones((10, 10, 10))
    origin = np.zeros(3)
    cell = np.eye(3)

    result = _extract_mesh(grid, origin, cell, isovalue=5.0, step_size=1)
    assert result is None


def test_extract_mesh_negative_isovalue():
    """Extract at negative isovalue from a grid with negative values."""
    from zndraw.routes.isosurface import _extract_mesh

    n = 20
    lin = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(lin, lin, lin, indexing="ij")
    # Field that crosses zero and has negative values
    grid = x + y + z

    origin = np.array([-1.0, -1.0, -1.0])
    cell = np.array([[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]])

    result = _extract_mesh(grid, origin, cell, isovalue=-0.5, step_size=1)
    assert result is not None
    vertices, faces = result
    assert len(vertices) > 0
    assert len(faces) > 0


def test_extract_mesh_step_size():
    """Coarser step_size produces fewer vertices."""
    from zndraw.routes.isosurface import _extract_mesh

    n = 40
    lin = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(lin, lin, lin, indexing="ij")
    grid = np.sqrt(x**2 + y**2 + z**2) - 0.5

    origin = np.zeros(3)
    cell = np.eye(3) * 2

    result_fine = _extract_mesh(grid, origin, cell, isovalue=0.0, step_size=1)
    result_coarse = _extract_mesh(grid, origin, cell, isovalue=0.0, step_size=2)

    assert result_fine is not None
    assert result_coarse is not None
    assert len(result_coarse[0]) < len(result_fine[0])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_isosurface.py::test_extract_mesh_sphere -v`
Expected: FAIL (import error — `_extract_mesh` doesn't exist yet).

- [ ] **Step 3: Implement the extraction helper**

Create `src/zndraw/routes/isosurface.py`:

```python
"""Isosurface REST endpoint — marching cubes extraction for volumetric data."""

from __future__ import annotations

import numpy as np


def _extract_mesh(
    grid: np.ndarray,
    origin: np.ndarray,
    cell: np.ndarray,
    isovalue: float,
    step_size: int,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Run marching cubes and transform vertices to world space.

    Parameters
    ----------
    grid : np.ndarray
        3D scalar field (Nx, Ny, Nz).
    origin : np.ndarray
        World-space origin (3,).
    cell : np.ndarray
        Axis vectors spanning the grid (3, 3).
    isovalue : float
        Scalar threshold for surface extraction.
    step_size : int
        Marching cubes step size (1 = full resolution).

    Returns
    -------
    tuple[np.ndarray, np.ndarray] | None
        (vertices, faces) as float32/uint32 arrays, or None if no surface.
    """
    from skimage.measure import marching_cubes

    try:
        verts, faces, _, _ = marching_cubes(grid, level=isovalue, step_size=step_size)
    except (ValueError, RuntimeError):
        return None

    if len(verts) == 0:
        return None

    shape = np.array(grid.shape, dtype=np.float64)
    fractional = verts / (shape - 1)
    world = origin + fractional @ cell

    return world.astype(np.float32), faces.astype(np.uint32)
```

- [ ] **Step 4: Run unit tests to verify they pass**

Run: `uv run pytest tests/test_isosurface.py -k "extract_mesh" -v`
Expected: All 4 extraction tests PASS.

- [ ] **Step 5: Commit**

```bash
git add src/zndraw/routes/isosurface.py tests/test_isosurface.py
git commit -m "feat(isosurface): add _extract_mesh helper with unit tests"
```

---

## Chunk 3: Backend — REST Endpoint + Integration Tests

### Task 4: Write integration tests for the isosurface endpoint

**Files:**
- Modify: `tests/test_isosurface.py` (add integration tests)

The integration tests follow the exact same fixture pattern as `tests/test_routes_frames.py`. Each test creates its own DB session, storage, user, room, stores cube data in a frame, then hits the endpoint.

- [ ] **Step 1: Write integration test fixtures and tests**

Append to `tests/test_isosurface.py`:

```python
from collections.abc import AsyncIterator
from unittest.mock import AsyncMock

import msgpack
import msgpack_numpy
import pytest_asyncio
from conftest import (
    MockSioServer,
    auth_header,
    create_test_room,
    create_test_user_in_db,
)
from httpx import ASGITransport, AsyncClient
from sqlalchemy.ext.asyncio import AsyncSession, async_sessionmaker, create_async_engine
from sqlalchemy.pool import StaticPool
from sqlmodel import SQLModel

from zndraw.storage import AsebytesStorage

msgpack_numpy.patch()


# =============================================================================
# Integration Test Fixtures
# =============================================================================


def _make_cube_data(
    n: int = 20,
) -> dict:
    """Create a sphere SDF cube data dict for testing."""
    lin = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(lin, lin, lin, indexing="ij")
    grid = np.sqrt(x**2 + y**2 + z**2) - 0.5
    return {
        "grid": grid,
        "origin": np.array([-1.0, -1.0, -1.0]),
        "cell": np.array([[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]),
    }


def _make_frame_with_cube(cube_key: str = "info.orbital_homo") -> dict[bytes, bytes]:
    """Create a raw frame containing volumetric cube data."""
    cube_data = _make_cube_data()
    return {
        cube_key.encode(): msgpack.packb(cube_data, default=msgpack_numpy.encode),
    }


@pytest_asyncio.fixture(name="iso_session")
async def iso_session_fixture() -> AsyncIterator[AsyncSession]:
    """Create a fresh in-memory async database session."""
    engine = create_async_engine(
        "sqlite+aiosqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
    )
    try:
        async with engine.begin() as conn:
            await conn.run_sync(SQLModel.metadata.create_all)
        factory = async_sessionmaker(
            bind=engine, class_=AsyncSession, expire_on_commit=False
        )
        async with factory() as session:
            yield session
    finally:
        await engine.dispose()


@pytest_asyncio.fixture(name="iso_storage")
async def iso_storage_fixture() -> AsyncIterator[AsebytesStorage]:
    """Create a fresh AsebytesStorage."""
    storage = AsebytesStorage("memory://")
    yield storage
    await storage.close()


@pytest_asyncio.fixture(name="iso_client")
async def iso_client_fixture(
    iso_session: AsyncSession, iso_storage: AsebytesStorage
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with storage + session overrides."""
    from contextlib import asynccontextmanager

    from zndraw_auth import get_session
    from zndraw_auth.settings import AuthSettings

    from zndraw.app import app
    from zndraw.dependencies import get_redis, get_storage, get_tsio

    mock_sio = MockSioServer()

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield iso_session

    @asynccontextmanager
    async def test_session_maker():
        yield iso_session

    app.state.auth_settings = AuthSettings()
    app.state.session_maker = test_session_maker
    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_storage] = lambda: iso_storage
    app.dependency_overrides[get_tsio] = lambda: mock_sio
    app.dependency_overrides[get_redis] = lambda: AsyncMock()

    async with AsyncClient(
        transport=ASGITransport(app=app), base_url="http://test"
    ) as client:
        yield client

    app.dependency_overrides.clear()


# =============================================================================
# Integration Tests
# =============================================================================


@pytest.mark.asyncio
async def test_isosurface_basic(
    iso_client: AsyncClient, iso_session: AsyncSession, iso_storage: AsebytesStorage
) -> None:
    """Store cube data, GET isosurface, verify 200 with vertices/faces."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    frame = _make_frame_with_cube("info.orbital_homo")
    await iso_storage.extend(room.id, [frame])  # type: ignore[arg-type]

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.orbital_homo", "isovalue": "0.0"},
        headers=auth_header(token),
    )
    assert response.status_code == 200
    assert response.headers["content-type"] == "application/x-msgpack"

    data = msgpack.unpackb(
        response.content, object_hook=msgpack_numpy.decode, raw=False
    )
    assert "vertices" in data
    assert "faces" in data
    assert len(data["vertices"]) > 0
    assert len(data["faces"]) > 0


@pytest.mark.asyncio
async def test_isosurface_missing_cube_key(
    iso_client: AsyncClient, iso_session: AsyncSession, iso_storage: AsebytesStorage
) -> None:
    """GET with nonexistent cube key returns 422."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    frame = _make_frame_with_cube("info.orbital_homo")
    await iso_storage.extend(room.id, [frame])  # type: ignore[arg-type]

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.nonexistent", "isovalue": "0.0"},
        headers=auth_header(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_isosurface_frame_not_found(
    iso_client: AsyncClient, iso_session: AsyncSession
) -> None:
    """GET with out-of-range frame index returns 404."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/99/isosurface",
        params={"cube_key": "info.foo", "isovalue": "0.0"},
        headers=auth_header(token),
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_isosurface_empty_surface(
    iso_client: AsyncClient, iso_session: AsyncSession, iso_storage: AsebytesStorage
) -> None:
    """Isovalue outside data range returns 200 with empty vertices/faces."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    frame = _make_frame_with_cube("info.orbital_homo")
    await iso_storage.extend(room.id, [frame])  # type: ignore[arg-type]

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.orbital_homo", "isovalue": "999.0"},
        headers=auth_header(token),
    )
    assert response.status_code == 200

    data = msgpack.unpackb(
        response.content, object_hook=msgpack_numpy.decode, raw=False
    )
    assert data["vertices"] == []
    assert data["faces"] == []


@pytest.mark.asyncio
async def test_isosurface_invalid_grid(
    iso_client: AsyncClient, iso_session: AsyncSession, iso_storage: AsebytesStorage
) -> None:
    """Cube data with non-3D grid returns 422."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    bad_cube = {
        "grid": np.ones((10, 10)),  # 2D, not 3D
        "origin": np.zeros(3),
        "cell": np.eye(3),
    }
    frame = {
        b"info.bad": msgpack.packb(bad_cube, default=msgpack_numpy.encode),
    }
    await iso_storage.extend(room.id, [frame])  # type: ignore[arg-type]

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.bad", "isovalue": "0.5"},
        headers=auth_header(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_isosurface_missing_dict_keys(
    iso_client: AsyncClient, iso_session: AsyncSession, iso_storage: AsebytesStorage
) -> None:
    """Cube data dict missing 'grid' key returns 422."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    bad_cube = {"origin": np.zeros(3), "cell": np.eye(3)}  # no 'grid'
    frame = {
        b"info.bad": msgpack.packb(bad_cube, default=msgpack_numpy.encode),
    }
    await iso_storage.extend(room.id, [frame])  # type: ignore[arg-type]

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.bad", "isovalue": "0.5"},
        headers=auth_header(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_isosurface_step_size(
    iso_client: AsyncClient, iso_session: AsyncSession, iso_storage: AsebytesStorage
) -> None:
    """Coarser step_size produces fewer vertices."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    # Use larger grid so step_size difference is visible
    cube_data = _make_cube_data(n=40)
    frame = {
        b"info.orb": msgpack.packb(cube_data, default=msgpack_numpy.encode),
    }
    await iso_storage.extend(room.id, [frame])  # type: ignore[arg-type]

    resp_fine = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.orb", "isovalue": "0.0", "step_size": "1"},
        headers=auth_header(token),
    )
    resp_coarse = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.orb", "isovalue": "0.0", "step_size": "2"},
        headers=auth_header(token),
    )

    assert resp_fine.status_code == 200
    assert resp_coarse.status_code == 200

    data_fine = msgpack.unpackb(
        resp_fine.content, object_hook=msgpack_numpy.decode, raw=False
    )
    data_coarse = msgpack.unpackb(
        resp_coarse.content, object_hook=msgpack_numpy.decode, raw=False
    )

    assert len(data_coarse["vertices"]) < len(data_fine["vertices"])


@pytest.mark.asyncio
async def test_isosurface_room_not_found(
    iso_client: AsyncClient, iso_session: AsyncSession
) -> None:
    """GET isosurface for non-existent room returns 404."""
    user, token = await create_test_user_in_db(iso_session)

    response = await iso_client.get(
        "/v1/rooms/nonexistent-room/frames/0/isosurface",
        params={"cube_key": "info.foo", "isovalue": "0.0"},
        headers=auth_header(token),
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_isosurface_pyscf_h2(
    iso_client: AsyncClient, iso_session: AsyncSession, iso_storage: AsebytesStorage
) -> None:
    """Use PySCF to generate a real H2 HOMO orbital, extract isosurface."""
    pyscf = pytest.importorskip("pyscf")
    from pyscf import gto, scf
    from pyscf.tools import cubegen

    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", unit="Angstrom")
    mf = scf.RHF(mol)
    mf.kernel()

    # Generate HOMO orbital on a grid
    homo_idx = mol.nelectron // 2 - 1
    cc = cubegen.Cube(mol, nx=30, ny=30, nz=30)
    orb_on_grid = cc.get_density(mf.mo_coeff[:, homo_idx])
    orb_on_grid = orb_on_grid.reshape(cc.nx, cc.ny, cc.nz)

    bohr_to_ang = 0.529177
    cube_data = {
        "grid": orb_on_grid,
        "origin": cc.boxorig * bohr_to_ang,
        "cell": cc.box * bohr_to_ang,
    }
    frame = {
        b"info.orbital_homo": msgpack.packb(
            cube_data, default=msgpack_numpy.encode
        ),
    }

    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)
    await iso_storage.extend(room.id, [frame])  # type: ignore[arg-type]

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.orbital_homo", "isovalue": "0.02"},
        headers=auth_header(token),
    )
    assert response.status_code == 200

    data = msgpack.unpackb(
        response.content, object_hook=msgpack_numpy.decode, raw=False
    )
    assert len(data["vertices"]) > 0
    assert len(data["faces"]) > 0

    # Vertices should be in a plausible coordinate range for H2 (~few Angstroms)
    verts = np.array(data["vertices"])
    assert np.all(np.abs(verts) < 10.0)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_isosurface.py::test_isosurface_basic -v`
Expected: FAIL (404 — route not registered yet).

- [ ] **Step 3: Implement the REST endpoint**

Replace `src/zndraw/routes/isosurface.py` with the full implementation (keeping `_extract_mesh` and adding the route):

```python
"""Isosurface REST endpoint — marching cubes extraction for volumetric data."""

from __future__ import annotations

from typing import Annotated

import msgpack
import msgpack_numpy
import numpy as np
from fastapi import APIRouter, Query, Response
from sqlalchemy.ext.asyncio import AsyncSession

from zndraw.dependencies import (
    CurrentUserFactoryDep,
    SessionMakerDep,
    StorageDep,
    verify_room,
)
from zndraw.exceptions import (
    FrameNotFound,
    RoomNotFound,
    UnprocessableContent,
    problem_responses,
)

router = APIRouter(
    prefix="/v1/rooms/{room_id}/frames/{index}/isosurface",
    tags=["isosurface"],
)


def _extract_mesh(
    grid: np.ndarray,
    origin: np.ndarray,
    cell: np.ndarray,
    isovalue: float,
    step_size: int,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Run marching cubes and transform vertices to world space.

    Parameters
    ----------
    grid : np.ndarray
        3D scalar field (Nx, Ny, Nz).
    origin : np.ndarray
        World-space origin (3,).
    cell : np.ndarray
        Axis vectors spanning the grid (3, 3).
    isovalue : float
        Scalar threshold for surface extraction.
    step_size : int
        Marching cubes step size (1 = full resolution).

    Returns
    -------
    tuple[np.ndarray, np.ndarray] | None
        (vertices, faces) as float32/uint32 arrays, or None if no surface.
    """
    from skimage.measure import marching_cubes

    try:
        verts, faces, _, _ = marching_cubes(grid, level=isovalue, step_size=step_size)
    except (ValueError, RuntimeError):
        return None

    if len(verts) == 0:
        return None

    shape = np.array(grid.shape, dtype=np.float64)
    fractional = verts / (shape - 1)
    world = origin + fractional @ cell

    return world.astype(np.float32), faces.astype(np.uint32)


@router.get(
    "",
    responses={
        **problem_responses(RoomNotFound, FrameNotFound),
        **problem_responses(UnprocessableContent),
    },
)
async def get_isosurface(
    session_maker: SessionMakerDep,
    storage: StorageDep,
    _current_user: CurrentUserFactoryDep,
    room_id: str,
    index: int,
    cube_key: Annotated[str, Query(description="Frame key for volumetric data dict")],
    isovalue: Annotated[float, Query(description="Scalar threshold")],
    step_size: Annotated[int, Query(ge=1, le=8, description="Marching cubes step")] = 1,
) -> Response:
    """Extract an isosurface mesh from volumetric frame data."""
    async with session_maker() as session:
        await verify_room(session, room_id)
        total = await storage.get_length(room_id)
        if index < 0 or index >= total:
            range_str = f"0-{total - 1}" if total > 0 else "no frames"
            raise FrameNotFound.exception(
                f"Frame index {index} out of range ({range_str})"
            )

    frame = await storage.get(room_id, index)
    if frame is None:
        raise FrameNotFound.exception(f"Frame at index {index} is empty")

    key_bytes = cube_key.encode()
    if key_bytes not in frame:
        raise UnprocessableContent.exception(
            f"Key '{cube_key}' not found in frame {index}"
        )

    cube_dict = msgpack.unpackb(
        frame[key_bytes], object_hook=msgpack_numpy.decode, raw=False
    )

    _REQUIRED_KEYS = {"grid", "origin", "cell"}
    if not isinstance(cube_dict, dict):
        raise UnprocessableContent.exception(
            f"Key '{cube_key}' is not a dict"
        )
    missing = _REQUIRED_KEYS - cube_dict.keys()
    if missing:
        raise UnprocessableContent.exception(
            f"Cube data missing keys: {', '.join(sorted(missing))}"
        )

    grid = np.asarray(cube_dict["grid"])
    if grid.ndim != 3:
        raise UnprocessableContent.exception(
            f"Grid must be 3D, got {grid.ndim}D"
        )

    origin = np.asarray(cube_dict["origin"], dtype=np.float64)
    cell = np.asarray(cube_dict["cell"], dtype=np.float64)

    result = _extract_mesh(grid, origin, cell, isovalue, step_size)

    if result is None:
        packed = msgpack.packb({"vertices": [], "faces": []})
    else:
        vertices, faces = result
        packed = msgpack.packb(
            {"vertices": vertices, "faces": faces},
            default=msgpack_numpy.encode,
        )

    return Response(content=packed, media_type="application/x-msgpack")
```

- [ ] **Step 4: Include the router in app.py**

Modify `src/zndraw/app.py`:

Add import (after the existing route imports, find the block that imports routers and add):

```python
from zndraw.routes.isosurface import router as isosurface_router
```

Add `app.include_router(isosurface_router)` in the router inclusion block (after `app.include_router(frames_router)`):

```python
app.include_router(isosurface_router)
```

- [ ] **Step 5: Run integration tests**

Run: `uv run pytest tests/test_isosurface.py -v`
Expected: All tests PASS (7 model + 4 unit + 9 integration = 20 total).

- [ ] **Step 6: Run linting**

Run: `uv run ruff check --select I --fix . && uv run ruff format .`
Expected: Clean.

- [ ] **Step 7: Run full test suite to verify no regressions**

Run: `uv run pytest tests/ -x -q`
Expected: All tests pass.

- [ ] **Step 8: Commit**

```bash
git add src/zndraw/routes/isosurface.py src/zndraw/app.py tests/test_isosurface.py
git commit -m "feat(isosurface): add REST endpoint with marching cubes extraction"
```

---

## Chunk 4: Frontend — Isosurface React Component + Canvas Registration

### Task 5: Create the Isosurface React component

The component fetches the mesh from the dedicated endpoint using `useQuery`, decodes the msgpack response, and renders a `BufferGeometry` with `MeshPhysicalMaterial`.

**Files:**
- Create: `frontend/src/components/three/Isosurface.tsx`

- [ ] **Step 1: Create the Isosurface component**

Create `frontend/src/components/three/Isosurface.tsx`:

```tsx
/**
 * Isosurface component for volumetric data visualization.
 *
 * Fetches mesh data from a dedicated server endpoint that runs marching cubes.
 * Does NOT use useRegisterFrameKeys or getFrameBatched — volumetric data stays
 * server-side. Only the extracted mesh (kilobytes) crosses the wire.
 */

import { useQuery } from "@tanstack/react-query";
import { useMemo } from "react";
import * as THREE from "three";
import { useAppStore } from "../../store";
import { getGeometryWithDefaults } from "../../utils/geometryDefaults";
import { unpackBinary } from "../../utils/msgpack-numpy";
import axios from "axios";
import { getToken } from "../../utils/auth";

interface IsosurfaceData {
	active: boolean;
	cube_key: string;
	isovalue: number;
	color: string;
	resolution: number;
	opacity: number;
}

async function fetchIsosurface(
	roomId: string,
	frame: number,
	cubeKey: string,
	isovalue: number,
	stepSize: number,
): Promise<{ vertices: number[][]; faces: number[][] }> {
	const params = new URLSearchParams({
		cube_key: cubeKey,
		isovalue: isovalue.toString(),
		step_size: stepSize.toString(),
	});
	const token = getToken();
	const response = await axios.get(
		`/v1/rooms/${roomId}/frames/${frame}/isosurface?${params.toString()}`,
		{
			responseType: "arraybuffer",
			headers: token ? { Authorization: `Bearer ${token}` } : {},
		},
	);
	return unpackBinary(new Uint8Array(response.data)) as {
		vertices: number[][] | Float32Array;
		faces: number[][] | Uint32Array;
	};
}

export default function Isosurface({
	data,
}: {
	data: IsosurfaceData;
	geometryKey: string;
}) {
	const geometryDefaults = useAppStore((state) => state.geometryDefaults);
	const roomId = useAppStore((state) => state.roomId);
	const currentFrame = useAppStore((state) => state.currentFrame);

	const fullData = getGeometryWithDefaults<IsosurfaceData>(
		data,
		"Isosurface",
		geometryDefaults,
	);

	const { data: meshData } = useQuery({
		queryKey: [
			"isosurface",
			roomId,
			currentFrame,
			fullData.cube_key,
			fullData.isovalue,
			fullData.resolution,
		],
		queryFn: () =>
			fetchIsosurface(
				roomId!,
				currentFrame,
				fullData.cube_key,
				fullData.isovalue,
				fullData.resolution,
			),
		enabled: !!roomId && !!fullData.cube_key && fullData.active,
		staleTime: 30000,
	});

	const geometry = useMemo(() => {
		if (!meshData || meshData.vertices.length === 0) return null;

		const geo = new THREE.BufferGeometry();

		// msgpack-numpy may return Float32Array or nested number[][]
		const vertices =
			meshData.vertices instanceof Float32Array
				? meshData.vertices
				: new Float32Array(
						Array.isArray(meshData.vertices[0])
							? (meshData.vertices as number[][]).flat()
							: (meshData.vertices as unknown as number[]),
					);
		geo.setAttribute("position", new THREE.BufferAttribute(vertices, 3));

		const indices =
			meshData.faces instanceof Uint32Array
				? meshData.faces
				: new Uint32Array(
						Array.isArray(meshData.faces[0])
							? (meshData.faces as number[][]).flat()
							: (meshData.faces as unknown as number[]),
					);
		geo.setIndex(new THREE.BufferAttribute(indices, 1));

		geo.computeVertexNormals();
		return geo;
	}, [meshData]);

	if (!fullData.active || !geometry) return null;

	const isTransparent = fullData.opacity < 1.0;

	return (
		<mesh geometry={geometry}>
			<meshPhysicalMaterial
				color={fullData.color}
				side={THREE.DoubleSide}
				transparent={isTransparent}
				opacity={fullData.opacity}
				depthWrite={!isTransparent}
				roughness={0.4}
			/>
		</mesh>
	);
}
```

- [ ] **Step 2: Verify the file compiles**

Run: `cd /Users/fzills/tools/zndraw-fastapi/frontend && bun run tsc --noEmit --pretty 2>&1 | head -20`
Expected: No errors in Isosurface.tsx (other pre-existing errors may appear).

- [ ] **Step 3: Commit**

```bash
git add frontend/src/components/three/Isosurface.tsx
git commit -m "feat(isosurface): add React/Three.js rendering component"
```

---

### Task 6: Register Isosurface in Canvas.tsx

**Files:**
- Modify: `frontend/src/components/Canvas.tsx`

- [ ] **Step 1: Add import**

In `frontend/src/components/Canvas.tsx`, add the import after the existing component imports (after the `import { Fog } from "./three/Fog";` line):

```typescript
import Isosurface from "./three/Isosurface";
```

- [ ] **Step 2: Add to SIMPLE_GEOMETRY_COMPONENTS**

In the `SIMPLE_GEOMETRY_COMPONENTS` object, add `Isosurface` (after `Fog: Fog,`):

```typescript
	Isosurface: Isosurface,
```

The full object should now be:

```typescript
const SIMPLE_GEOMETRY_COMPONENTS = {
	Camera: Camera,
	DirectionalLight: DirectionalLight,
	AmbientLight: AmbientLight,
	HemisphereLight: HemisphereLight,
	Fog: Fog,
	Isosurface: Isosurface,
} as const;
```

- [ ] **Step 3: Verify compilation**

Run: `cd /Users/fzills/tools/zndraw-fastapi/frontend && bun run tsc --noEmit --pretty 2>&1 | head -20`
Expected: No new errors.

- [ ] **Step 4: Run linting and format**

Run: `cd /Users/fzills/tools/zndraw-fastapi/frontend && bun run lint 2>&1 | tail -5`
Expected: Clean (or pre-existing issues only).

- [ ] **Step 5: Commit**

```bash
git add frontend/src/components/Canvas.tsx
git commit -m "feat(isosurface): register component in Canvas geometry dispatch"
```

---

## Chunk 5: Final Verification

### Task 7: Full test suite and type checking

- [ ] **Step 1: Run full Python test suite**

Run: `uv run pytest tests/ -x -q`
Expected: All tests pass.

- [ ] **Step 2: Run pyright**

Run: `uv run pyright src/zndraw/geometries/isosurface.py src/zndraw/routes/isosurface.py`
Expected: No errors (or only pre-existing suppressions).

- [ ] **Step 3: Run ruff**

Run: `uv run ruff check . && uv run ruff format --check .`
Expected: Clean.

- [ ] **Step 4: Run frontend type check**

Run: `cd /Users/fzills/tools/zndraw-fastapi/frontend && bun run tsc --noEmit`
Expected: No new errors in Isosurface.tsx or Canvas.tsx.

- [ ] **Step 5: Final commit (if any fixes needed)**

If any fixes were required, commit them:
```bash
git add -A
git commit -m "fix(isosurface): address type/lint issues"
```
