"""Tests for Isosurface geometry model and endpoint."""

from collections.abc import AsyncIterator
from unittest.mock import AsyncMock

import msgpack
import msgpack_numpy
import numpy as np
import pytest
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

from zndraw.storage import FrameStorage

# =============================================================================
# Model Validation Tests
# =============================================================================


def test_isosurface_defaults():
    """Default construction produces expected field values."""
    from zndraw.geometries.isosurface import Isosurface

    iso = Isosurface()
    assert iso.cube_key == ""
    assert iso.isovalue == 0.02
    assert iso.isovalue_min == -0.25
    assert iso.isovalue_max == 0.25
    assert iso.sigma == 0.0
    assert iso.color == "#2244CC"
    assert iso.resolution == 1.0
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
        resolution=0.5,
        opacity=0.8,
        owner="user123",
    )
    assert iso.cube_key == "info.orbital_homo"
    assert iso.isovalue == -0.05
    assert iso.color == "#FF0000"
    assert iso.resolution == 0.5
    assert iso.opacity == 0.8
    assert iso.owner == "user123"


def test_isosurface_custom_range():
    """Construction with custom isovalue range for SDF data."""
    from zndraw.geometries.isosurface import Isosurface

    iso = Isosurface(
        isovalue=5.0,
        isovalue_min=0.0,
        isovalue_max=30.0,
    )
    assert iso.isovalue == 5.0
    assert iso.isovalue_min == 0.0
    assert iso.isovalue_max == 30.0


@pytest.mark.parametrize(
    ("isovalue", "isovalue_min", "isovalue_max"),
    [
        (-1.0, -0.25, 0.25),  # isovalue below min
        (1.0, -0.25, 0.25),  # isovalue above max
        (0.0, 0.5, -0.5),  # min > max
    ],
)
def test_isosurface_range_validation(
    isovalue: float, isovalue_min: float, isovalue_max: float
):
    """Isovalue must be within [isovalue_min, isovalue_max], and min <= max."""
    from pydantic import ValidationError

    from zndraw.geometries.isosurface import Isosurface

    with pytest.raises(ValidationError):
        Isosurface(
            isovalue=isovalue,
            isovalue_min=isovalue_min,
            isovalue_max=isovalue_max,
        )


def test_isosurface_sigma_positive():
    """Sigma must be non-negative."""
    from pydantic import ValidationError

    from zndraw.geometries.isosurface import Isosurface

    Isosurface(sigma=1.5)  # valid
    with pytest.raises(ValidationError):
        Isosurface(sigma=-0.1)


def test_isosurface_frozen():
    """Model is immutable (frozen=True)."""
    from pydantic import ValidationError

    from zndraw.geometries.isosurface import Isosurface

    iso = Isosurface()
    with pytest.raises(ValidationError):
        iso.isovalue = 0.5  # type: ignore[misc]


def test_isosurface_resolution_bounds():
    """Resolution must be between 0.0 and 1.0."""
    from pydantic import ValidationError

    from zndraw.geometries.isosurface import Isosurface

    with pytest.raises(ValidationError):
        Isosurface(resolution=-0.1)
    with pytest.raises(ValidationError):
        Isosurface(resolution=1.1)


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


def test_isosurface_schema_editable_range():
    """isovalue field schema uses editable-range with min/max field refs."""
    from zndraw.geometries.isosurface import Isosurface

    schema = Isosurface.model_json_schema()
    iso_props = schema["properties"]["isovalue"]
    assert iso_props["x-custom-type"] == "editable-range"
    assert iso_props["x-min-field"] == "isovalue_min"
    assert iso_props["x-max-field"] == "isovalue_max"
    assert iso_props["step"] == 0.001
    # No ge/le constraints — validation is via model validator
    assert "minimum" not in iso_props
    assert "maximum" not in iso_props


def test_isosurface_schema_hidden_fields():
    """isovalue_min and isovalue_max are hidden from the form."""
    from zndraw.geometries.isosurface import Isosurface

    schema = Isosurface.model_json_schema()
    for field_name in ("isovalue_min", "isovalue_max"):
        assert schema["properties"][field_name]["x-hidden"] is True


def test_isosurface_in_geometries_dict():
    """Isosurface is registered in the geometries dict."""
    from zndraw.geometries import geometries

    assert "Isosurface" in geometries


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


@pytest_asyncio.fixture(name="iso_client")
async def iso_client_fixture(
    iso_session: AsyncSession, frame_storage: FrameStorage
) -> AsyncIterator[AsyncClient]:
    """Create an async test client with storage + session overrides."""
    from contextlib import asynccontextmanager

    from zndraw_auth import get_session
    from zndraw_auth.settings import AuthSettings
    from zndraw_joblib.settings import JobLibSettings

    from zndraw.app import app
    from zndraw.dependencies import (
        get_frame_storage,
        get_joblib_settings,
        get_redis,
        get_result_backend,
        get_tsio,
    )

    mock_sio = MockSioServer()

    async def get_session_override() -> AsyncIterator[AsyncSession]:
        yield iso_session

    @asynccontextmanager
    async def test_session_maker():
        yield iso_session

    app.state.auth_settings = AuthSettings()
    app.state.session_maker = test_session_maker
    app.dependency_overrides[get_session] = get_session_override
    app.dependency_overrides[get_frame_storage] = lambda: frame_storage
    app.dependency_overrides[get_tsio] = lambda: mock_sio
    app.dependency_overrides[get_redis] = lambda: AsyncMock()
    app.dependency_overrides[get_result_backend] = lambda: AsyncMock()
    app.dependency_overrides[get_joblib_settings] = lambda: JobLibSettings()

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
    iso_client: AsyncClient, iso_session: AsyncSession, frame_storage: FrameStorage
) -> None:
    """Store cube data, GET isosurface, verify 200 with vertices/faces."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    frame = _make_frame_with_cube("info.orbital_homo")
    await frame_storage[room.id].extend([frame])

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
    iso_client: AsyncClient, iso_session: AsyncSession, frame_storage: FrameStorage
) -> None:
    """GET with nonexistent cube key returns 422."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    frame = _make_frame_with_cube("info.orbital_homo")
    await frame_storage[room.id].extend([frame])

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
    iso_client: AsyncClient, iso_session: AsyncSession, frame_storage: FrameStorage
) -> None:
    """Isovalue outside data range returns 200 with empty vertices/faces."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    frame = _make_frame_with_cube("info.orbital_homo")
    await frame_storage[room.id].extend([frame])

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
    iso_client: AsyncClient, iso_session: AsyncSession, frame_storage: FrameStorage
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
    await frame_storage[room.id].extend([frame])

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.bad", "isovalue": "0.5"},
        headers=auth_header(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_isosurface_missing_dict_keys(
    iso_client: AsyncClient, iso_session: AsyncSession, frame_storage: FrameStorage
) -> None:
    """Cube data dict missing 'grid' key returns 422."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    bad_cube = {"origin": np.zeros(3), "cell": np.eye(3)}  # no 'grid'
    frame = {
        b"info.bad": msgpack.packb(bad_cube, default=msgpack_numpy.encode),
    }
    await frame_storage[room.id].extend([frame])

    response = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.bad", "isovalue": "0.5"},
        headers=auth_header(token),
    )
    assert response.status_code == 422


@pytest.mark.asyncio
async def test_isosurface_resolution(
    iso_client: AsyncClient, iso_session: AsyncSession, frame_storage: FrameStorage
) -> None:
    """Higher resolution produces more vertices."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    # Use larger grid so resolution difference is visible
    cube_data = _make_cube_data(n=40)
    frame = {
        b"info.orb": msgpack.packb(cube_data, default=msgpack_numpy.encode),
    }
    await frame_storage[room.id].extend([frame])

    resp_fine = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.orb", "isovalue": "0.0", "resolution": "1.0"},
        headers=auth_header(token),
    )
    resp_coarse = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.orb", "isovalue": "0.0", "resolution": "0.0"},
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
async def test_isosurface_sigma_smoothing(
    iso_client: AsyncClient, iso_session: AsyncSession, frame_storage: FrameStorage
) -> None:
    """Gaussian smoothing with sigma > 0 produces a valid mesh from noisy data."""
    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)

    # Create a noisy sphere grid that fragments without smoothing
    rng = np.random.default_rng(42)
    n = 30
    lin = np.linspace(-1, 1, n)
    x, y, z = np.meshgrid(lin, lin, lin, indexing="ij")
    grid = np.sqrt(x**2 + y**2 + z**2) - 0.5 + rng.normal(0, 0.3, (n, n, n))
    cube_data = {
        "grid": grid,
        "origin": np.array([-1.0, -1.0, -1.0]),
        "cell": np.array([[2.0, 0, 0], [0, 2.0, 0], [0, 0, 2.0]]),
    }
    frame = {
        b"info.noisy": msgpack.packb(cube_data, default=msgpack_numpy.encode),
    }
    await frame_storage[room.id].extend([frame])

    # Without smoothing
    resp_raw = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.noisy", "isovalue": "0.0"},
        headers=auth_header(token),
    )
    assert resp_raw.status_code == 200

    # With smoothing
    resp_smooth = await iso_client.get(
        f"/v1/rooms/{room.id}/frames/0/isosurface",
        params={"cube_key": "info.noisy", "isovalue": "0.0", "sigma": "1.0"},
        headers=auth_header(token),
    )
    assert resp_smooth.status_code == 200

    data_raw = msgpack.unpackb(
        resp_raw.content, object_hook=msgpack_numpy.decode, raw=False
    )
    data_smooth = msgpack.unpackb(
        resp_smooth.content, object_hook=msgpack_numpy.decode, raw=False
    )
    # Smoothed mesh should have fewer vertices (less fragmented)
    assert len(data_smooth["vertices"]) > 0
    assert len(data_smooth["vertices"]) < len(data_raw["vertices"])


@pytest.mark.asyncio
async def test_isosurface_room_not_found(
    iso_client: AsyncClient, iso_session: AsyncSession
) -> None:
    """GET isosurface for non-existent room returns 404."""
    _, token = await create_test_user_in_db(iso_session)

    response = await iso_client.get(
        "/v1/rooms/nonexistent-room/frames/0/isosurface",
        params={"cube_key": "info.foo", "isovalue": "0.0"},
        headers=auth_header(token),
    )
    assert response.status_code == 404


@pytest.mark.asyncio
async def test_isosurface_pyscf_h2(
    iso_client: AsyncClient, iso_session: AsyncSession, frame_storage: FrameStorage
) -> None:
    """Use PySCF to generate a real H2 HOMO orbital, extract isosurface."""
    pytest.importorskip("pyscf")
    from pyscf import gto, lib, scf
    from pyscf.tools.cubegen import Cube

    mol = gto.M(atom="H 0 0 0; H 0 0 0.74", basis="sto-3g", unit="Angstrom")
    mf = scf.RHF(mol)
    mf.kernel()

    # Generate HOMO orbital on a grid
    homo_idx = mol.nelectron // 2 - 1
    cc = Cube(mol, nx=30, ny=30, nz=30)
    coeff = mf.mo_coeff[:, homo_idx]
    coords = cc.get_coords()
    ngrids = cc.get_ngrids()
    blksize = min(8000, ngrids)
    orb_on_grid = np.empty(ngrids)
    for ip0, ip1 in lib.prange(0, ngrids, blksize):
        ao = mol.eval_gto("GTOval", coords[ip0:ip1])
        orb_on_grid[ip0:ip1] = ao @ coeff
    orb_on_grid = orb_on_grid.reshape(cc.nx, cc.ny, cc.nz)

    bohr_to_ang = 0.529177
    cube_data = {
        "grid": orb_on_grid,
        "origin": cc.boxorig * bohr_to_ang,
        "cell": cc.box * bohr_to_ang,
    }
    frame = {
        b"info.orbital_homo": msgpack.packb(cube_data, default=msgpack_numpy.encode),
    }

    user, token = await create_test_user_in_db(iso_session)
    room = await create_test_room(iso_session, user)
    await frame_storage[room.id].extend([frame])

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
