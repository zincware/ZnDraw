"""Tests for ASE constraint serialization through encode/decode, storage, and client."""

import uuid

import ase
import numpy as np
import pytest
import pytest_asyncio
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms, FixBondLengths, FixedLine, FixedPlane
from asebytes import decode, encode

from zndraw.storage import InMemoryStorage, LMDBStorage

# =============================================================================
# Section 1: Unit tests (encode/decode roundtrip, no server needed)
# =============================================================================


def test_fixatoms_roundtrip():
    """FixAtoms indices survive encode/decode roundtrip."""
    atoms = ase.Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms.set_constraint(FixAtoms(indices=[0, 2]))

    decoded = decode(encode(atoms))

    assert len(decoded.constraints) == 1
    constraint = decoded.constraints[0]
    assert isinstance(constraint, FixAtoms)
    np.testing.assert_array_equal(constraint.index, [0, 2])


def test_fixbondlengths_roundtrip():
    """FixBondLengths pairs survive encode/decode roundtrip."""
    atoms = ase.Atoms("H4", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
    atoms.set_constraint(FixBondLengths(pairs=[[0, 1], [2, 3]]))

    decoded = decode(encode(atoms))

    assert len(decoded.constraints) == 1
    constraint = decoded.constraints[0]
    assert isinstance(constraint, FixBondLengths)
    np.testing.assert_array_equal(constraint.pairs, [[0, 1], [2, 3]])


def test_fixedline_roundtrip():
    """FixedLine direction and index survive encode/decode roundtrip."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms.set_constraint(FixedLine(0, direction=[1, 0, 0]))

    decoded = decode(encode(atoms))

    assert len(decoded.constraints) == 1
    constraint = decoded.constraints[0]
    assert isinstance(constraint, FixedLine)
    np.testing.assert_array_equal(constraint.index, [0])
    np.testing.assert_allclose(constraint.dir, [1.0, 0.0, 0.0])


def test_fixedplane_roundtrip():
    """FixedPlane normal and index survive encode/decode roundtrip."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms.set_constraint(FixedPlane(0, direction=[0, 0, 1]))

    decoded = decode(encode(atoms))

    assert len(decoded.constraints) == 1
    constraint = decoded.constraints[0]
    assert isinstance(constraint, FixedPlane)
    np.testing.assert_array_equal(constraint.index, [0])
    np.testing.assert_allclose(constraint.dir, [0.0, 0.0, 1.0])


def test_multiple_constraints():
    """Multiple different constraints on the same Atoms object are preserved."""
    atoms = ase.Atoms("H4", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
    atoms.set_constraint(
        [
            FixAtoms(indices=[0, 3]),
            FixedLine(1, direction=[0, 1, 0]),
        ]
    )

    decoded = decode(encode(atoms))

    assert len(decoded.constraints) == 2
    types = {type(c) for c in decoded.constraints}
    assert types == {FixAtoms, FixedLine}

    # Verify each constraint's data
    for c in decoded.constraints:
        if isinstance(c, FixAtoms):
            np.testing.assert_array_equal(c.index, [0, 3])
        elif isinstance(c, FixedLine):
            np.testing.assert_array_equal(c.index, [1])
            np.testing.assert_allclose(c.dir, [0.0, 1.0, 0.0])


def test_no_constraints():
    """Atoms without constraints roundtrip with an empty constraints list."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])

    decoded = decode(encode(atoms))

    assert decoded.constraints == []


def test_constraints_with_calculator():
    """Constraints and SinglePointCalculator both survive encode/decode."""
    atoms = ase.Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms.set_constraint(FixAtoms(indices=[1]))
    energy = -5.0
    calc = SinglePointCalculator(
        atoms,
        energy=energy,
        forces=np.array([[0.1, 0, 0], [0, 0.2, 0], [0, 0, 0.3]]),
    )
    atoms.calc = calc

    decoded = decode(encode(atoms))

    # Constraints preserved
    assert len(decoded.constraints) == 1
    assert isinstance(decoded.constraints[0], FixAtoms)
    np.testing.assert_array_equal(decoded.constraints[0].index, [1])

    # Calculator preserved
    assert decoded.calc is not None
    assert decoded.get_potential_energy() == energy


# =============================================================================
# Section 2: Storage integration tests (async)
# =============================================================================


@pytest_asyncio.fixture
async def memory_storage():
    """Create a fresh InMemoryStorage instance."""
    storage = InMemoryStorage()
    yield storage
    await storage.close()


@pytest_asyncio.fixture
async def lmdb_storage(tmp_path):
    """Create a fresh LMDBStorage instance with a temp directory."""
    storage = LMDBStorage(path=tmp_path / "test.lmdb")
    yield storage
    await storage.close()


@pytest.mark.asyncio
async def test_constraints_through_memory_storage(memory_storage: InMemoryStorage):
    """Constraints survive InMemoryStorage extend/get roundtrip."""
    atoms = ase.Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms.set_constraint(FixAtoms(indices=[0, 2]))

    frame = encode(atoms)
    await memory_storage.extend(room_id="room1", frames=[frame])

    raw = await memory_storage.get(room_id="room1", index=0)
    assert raw is not None
    decoded = decode(raw)

    assert len(decoded.constraints) == 1
    assert isinstance(decoded.constraints[0], FixAtoms)
    np.testing.assert_array_equal(decoded.constraints[0].index, [0, 2])


@pytest.mark.asyncio
async def test_constraints_through_lmdb_storage(lmdb_storage: LMDBStorage):
    """Constraints survive LMDBStorage extend/get roundtrip."""
    atoms = ase.Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms.set_constraint(FixAtoms(indices=[0, 2]))

    frame = encode(atoms)
    await lmdb_storage.extend(room_id="room1", frames=[frame])

    raw = await lmdb_storage.get(room_id="room1", index=0)
    assert raw is not None
    decoded = decode(raw)

    assert len(decoded.constraints) == 1
    assert isinstance(decoded.constraints[0], FixAtoms)
    np.testing.assert_array_equal(decoded.constraints[0].index, [0, 2])


async def _assert_variable_constraints_across_frames(
    storage: InMemoryStorage | LMDBStorage,
) -> None:
    """Verify different constraint types across frames are each preserved."""
    # Frame 0: FixAtoms
    atoms0 = ase.Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms0.set_constraint(FixAtoms(indices=[1]))

    # Frame 1: FixedLine
    atoms1 = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms1.set_constraint(FixedLine(0, direction=[1, 0, 0]))

    # Frame 2: FixedPlane
    atoms2 = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    atoms2.set_constraint(FixedPlane(1, direction=[0, 0, 1]))

    frames = [encode(a) for a in [atoms0, atoms1, atoms2]]
    await storage.extend(room_id="room1", frames=frames)

    # Verify frame 0
    raw0 = await storage.get(room_id="room1", index=0)
    assert raw0 is not None
    dec0 = decode(raw0)
    assert len(dec0.constraints) == 1
    assert isinstance(dec0.constraints[0], FixAtoms)
    np.testing.assert_array_equal(dec0.constraints[0].index, [1])

    # Verify frame 1
    raw1 = await storage.get(room_id="room1", index=1)
    assert raw1 is not None
    dec1 = decode(raw1)
    assert len(dec1.constraints) == 1
    assert isinstance(dec1.constraints[0], FixedLine)
    np.testing.assert_allclose(dec1.constraints[0].dir, [1.0, 0.0, 0.0])

    # Verify frame 2
    raw2 = await storage.get(room_id="room1", index=2)
    assert raw2 is not None
    dec2 = decode(raw2)
    assert len(dec2.constraints) == 1
    assert isinstance(dec2.constraints[0], FixedPlane)
    np.testing.assert_array_equal(dec2.constraints[0].index, [1])
    np.testing.assert_allclose(dec2.constraints[0].dir, [0.0, 0.0, 1.0])


@pytest.mark.asyncio
async def test_variable_constraints_across_frames_memory(
    memory_storage: InMemoryStorage,
):
    """Different constraint types across InMemoryStorage frames are preserved."""
    await _assert_variable_constraints_across_frames(memory_storage)


@pytest.mark.asyncio
async def test_variable_constraints_across_frames_lmdb(lmdb_storage: LMDBStorage):
    """Different constraint types across LMDBStorage frames are preserved."""
    await _assert_variable_constraints_across_frames(lmdb_storage)


# =============================================================================
# Section 3: Client integration tests (sync, real server)
# =============================================================================


def test_client_constraint_roundtrip(server: str):
    """Constraints survive the full client append/retrieve cycle."""
    from zndraw import ZnDraw

    room_id = uuid.uuid4().hex
    client = ZnDraw(url=server, room=room_id)

    atoms = ase.Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    atoms.set_constraint(FixAtoms(indices=[0, 2]))
    client.append(atoms)

    retrieved = client[0]
    assert len(retrieved.constraints) == 1
    constraint = retrieved.constraints[0]
    assert isinstance(constraint, FixAtoms)
    np.testing.assert_array_equal(constraint.index, [0, 2])

    client.disconnect()


def test_client_extend_with_constraints(server: str):
    """Multiple constrained frames survive client extend/retrieve."""
    from zndraw import ZnDraw

    room_id = uuid.uuid4().hex
    client = ZnDraw(url=server, room=room_id)

    # Frame 0: FixAtoms
    a0 = ase.Atoms("H3", positions=[[0, 0, 0], [1, 0, 0], [2, 0, 0]])
    a0.set_constraint(FixAtoms(indices=[1]))

    # Frame 1: FixedLine
    a1 = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])
    a1.set_constraint(FixedLine(0, direction=[0, 1, 0]))

    # Frame 2: no constraints
    a2 = ase.Atoms("H2", positions=[[0, 0, 0], [1, 0, 0]])

    client.extend([a0, a1, a2])
    assert len(client) == 3

    # Verify frame 0
    r0 = client[0]
    assert len(r0.constraints) == 1
    assert isinstance(r0.constraints[0], FixAtoms)
    np.testing.assert_array_equal(r0.constraints[0].index, [1])

    # Verify frame 1
    r1 = client[1]
    assert len(r1.constraints) == 1
    assert isinstance(r1.constraints[0], FixedLine)
    np.testing.assert_allclose(r1.constraints[0].dir, [0.0, 1.0, 0.0])

    # Verify frame 2
    r2 = client[2]
    assert r2.constraints == []

    client.disconnect()


# =============================================================================
# Section 4: Default geometry tests
# =============================================================================


def test_default_constraint_geometry_created_on_room_creation(server: str):
    """New rooms include a constraints-fixed-atoms geometry."""
    from zndraw import ZnDraw
    from zndraw.geometries import Sphere
    from zndraw.transformations import InArrayTransform

    room_id = uuid.uuid4().hex
    client = ZnDraw(url=server, room=room_id)

    geo = client.geometries["constraints-fixed-atoms"]
    assert geo is not None

    # Verify it's a Sphere with InArrayTransform position
    assert isinstance(geo, Sphere)
    assert isinstance(geo.position, InArrayTransform)
    assert geo.position.source == "constraints"
    assert geo.position.path == "0.kwargs.indices"
    assert geo.position.filter == "arrays.positions"

    assert isinstance(geo.radius, InArrayTransform)
    assert geo.radius.source == "constraints"
    assert geo.radius.filter == "arrays.radii"

    assert geo.color == ["#FF0000"]
    assert geo.selecting.enabled is False
    assert geo.hovering.enabled is False

    client.disconnect()


def test_constraint_geometry_renders_fixed_atoms(server: str):
    """Constraint geometry filters positions to only fixed atoms."""
    from zndraw import ZnDraw
    from zndraw.transformations import InArrayTransform

    room_id = uuid.uuid4().hex
    client = ZnDraw(url=server, room=room_id)

    atoms = ase.Atoms("H5", positions=[[i, 0, 0] for i in range(5)])
    atoms.set_constraint(FixAtoms(indices=[1, 3]))
    client.append(atoms)

    # The constraint geometry exists and has correct transform config
    geo = client.geometries["constraints-fixed-atoms"]
    assert isinstance(geo.position, InArrayTransform)

    # Verify constraint data roundtrip
    retrieved = client[0]
    assert len(retrieved.constraints) == 1
    assert isinstance(retrieved.constraints[0], FixAtoms)
    np.testing.assert_array_equal(retrieved.constraints[0].index, [1, 3])

    client.disconnect()
