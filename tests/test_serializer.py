import ase
import matplotlib.pyplot as plt
import numpy as np
import numpy.testing as npt
import znjson
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms
import rdkit2ase
import plotly.graph_objects as go

from zndraw import Figure
from zndraw.converter import ASEConverter
import pytest


@pytest.fixture
def atoms():
    """Create a simple ASE Atoms object for testing."""
    return ase.Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])


@pytest.fixture
def atoms_x_info_str(atoms):
    """Create an ASE Atoms object with additional info."""
    atoms.info["key"] = "value"
    return atoms


@pytest.fixture
def atoms_x_info_nested_dict(atoms):
    """Create an ASE Atoms object with additional info."""
    atoms.info["key"] = {"subkey": "subvalue"}
    atoms.info["nested"] = {"key1": "value1", "key2": {"subkey": "subvalue"}, "key3": [1, 2, 3]}
    return atoms


@pytest.fixture
def atoms_x_info_nested_list(atoms):
    """Create an ASE Atoms object with a list in info."""
    atoms.info["key"] = ["value1", "value2"]
    atoms.info["nested"] = [{"subkey": "subvalue"}, [1, 2, 3]]
    return atoms


@pytest.fixture
def atoms_x_info_nested_objects(atoms):
    """Create an ASE Atoms object with nested objects in info."""
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 2, 3])

    atoms.info["nested"] = {
        "lvl1": {
            "lvl2": {
                "lvl3": {
                    "array": np.array([1, 2, 3]),
                    "figure": Figure.from_matplotlib(fig),
                }
            }
        }
    }
    return atoms


@pytest.fixture
def atoms_x_info_figure(atoms):
    """Create an ASE Atoms object with a Figure in info."""
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 2, 3])
    atoms.info["figure"] = Figure.from_matplotlib(fig)
    return atoms

@pytest.fixture
def atoms_x_info_plotly(atoms):
    """Create an ASE Atoms object with a Plotly figure in info."""
    fig = go.Figure(data=go.Scatter(x=[1, 2, 3], y=[1, 2, 3]))
    atoms.info["plotly_figure"] = fig
    return atoms


@pytest.fixture
def atoms_x_info_array(atoms):
    """Create an ASE Atoms object with an array in info."""
    atoms.info["array"] = np.array([1, 2, 3])
    return atoms

@pytest.fixture
def atoms_x_info_np_generic(atoms):
    """Create an ASE Atoms object with a numpy generic type in info."""
    atoms.info["generic"] = np.int64(42)
    return atoms

@pytest.fixture
def atoms_x_info_znsocket(atoms):
    """Create an ASE Atoms object with a znsocket in info."""
    

@pytest.fixture
def atoms_x_connectivity():
    """Create an ASE Atoms object with connectivity info."""
    atoms = rdkit2ase.smiles2atoms("CCO")
    assert "connectivity" in atoms.info
    return atoms


@pytest.fixture
def atoms_x_calc(atoms):
    """Create an ASE Atoms object with a calculator."""
    atoms.calc = SinglePointCalculator(atoms, energy=0.0, forces=np.zeros((2, 3)))
    return atoms


@pytest.fixture
def atoms_x_constraints(atoms):
    """Create an ASE Atoms object with constraints."""
    atoms.set_constraint(FixAtoms([0]))
    return atoms


@pytest.fixture
def atoms_x_arrays(atoms):
    """Create an ASE Atoms object with additional arrays."""
    atoms.arrays["colors"] = np.array(["#ff0000", "#00ff00"])
    atoms.arrays["radii"] = np.array([0.3, 0.4])
    return atoms


@pytest.fixture
def atoms_x_cell(atoms):
    """Create an ASE Atoms object with a cell."""
    atoms.set_cell([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    atoms.set_pbc([True, True, False])
    return atoms


@pytest.mark.parametrize(
    "myatoms",
    [
        "atoms_x_info_str",
        "atoms_x_info_nested_dict",
        "atoms_x_info_nested_list",
        "atoms_x_info_nested_objects",
        "atoms_x_info_figure",
        "atoms_x_info_array",
        "atoms_x_info_np_generic",
        "atoms_x_info_plotly",
        "atoms_x_connectivity",
        "atoms_x_calc",
        "atoms_x_constraints",
        "atoms_x_arrays",
        "atoms_x_cell",
    ],
)
def test_serialization(myatoms, request):
    """Test serialization and deserialization of ASE Atoms objects."""
    atoms = request.getfixturevalue(myatoms)
    serialized = znjson.dumps(atoms, cls=znjson.ZnEncoder.from_converters([ASEConverter]))
    deserialized = znjson.loads(serialized, cls=znjson.ZnDecoder.from_converters([ASEConverter]))

    npt.assert_array_equal(atoms.get_atomic_numbers(), deserialized.get_atomic_numbers())
    npt.assert_array_equal(atoms.get_positions(), deserialized.get_positions())
    npt.assert_array_equal(atoms.get_pbc(), deserialized.get_pbc())
    npt.assert_array_equal(atoms.get_cell(), deserialized.get_cell())

    for key in atoms.info:
        npt.assert_equal(atoms.info[key], deserialized.info[key])
    for key in deserialized.info:
        npt.assert_equal(atoms.info[key], deserialized.info[key])
    for key in atoms.arrays:
        npt.assert_equal(atoms.arrays[key], deserialized.arrays[key])
    for key in deserialized.arrays:
        npt.assert_equal(atoms.arrays[key], deserialized.arrays[key])
    if atoms.calc is not None:
        for key in atoms.calc.results:
            npt.assert_equal(
                atoms.calc.results[key], deserialized.calc.results[key]
            )
    if deserialized.calc is not None:
        for key in deserialized.calc.results:
            npt.assert_equal(
                atoms.calc.results[key], deserialized.calc.results[key]
            )
    if atoms.constraints:
        for a, b in zip(atoms.constraints, deserialized.constraints):
            assert repr(a) == repr(b)

    if deserialized.constraints:
        for a, b in zip(atoms.constraints, deserialized.constraints):
            assert repr(a) == repr(b)


def test_unsupported_type():
    """Test serialization of an unsupported type."""
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 2, 3])
    atoms.info["unsupported"] = fig  # Matplotlib figure is not supported
    with pytest.raises(TypeError):
        znjson.dumps(atoms, cls=znjson.ZnEncoder.from_converters([ASEConverter]))
