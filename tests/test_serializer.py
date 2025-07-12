import json

import ase
import numpy.testing as npt
import znjson
from ase.calculators.singlepoint import SinglePointCalculator
from ase.constraints import FixAtoms

from zndraw.converter import ASEConverter


def test_ase_converter(s22):
    s22[0].info["connectivity"] = [[0, 1, 1], [1, 2, 1], [2, 3, 1]]
    s22[3].calc = SinglePointCalculator(s22[3])
    s22[3].calc.results = {"energy": 0.0, "predicted_energy": 1.0}
    s22[4].info = {"key": "value"}

    structures_json = znjson.dumps(
        s22, cls=znjson.ZnEncoder.from_converters([ASEConverter])
    )

    non_json = json.loads(structures_json)
    assert "numbers" not in non_json[0]["value"]["arrays"]
    assert "positions" not in non_json[0]["value"]["arrays"]
    assert "pbc" not in non_json[0]["value"]["info"]
    assert "cell" not in non_json[0]["value"]["info"]

    structures = znjson.loads(
        structures_json, cls=znjson.ZnDecoder.from_converters([ASEConverter])
    )
    for s1, s2 in zip(s22, structures):
        assert s1 == s2

    assert structures[0].info["connectivity"] == [[0, 1, 1], [1, 2, 1], [2, 3, 1]]
    assert "connectivity" not in structures[1].info

    assert structures[3].calc.results == {"energy": 0.0, "predicted_energy": 1.0}

    assert "colors" not in structures[0].arrays
    assert "radii" not in structures[0].arrays

    assert structures[4].info == {"key": "value"}


def test_exotic_atoms():
    atoms = ase.Atoms("X", positions=[[0, 0, 0]])
    atoms.arrays["colors"] = ["#ff0000"]
    atoms.arrays["radii"] = [0.3]

    new_atoms = znjson.loads(
        znjson.dumps(atoms, cls=znjson.ZnEncoder.from_converters([ASEConverter])),
        cls=znjson.ZnDecoder.from_converters([ASEConverter]),
    )
    npt.assert_array_equal(new_atoms.arrays["colors"], ["#ff0000"])
    npt.assert_array_equal(new_atoms.arrays["radii"], [0.3])


def test_modified_atoms():
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    new_atoms = znjson.loads(
        znjson.dumps(atoms, cls=znjson.ZnEncoder.from_converters([ASEConverter])),
        cls=znjson.ZnDecoder.from_converters([ASEConverter]),
    )
    npt.assert_array_equal(new_atoms.get_atomic_numbers(), [1, 1])

    # subtract
    atoms = new_atoms[:1]
    new_atoms = znjson.loads(
        znjson.dumps(atoms, cls=znjson.ZnEncoder.from_converters([ASEConverter])),
        cls=znjson.ZnDecoder.from_converters([ASEConverter]),
    )

    npt.assert_array_equal(new_atoms.get_atomic_numbers(), [1])

    # add
    atoms = new_atoms + ase.Atoms("H", positions=[[0, 0, 1]])

    new_atoms = znjson.loads(
        znjson.dumps(atoms, cls=znjson.ZnEncoder.from_converters([ASEConverter])),
        cls=znjson.ZnDecoder.from_converters([ASEConverter]),
    )

    npt.assert_array_equal(new_atoms.get_atomic_numbers(), [1, 1])


def test_constraints_fixed_atoms():
    atoms = ase.Atoms("H2", positions=[[0, 0, 0], [0, 0, 1]])
    atoms.set_constraint(FixAtoms([0]))
    new_atoms = znjson.loads(
        znjson.dumps(atoms, cls=znjson.ZnEncoder.from_converters([ASEConverter])),
        cls=znjson.ZnDecoder.from_converters([ASEConverter]),
    )
    assert isinstance(new_atoms.constraints[0], FixAtoms)
    assert new_atoms.constraints[0].index == [0]
