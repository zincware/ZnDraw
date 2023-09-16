import pytest
from ase.build import bulk, molecule

from zndraw.data import atoms_from_json, atoms_to_json


@pytest.mark.parametrize(
    "atoms", [(molecule("C6H6")), (molecule("H2O")), (bulk("Cu", "fcc", a=3.6))]
)
def test_atoms_from_and_to_json_return_unchanged_atoms_object(atoms):
    atom_copy = atoms.copy()
    assert atom_copy == atoms_from_json(atoms_to_json(atoms))
