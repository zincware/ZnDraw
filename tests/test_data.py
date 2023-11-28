import pytest
from ase.build import bulk, molecule

from zndraw.frame import Frame


@pytest.mark.parametrize(
    "atoms", [(molecule("C6H6")), (molecule("H2O")), (bulk("Cu", "fcc", a=3.6))]
)
def test_atoms_from_and_to_json_return_unchanged_atoms_object(atoms):
    frame = Frame.from_atoms(atoms)
    assert frame == Frame(positions=atoms.positions,
                          numbers=atoms.numbers,
                          cell = atoms.cell,
                          pbc = atoms.pbc)
    assert frame.to_atoms() == atoms
