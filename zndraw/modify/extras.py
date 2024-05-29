import typing as t

import numpy as np
from pydantic import BaseModel, Field
from rdkit2ase import pack, smiles2atoms

from zndraw.settings import _MODIFY_FUNCTIONS

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw

_MODIFY_FUNCTIONS.append("zndraw.modify.extras.AddFromSMILES")
_MODIFY_FUNCTIONS.append("zndraw.modify.extras.Solvate")


class AddFromSMILES(BaseModel):
    """Place a molecule from a SMILES at all points."""

    discriminator: t.Literal["AddFromSMILES"] = Field("AddFromSMILES")
    SMILES: str = Field(..., description="SMILES string of the molecule to add")

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        molecule = smiles2atoms(self.SMILES)

        scene = vis.atoms

        points = vis.points
        if len(points) == 0:
            points = [np.array([0, 0, 0])]

        for point in points:
            molecule_copy = molecule.copy()
            molecule_copy.translate(point)
            scene.extend(molecule_copy)

        if hasattr(scene, "connectivity"):
            del scene.connectivity

        vis.append(scene)


class Solvate(BaseModel):
    """Solvate the current scene."""

    discriminator: t.Literal["Solvate"] = Field("Solvate")
    solvent: str = Field(..., description="Solvent to use (SMILES)")
    count: int = Field(10, ge=1, le=500, description="Number of solvent molecules to add")
    density: float = Field(789, description="Density of the solvent")
    pbc: bool = Field(True, description="Whether to use periodic boundary conditions")
    tolerance: float = Field(2.0, description="Tolerance for the solvent")

    def run(self, vis: "ZnDraw", **kwargs) -> None:
        scene = vis.atoms
        scene = pack(
            [(scene, 1), (self.solvent, self.count)],
            density=self.density,
            pbc=self.pbc,
            tolerance=self.tolerance,
        )
        vis.append(scene)
