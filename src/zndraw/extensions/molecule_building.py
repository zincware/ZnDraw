"""Extensions requiring the molify package."""

import time
import typing as t

import molify
from pydantic import BaseModel, Field

from zndraw.extensions.abc import Category, Extension


class AddFromSMILES(Extension):
    """Add a molecule from SMILES notation."""

    category = Category.MODIFIER
    smiles: str = "CCO"

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        """Generate JSON schema with custom UI hints."""
        schema = super().model_json_schema(**kwargs)
        schema["properties"]["smiles"]["x-custom-type"] = "smiles"
        return schema

    def run(self, vis, **kwargs):
        atoms = molify.smiles2atoms(self.smiles)
        vis.append(atoms)
        vis.log(f"Added molecule from SMILES: {self.smiles}")
        vis.step = len(vis) - 1


class MoleculeSpec(BaseModel):
    """Specification for a molecule in PackBox."""

    smiles: str = Field(
        ...,
        json_schema_extra={
            "x-custom-type": "smiles",
            "description": "SMILES notation for molecule",
        },
    )
    count: int = Field(..., ge=1, description="Number of molecules")


class PackBox(Extension):
    """Pack molecules into a box at specified density."""

    category = Category.MODIFIER
    molecules: list[MoleculeSpec] = Field(
        default=[],
        json_schema_extra={"x-custom-type": "smiles-pack-array"},
    )
    density: float = Field(1.0, ge=0.0)

    def run(self, vis, **kwargs):
        conformers = [
            molify.smiles2conformers(mol.smiles, numConfs=mol.count)
            for mol in self.molecules
        ]
        box = molify.pack(
            data=conformers,
            counts=[mol.count for mol in self.molecules],
            density=self.density,
        )
        vis.append(box)
        vis.step = len(vis) - 1
        vis.log(
            f"Packed box with {len(self.molecules)} molecule types "
            f"at density {self.density} kg/m³"
        )


class Wait(Extension):
    """Wait for a specified time."""

    category = Category.MODIFIER
    time: float = 1.0

    def run(self, vis, **kwargs):
        with vis.progress_tracker(
            description=f"Waiting for {self.time} seconds"
        ) as tracker:
            for idx in range(100):
                time.sleep(self.time / 100)
                tracker.update(progress=idx + 1)


molify_modifiers: dict[str, type[Extension]] = {
    AddFromSMILES.__name__: AddFromSMILES,
    PackBox.__name__: PackBox,
}
