"""Extensions requiring the molify package."""

import typing as t

from pydantic import BaseModel, Field

from zndraw.extensions.abc import Category, Extension


class AddFromSMILES(Extension):
    """Add a molecule from SMILES notation."""

    category: t.ClassVar[Category] = Category.MODIFIER
    smiles: str = Field(
        ...,
        json_schema_extra={
            "x-custom-type": "smiles",
            "description": "SMILES notation for molecule",
        },
    )

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        import molify

        atoms = molify.smiles2atoms(self.smiles)
        vis.append(atoms)
        vis.log(f"""Added molecule {self.smiles}
```smiles
{self.smiles}
```""")
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

    category: t.ClassVar[Category] = Category.MODIFIER
    molecules: list[MoleculeSpec] = Field(
        default=[],
        json_schema_extra={"x-custom-type": "smiles-pack-array"},
    )
    density: float = Field(1.0, ge=0.0)

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        import molify

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
            f"at density {self.density} kg/m\u00b3"
        )


molify_modifiers: dict[str, type[Extension]] = {
    AddFromSMILES.__name__: AddFromSMILES,
    PackBox.__name__: PackBox,
}
