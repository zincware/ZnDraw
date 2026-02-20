"""Modifier extensions for manipulating atoms in the scene.

These extensions modify the atomic structure (add, delete, move atoms, etc.).
"""

import enum
import typing as t

from ase.data import chemical_symbols
from pydantic import Field

from zndraw.extensions.abc import Category, Extension
from zndraw.extensions.molecule_building import molify_modifiers

# Create an enum of all chemical symbols
Symbols = enum.Enum("Symbols", {symbol: symbol for symbol in chemical_symbols})


class UpdateScene(Extension):
    """Base class for all modifier extensions."""

    category: t.ClassVar[Category] = Category.MODIFIER


class Delete(UpdateScene):
    """Delete the selected atoms."""

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        import ase

        atom_ids = vis.selection
        atoms = vis[vis.step]

        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]
        if len(atom_ids) == len(atoms):
            vis.append(ase.Atoms())
        else:
            for idx, atom_id in enumerate(sorted(atom_ids)):
                atoms.pop(atom_id - idx)
            atoms.info.pop("connectivity", None)
        vis.append(atoms)
        vis.selection = []
        vis.step += 1


class Duplicate(UpdateScene):
    """Duplicate selected atoms with an offset."""

    x: float = Field(0.5, le=5, ge=0, description="X offset")
    y: float = Field(0.5, le=5, ge=0, description="Y offset")
    z: float = Field(0.5, le=5, ge=0, description="Z offset")
    symbol: Symbols = Field(Symbols.X, description="Symbol of the new atoms")

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        import ase
        import numpy as np

        atoms = vis.atoms
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        for atom_id in vis.selection:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atom.symbol = self.symbol.name if self.symbol.name != "X" else atom.symbol
            atoms += atom
            atoms.arrays.pop("colors", None)
            atoms.arrays.pop("radii", None)
            atoms.info.pop("connectivity", None)

        vis.append(atoms)
        vis.selection = []


class ChangeType(UpdateScene):
    """Change the type of selected atoms."""

    symbol: Symbols = Field(..., description="New symbol for selected atoms")

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms.copy()
        for atom_id in vis.selection:
            atoms[atom_id].symbol = self.symbol.name

        atoms.arrays.pop("colors", None)
        atoms.arrays.pop("radii", None)
        atoms.info.pop("connectivity", None)

        vis.append(atoms)
        vis.selection = []


class Wrap(UpdateScene):
    """Wrap atoms to the periodic cell."""

    recompute_bonds: bool = Field(True, description="Recompute bonds after wrapping")
    all: bool = Field(False, description="Apply to the full trajectory")

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        if self.all:
            for idx, atoms in enumerate(vis):
                atoms.wrap()
                if self.recompute_bonds:
                    atoms.info.pop("connectivity", None)
                vis[idx] = atoms
        else:
            atoms = vis.atoms
            atoms.wrap()
            if self.recompute_bonds:
                atoms.info.pop("connectivity", None)
            vis[vis.step] = atoms


class Center(UpdateScene):
    """Center atoms around the selected atom(s) center of mass."""

    recompute_bonds: bool = Field(True, description="Recompute bonds after centering")
    dynamic: bool = Field(
        False, description="Recalculate center at each step for trajectories"
    )
    wrap: bool = Field(True, description="Wrap atoms to the cell after centering")
    all: bool = Field(False, description="Apply to the full trajectory")

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        import numpy as np

        selection = vis.selection
        if len(selection) < 1:
            vis.log("Please select at least one atom.")
            raise ValueError("No atoms selected for centering.")

        if not self.dynamic:
            center = vis.atoms[selection].get_center_of_mass()
        else:
            center = None

        if self.all:
            for idx, atoms in enumerate(vis):
                if self.dynamic:
                    center = atoms[selection].get_center_of_mass()
                atoms.positions -= center
                atoms.positions += np.diag(atoms.cell) / 2
                if self.wrap:
                    atoms.wrap()
                if self.recompute_bonds:
                    atoms.info.pop("connectivity", None)
                vis[idx] = atoms
        else:
            atoms = vis.atoms
            center = atoms[selection].get_center_of_mass()
            atoms.positions -= center
            atoms.positions += np.diag(atoms.cell) / 2
            if self.wrap:
                atoms.wrap()
            if self.recompute_bonds:
                atoms.info.pop("connectivity", None)
            vis[vis.step] = atoms


class Replicate(UpdateScene):
    """Replicate the unit cell."""

    x: int = Field(2, ge=1, le=10, json_schema_extra={"format": "range"})
    y: int = Field(2, ge=1, le=10, json_schema_extra={"format": "range"})
    z: int = Field(2, ge=1, le=10, json_schema_extra={"format": "range"})
    keep_box: bool = Field(False, description="Keep the original box size")
    all: bool = Field(False, description="Apply to the full trajectory")

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        if self.all:
            for idx, atoms in enumerate(vis):
                original_cell = atoms.cell.copy()
                atoms = atoms.repeat((self.x, self.y, self.z))
                if self.keep_box:
                    atoms.cell = original_cell
                vis[idx] = atoms
        else:
            atoms = vis.atoms
            original_cell = atoms.cell.copy()
            atoms = atoms.repeat((self.x, self.y, self.z))
            if self.keep_box:
                atoms.cell = original_cell
            vis[vis.step] = atoms


class NewCanvas(UpdateScene):
    """Clear the scene and start fresh."""

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        import ase

        from zndraw.geometries import Curve, Plane

        # Clear drawing curve
        vis.geometries["curve"] = Curve(position=[])
        vis.append(ase.Atoms())
        vis.selection = []
        step = len(vis) - 1
        vis.step = step
        vis.bookmarks[step] = "New Scene"
        vis.geometries["plane"] = Plane(
            position=[(0.0, 0.0, 0.0)],
            rotation=[(0.0, 0.0, 0.0)],
            scale=[(1.0, 1.0, 1.0)],
            size=[(10.0, 10.0)],
        )


class RemoveAtoms(UpdateScene):
    """Remove the current frame from the trajectory."""

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        del vis[vis.step]


class Empty(UpdateScene):
    """Add an empty frame (no atoms)."""

    def run(self, vis: t.Any, **kwargs: t.Any) -> None:
        import ase

        vis.append(ase.Atoms())


# Registry of all modifier extensions
modifiers: dict[str, type[Extension]] = {
    Delete.__name__: Delete,
    Duplicate.__name__: Duplicate,
    ChangeType.__name__: ChangeType,
    Wrap.__name__: Wrap,
    Center.__name__: Center,
    Replicate.__name__: Replicate,
    NewCanvas.__name__: NewCanvas,
    RemoveAtoms.__name__: RemoveAtoms,
    Empty.__name__: Empty,
}

modifiers.update(molify_modifiers)
