import abc
import enum
import logging
import typing as t

import ase
import numpy as np
from ase.data import chemical_symbols
from pydantic import BaseModel, ConfigDict, Field

log = logging.getLogger("zndraw")

Symbols = enum.Enum("Symbols", {symbol: symbol for symbol in chemical_symbols})

if t.TYPE_CHECKING:
    from zndraw.zndraw import ZnDraw


class UpdateScene(BaseModel, abc.ABC):
    @abc.abstractmethod
    def run(self, vis: "ZnDraw") -> None:
        """Method called when running the modifier."""
        pass

    def apply_selection(
        self, atom_ids: list[int], atoms: ase.Atoms
    ) -> t.Tuple[ase.Atoms, ase.Atoms]:
        """Split the atoms object into the selected and remaining atoms."""
        atoms_selected = atoms[atom_ids]
        atoms_remaining_ids = [x for x in range(len(atoms)) if x not in atom_ids]
        if len(atoms_remaining_ids) > 0:
            atoms_remaining = atoms[atoms_remaining_ids]
        else:
            atoms_remaining = ase.Atoms()
        return atoms_selected, atoms_remaining


class Rotate(UpdateScene):
    """Rotate the selected atoms around a the line (2 points only)."""

    discriminator: t.Literal["Rotate"] = Field("Rotate")

    angle: float = Field(90, le=360, ge=0, description="Angle in degrees")
    direction: t.Literal["left", "right"] = Field(
        "left", description="Direction of rotation"
    )
    steps: int = Field(
        30, ge=1, description="Number of steps to take to complete the rotation"
    )

    def run(self, vis: "ZnDraw") -> None:
        # split atoms object into the selected from atoms_ids and the remaining
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        points = vis.points
        atom_ids = vis.selection
        atoms = vis.atoms
        assert len(points) == 2

        angle = self.angle if self.direction == "left" else -self.angle
        angle = angle / self.steps

        atoms_selected, atoms_remaining = self.apply_selection(atom_ids, atoms)
        # create a vector from the two points
        vector = points[1] - points[0]
        for _ in range(self.steps):
            # rotate the selected atoms around the vector
            atoms_selected.rotate(angle, vector, center=points[0])
            # merge the selected and remaining atoms
            atoms = atoms_selected + atoms_remaining
            vis.append(atoms)
        vis.selection = []


class Delete(UpdateScene):
    """Delete the selected atoms."""

    discriminator: t.Literal["Delete"] = Field("Delete")

    def run(self, vis: "ZnDraw") -> None:
        atom_ids = vis.selection
        atoms = vis.atoms

        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]
        vis.log(f"Deleting atoms {atom_ids}")
        for idx, atom_id in enumerate(sorted(atom_ids)):
            atoms.pop(atom_id - idx)  # we remove the atom and shift the index
        del atoms.connectivity
        vis.append(atoms)
        vis.selection = []


class Move(UpdateScene):
    """Move the selected atoms along the line."""

    discriminator: t.Literal["Move"] = Field("Move")

    steps: int = Field(10, ge=1)

    def run(self, vis: "ZnDraw") -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        atoms_selected, atoms_remaining = self.apply_selection(vis.selection, atoms)
        if self.steps > len(vis.segments):
            raise ValueError(
                "The number of steps must be less than the number of segments. You can add more points to increase the number of segments."
            )

        for idx in range(1, self.steps):
            # get the vector between the two points
            start_idx = int((idx - 1) * len(vis.segments) / self.steps)
            end_idx = int(idx * len(vis.segments) / self.steps)

            vector = vis.segments[end_idx] - vis.segments[start_idx]
            # move the selected atoms along the vector
            atoms_selected.positions += vector
            # merge the selected and remaining atoms
            atoms = atoms_selected + atoms_remaining
            vis.append(atoms)
        vis.selection = []


class Duplicate(UpdateScene):
    discriminator: t.Literal["Duplicate"] = Field("Duplicate")

    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)
    symbol: Symbols

    def run(self, vis: "ZnDraw") -> None:
        atoms = vis.atoms
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        for atom_id in vis.selection:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atom.symbol = self.symbol.name if self.symbol.name != "X" else atom.symbol
            atoms += atom

        vis.append(atoms)
        vis.selection = []


class ChangeType(UpdateScene):
    discriminator: t.Literal["ChangeType"] = Field("ChangeType")

    symbol: Symbols

    def run(self, vis: "ZnDraw") -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        for atom_id in vis.selection:
            atoms[atom_id].symbol = self.symbol.name

        vis.append(atoms)
        vis.selection = []


class AddLineParticles(UpdateScene):
    discriminator: t.Literal["AddLineParticles"] = Field("AddLineParticles")

    symbol: Symbols
    steps: int = Field(10, le=100, ge=1)

    def run(self, vis: "ZnDraw") -> None:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        for point in vis.points:
            atoms += ase.Atom(self.symbol.name, position=point)

        for _ in range(self.steps):
            vis.append(atoms)


class Wrap(UpdateScene):
    """Wrap the atoms to the cell."""

    discriminator: t.Literal["Wrap"] = Field("Wrap")
    recompute_bonds: bool = True

    def run(self, vis: "ZnDraw") -> None:
        vis.log("Downloading atoms...")
        atoms_list = list(vis)
        vis.step = 0

        del vis[1:]

        for idx, atoms in enumerate(atoms_list):
            atoms.wrap()
            if self.recompute_bonds:
                delattr(atoms, "connectivity")
            vis[idx] = atoms


class Center(UpdateScene):
    """Move the atoms, such that the selected atom is in the center of the cell."""

    discriminator: t.Literal["Center"] = Field("Center")
    recompute_bonds: bool = True
    dynamic: bool = Field(
        False, description="Move the atoms to the center of the cell at each step"
    )
    wrap: bool = Field(True, description="Wrap the atoms to the cell")

    def run(self, vis: "ZnDraw") -> None:
        selection = vis.selection
        if len(selection) != 1:
            vis.log("Please select exactly one atom to center on.")
            return

        vis.log("Downloading atoms...")
        atoms_list = list(vis)

        if not self.dynamic:
            center = atoms_list[vis.step][selection[0]].position
        else:
            center = None

        vis.step = 0

        del vis[1:]

        for idx, atoms in enumerate(atoms_list):
            if self.dynamic:
                center = atoms[selection[0]].position
            atoms.positions -= center
            atoms.positions += np.diag(atoms.cell) / 2
            if self.wrap:
                atoms.wrap()
            if self.recompute_bonds:
                delattr(atoms, "connectivity")

            vis[idx] = atoms


class Replicate(UpdateScene):
    discriminator: t.Literal["Replicate"] = Field("Replicate")
    x: int = Field(2, ge=1)
    y: int = Field(2, ge=1)
    z: int = Field(2, ge=1)

    keep_box: bool = Field(False, description="Keep the original box size")

    def run(self, vis: "ZnDraw") -> None:
        vis.log("Downloading atoms...")
        atoms_list = list(vis)
        vis.step = 0

        del vis[1:]

        for idx, atoms in enumerate(atoms_list):
            atoms = atoms.repeat((self.x, self.y, self.z))
            if self.keep_box:
                atoms.cell = atoms_list[idx].cell
            vis[idx] = atoms


# class CustomModifier(UpdateScene):
#     discriminator: t.Literal["CustomModifier"] = Field("CustomModifier")

#     methods: t.Union[None, AddLineParticles, Rotate, Explode, Delete] = None


def get_modify_class(methods):
    class Modifier(UpdateScene):
        method: methods = Field(
            ..., description="Modify method", discriminator="discriminator"
        )

        model_config = ConfigDict(json_schema_extra=None)  # disable method hiding

        def run(self, vis, **kwargs) -> None:
            self.method.run(vis, **kwargs)

    return Modifier
