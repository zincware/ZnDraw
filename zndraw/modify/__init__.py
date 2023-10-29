import abc
import enum
import logging
import time
import typing as t

import ase
import numpy as np
from ase.data import chemical_symbols
from pydantic import BaseModel, ConfigDict, Field

log = logging.getLogger("zndraw")

Symbols = enum.Enum("Symbols", {symbol: symbol for symbol in chemical_symbols})


class UpdateScene(BaseModel, abc.ABC):
    @abc.abstractmethod
    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
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

    def run(self, vis) -> list[ase.Atoms]:
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
            vis.step += 1
        vis.selection = []


class Explode(UpdateScene):
    discriminator: t.Literal["Explode"] = Field("Explode")

    steps: int = Field(100, le=1000, ge=1)
    particles: int = Field(10, le=20, ge=1)
    delay: int = Field(0, le=60000, ge=0, description="Delay between each step in ms")

    def run(self, vis) -> list[ase.Atoms]:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atom_ids = vis.selection
        atoms = vis.atoms
        particles = []
        for _atom_id in atom_ids:
            for _ in range(self.particles):
                particles.append(ase.Atoms("Na", positions=[atoms.positions[_atom_id]]))

        for _ in range(self.steps):
            struct = atoms.copy()
            for particle in particles:
                particle.positions += np.random.normal(scale=0.1, size=(1, 3))
                struct += particle
            time.sleep(self.delay / 1000)
            vis.append(struct)
            vis.step += 1
        vis.selection = []


class Delete(UpdateScene):
    """Delete the selected atoms."""

    discriminator: t.Literal["Delete"] = Field("Delete")

    def run(self, vis) -> list[ase.Atoms]:
        atom_ids = vis.selection
        atoms = vis.atoms

        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]
        print(f"Deleting atoms {atom_ids}")
        vis.log(f"Deleting atoms {atom_ids}")
        for idx, atom_id in enumerate(sorted(atom_ids)):
            atoms.pop(atom_id - idx)  # we remove the atom and shift the index
        del atoms.connectivity
        vis.append(atoms)
        vis.selection = []
        vis.step += 1


class Move(UpdateScene):
    """Move the selected atoms along the line."""

    discriminator: t.Literal["Move"] = Field("Move")

    steps: int = Field(10, ge=1)

    def run(self, vis) -> list[ase.Atoms]:
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
            vis.step += 1
        vis.selection = []


class Duplicate(UpdateScene):
    discriminator: t.Literal["Duplicate"] = Field("Duplicate")

    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)
    symbol: Symbols

    def run(self, vis) -> list[ase.Atoms]:
        atoms = vis.atoms
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        for atom_id in vis.selection:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atom.symbol = self.symbol.name if self.symbol.name != "X" else atom.symbol
            atoms += atom
        vis.append(atoms)
        vis.step += 1
        vis.selection = []


class ChangeType(UpdateScene):
    discriminator: t.Literal["ChangeType"] = Field("ChangeType")

    symbol: Symbols

    def run(self, vis) -> list[ase.Atoms]:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        for atom_id in vis.selection:
            atoms[atom_id].symbol = self.symbol.name
        vis.append(atoms)
        vis.step += 1
        vis.selection = []


class AddLineParticles(UpdateScene):
    discriminator: t.Literal["AddLineParticles"] = Field("AddLineParticles")

    symbol: Symbols
    steps: int = Field(10, le=100, ge=1)

    def run(self, vis) -> list[ase.Atoms]:
        if len(vis) > vis.step + 1:
            del vis[vis.step + 1 :]

        atoms = vis.atoms
        for point in vis.points:
            atoms += ase.Atom(self.symbol.name, position=point)
        for _ in range(self.steps):
            vis.append(atoms)
            vis.step += 1


# class CustomModifier(UpdateScene):
#     discriminator: t.Literal["CustomModifier"] = Field("CustomModifier")

#     methods: t.Union[None, AddLineParticles, Rotate, Explode, Delete] = None


def get_modify_class(methods):
    class Modifier(UpdateScene):
        method: methods = Field(
            ..., description="Modify method", discriminator="discriminator"
        )

        model_config = ConfigDict(json_schema_extra=None)  # disable method hiding

        def run(self, vis) -> None:
            self.method.run(vis)

    return Modifier
