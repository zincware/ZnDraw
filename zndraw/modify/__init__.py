import abc
import enum
import typing as t

import ase
import numpy as np
from ase.data import chemical_symbols
from pydantic import BaseModel, Field

Symbols = enum.Enum("Symbols", {symbol: symbol for symbol in chemical_symbols})


class UpdateScene(BaseModel, abc.ABC):
    @abc.abstractmethod
    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        pass


class Rotate(UpdateScene):
    """Rotate the selected atoms around a the line (2 points only)."""

    angle: float = Field(90, le=360, ge=0, description="Angle in degrees")
    direction: t.Literal["left", "right"] = Field(
        "left", description="Direction of rotation"
    )
    steps: int = Field(
        30, ge=1, description="Number of steps to take to complete the rotation"
    )

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        # split atoms object into the selected from atoms_ids and the remaining
        points = kwargs["points"]
        assert len(points) == 2

        angle = self.angle if self.direction == "left" else -self.angle
        angle = angle / self.steps

        atoms_selected = atoms[atom_ids]
        atoms_remaining = atoms[[x for x in range(len(atoms)) if x not in atom_ids]]
        # create a vector from the two points
        vector = points[1] - points[0]
        for _ in range(self.steps):
            # rotate the selected atoms around the vector
            atoms_selected.rotate(angle, vector, center=points[0])
            # merge the selected and remaining atoms
            atoms = atoms_selected + atoms_remaining
            yield atoms


class Explode(UpdateScene):
    steps: int = Field(100, le=1000, ge=1)
    particles: int = Field(10, le=20, ge=1)

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        particles = []
        for _atom_id in atom_ids:
            for _ in range(self.particles):
                particles.append(ase.Atoms("Na", positions=[atoms.positions[_atom_id]]))

        for _ in range(self.steps):
            struct = atoms.copy()
            for particle in particles:
                particle.positions += np.random.normal(scale=0.1, size=(1, 3))
                struct += particle
            yield struct


class Delete(UpdateScene):
    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        for idx, atom_id in enumerate(sorted(atom_ids)):
            atoms.pop(atom_id - idx)  # we remove the atom and shift the index
        return [atoms]


class Move(UpdateScene):
    """Move the selected atoms along the line."""

    steps: int = Field(10, ge=1)

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        segments = kwargs["segments"]
        atoms_selected = atoms[atom_ids]
        atoms_remaining = atoms[[x for x in range(len(atoms)) if x not in atom_ids]]
        if self.steps > len(segments):
            raise ValueError(
                "The number of steps must be less than the number of segments. You can add more points to increase the number of segments."
            )

        # atoms_selected.positions = segments[0]
        # yield atoms_selected + atoms_remaining

        for idx in range(1, self.steps):
            # get the vector between the two points

            start_idx = int((idx - 1) * len(segments) / self.steps)
            end_idx = int(idx * len(segments) / self.steps)

            vector = segments[end_idx] - segments[start_idx]
            # move the selected atoms along the vector
            atoms_selected.positions += vector
            # merge the selected and remaining atoms
            atoms = atoms_selected + atoms_remaining
            yield atoms


class Duplicate(UpdateScene):
    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)
    symbol: Symbols

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        for atom_id in atom_ids:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atom.symbol = self.symbol.name if self.symbol.name != "X" else atom.symbol
            atoms += atom
        return [atoms]


class ChangeType(UpdateScene):
    symbol: Symbols

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        for atom_id in atom_ids:
            atoms[atom_id].symbol = self.symbol.name
        return [atoms]


class AddLineParticles(UpdateScene):
    symbol: Symbols
    steps: int = Field(10, le=100, ge=1)

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        points = kwargs["points"]
        for point in points:
            atoms += ase.Atom(self.symbol.name, position=point)
        for _ in range(self.steps):
            yield atoms


class Demo(UpdateScene):
    """Scene update for testing purposes."""

    z: float = Field(0.5, le=5, ge=0)
    symbol: Symbols

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        for atom_id in atom_ids:
            atoms[atom_id].symbol = self.symbol.name
        return [atoms]
