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

    method: t.Literal["Rotate"] = Field("Rotate")

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

        atoms_selected, atoms_remaining = self.apply_selection(atom_ids, atoms)
        # create a vector from the two points
        vector = points[1] - points[0]
        for _ in range(self.steps):
            # rotate the selected atoms around the vector
            atoms_selected.rotate(angle, vector, center=points[0])
            # merge the selected and remaining atoms
            atoms = atoms_selected + atoms_remaining
            yield atoms


class Explode(UpdateScene):
    method: t.Literal["Explode"] = Field("Explode")

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
    """Delete the selected atoms."""

    method: t.Literal["Delete"] = Field("Delete")

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        for idx, atom_id in enumerate(sorted(atom_ids)):
            atoms.pop(atom_id - idx)  # we remove the atom and shift the index
        del atoms.connectivity
        return [atoms]


class Move(UpdateScene):
    """Move the selected atoms along the line."""

    method: t.Literal["Move"] = Field("Move")

    steps: int = Field(10, ge=1)

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        segments = kwargs["segments"]
        atoms_selected, atoms_remaining = self.apply_selection(atom_ids, atoms)

        if self.steps > len(segments):
            raise ValueError(
                "The number of steps must be less than the number of segments. You can add more points to increase the number of segments."
            )

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
    method: t.Literal["Duplicate"] = Field("Duplicate")

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
    method: t.Literal["ChangeType"] = Field("ChangeType")

    symbol: Symbols

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        for atom_id in atom_ids:
            atoms[atom_id].symbol = self.symbol.name
        return [atoms]


class AddLineParticles(UpdateScene):
    method: t.Literal["AddLineParticles"] = Field("AddLineParticles")

    symbol: Symbols
    steps: int = Field(10, le=100, ge=1)

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        points = kwargs["points"]
        for point in points:
            atoms += ase.Atom(self.symbol.name, position=point)
        for _ in range(self.steps):
            yield atoms


def get_modify_class(methods):
    class Modifier(UpdateScene):
        method: methods = Field(
            ..., description="Modify method", discriminator="method"
        )

        def run(self, *args, **kwargs) -> list[ase.Atoms]:
            return self.method.run(*args, **kwargs)

        @classmethod
        def model_json_schema(cls, *args, **kwargs) -> dict[str, t.Any]:
            schema = super().model_json_schema(*args, **kwargs)
            for prop in [x.__name__ for x in t.get_args(methods)]:
                schema["$defs"][prop]["properties"]["method"]["options"] = {
                    "hidden": True
                }
                schema["$defs"][prop]["properties"]["method"]["type"] = "string"

            return schema

    return Modifier
