import abc

import ase
import numpy as np
from pydantic import BaseModel, Field


class UpdateFunction(BaseModel, abc.ABC):
    @abc.abstractmethod
    def run(self, atom_ids: list[int], atoms: ase.Atoms) -> list[ase.Atoms]:
        pass


class Explode(UpdateFunction):
    steps: int = Field(100, le=1000, ge=1)
    particles: int = Field(10, le=100, ge=1)

    def run(self, atom_ids: list[int], atoms: ase.Atoms) -> list[ase.Atoms]:
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


class Delete(UpdateFunction):
    def run(self, atom_ids: list[int], atoms: ase.Atoms) -> list[ase.Atoms]:
        for idx, atom_id in enumerate(sorted(atom_ids)):
            atoms.pop(atom_id - idx)  # we remove the atom and shift the index
        return [atoms]


class Move(UpdateFunction):
    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)

    def run(self, atom_ids: list[int], atoms: ase.Atoms) -> list[ase.Atoms]:
        for atom_id in atom_ids:
            atom = atoms[atom_id]
            atom.position += np.array([self.x, self.y, self.z])
            atoms += atom
        return [atoms]


class Duplicate(UpdateFunction):
    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)

    def run(self, atom_ids: list[int], atoms: ase.Atoms) -> list[ase.Atoms]:
        for atom_id in atom_ids:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atoms += atom
        return [atoms]


class ChangeType(UpdateFunction):
    symbol: str = Field("")

    def run(self, atom_ids: list[int], atoms: ase.Atoms) -> list[ase.Atoms]:
        raise NotImplementedError("Chaning Atom Type is currently not supported")
        for atom_id in atom_ids:
            atom = atoms.pop(atom_id)
            atom.symbol = self.symbol
            atoms += atom
        print(atoms)
        return [atoms]
