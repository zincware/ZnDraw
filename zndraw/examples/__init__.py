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
        for atom_id in atom_ids:
            atoms.pop(atom_id)
        return [atoms]
