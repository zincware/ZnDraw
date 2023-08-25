import enum
from typing import Union

import ase
import numpy as np
from ase.data import chemical_symbols
from pydantic import BaseModel, Field
from typing_extensions import Annotated

Symbols = enum.Enum("Symbols", {symbol: symbol for symbol in chemical_symbols})

class Explode(BaseModel):
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

class Duplicate(BaseModel):
    x: float = Field(0.5, le=5, ge=0)
    y: float = Field(0.5, le=5, ge=0)
    z: float = Field(0.5, le=5, ge=0)
    symbol: str

    def run(self, atom_ids: list[int], atoms: ase.Atoms, **kwargs) -> list[ase.Atoms]:
        for atom_id in atom_ids:
            atom = ase.Atom(atoms[atom_id].symbol, atoms[atom_id].position)
            atom.position += np.array([self.x, self.y, self.z])
            atom.symbol = self.symbol.name if self.symbol.name != "X" else atom.symbol
            atoms += atom
        return [atoms]
    
class Modifier(BaseModel):
    method: Annotated[Union[Duplicate, Explode], Field(alias="Method")] 