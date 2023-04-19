import ase
import numpy as np


def explode(atom_id: list[int], atoms: ase.Atoms) -> list[ase.Atoms]:
    particles = []
    for _atom_id in atom_id:
        for _ in range(5):
            particles.append(ase.Atoms("Na", positions=[atoms.positions[_atom_id]]))

    for _ in range(102):
        struct = atoms.copy()
        for particle in particles:
            particle.positions += np.random.normal(scale=0.1, size=(1, 3))
            struct += particle
        yield struct
