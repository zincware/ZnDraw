import dataclasses
import numpy as np
import networkx as nx
import typing as t

import ase
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data.colors import jmol_colors

from zndraw.bonds import ASEComputeBonds

@dataclasses.dataclass
class Frame:
    positions: t.Union[np.ndarray, list] = None
    cell: t.Union[np.ndarray, list] = None
    numbers: t.Union[np.ndarray, list, int] = None
    colors: t.Union[np.ndarray, list] = None
    radii: t.Union[np.ndarray, list] = None
    pbc: bool = False
    connectivity: nx.Graph() = None
    calc: t.Any = None

    @classmethod
    def from_atoms(cls, atoms: ase.Atoms):
        frame = cls(**atoms.arrays)

        frame.connectivity = atoms.connectivity

        try:
            calc_data = {}
            for key in atoms.calc.results:
                value = atoms.calc.results[key]
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                calc_data[key] = value

            frame.calc = calc_data
        except (RuntimeError, AttributeError):
            pass
        
        return frame

    def to_atoms(self) -> ase.Atoms:
        atoms = ase.Atoms(positions = self.positions, 
                          numbers = self.numbers, 
                          cell = None, #self.cell funktioniert wieso auch immer nicht. TODO
                          pbc = self.pbc)
        
        atoms.arrays["colors"] = self.colors
        
        atoms.connectivity = nx.Graph()
        for edge in self.connectivity:
            atoms.connectivity.add_edge(edge[0], edge[1], weight=edge[2])

        if self.calc is not None:
            atoms.calc = SinglePointCalculator(atoms)
            atoms.calc.results = {
                key: np.array(val) if isinstance(val, list) else val
                for key, val in self.calc.items()
            }

        return atoms

    def __len__(self):
        if isinstance(self.numbers, np.ndarray):
            return self.numbers.size
        elif isinstance(self.numbers, list):
            return len(self.numbers)
        elif isinstance(self.numbers, int):
            return 1