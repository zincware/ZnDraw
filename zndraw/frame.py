import copy
import dataclasses
import typing as t

import ase
import networkx as nx
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data.colors import jmol_colors
from ase.neighborlist import natural_cutoffs
from pydantic import Field

from zndraw.bonds import ASEComputeBonds
from zndraw.data import _get_radius, _rgb2hex


@dataclasses.dataclass
class Frame:
    # Main Attributes
    positions: t.Union[np.ndarray, list] = None
    cell: t.Union[np.ndarray, list] = np.array([0.0, 0.0, 0.0])
    numbers: t.Union[np.ndarray, list, int] = None
    colors: t.Union[np.ndarray, list] = None
    radii: t.Union[np.ndarray, list] = None
    pbc: bool = False
    connectivity: nx.Graph() = nx.Graph()
    calc: t.Any = None

    vecField: t.Any = np.array(["TestArray", "0"])
    # Secondary Attributes
    # define calculations in this class
    bonds: bool = True
    auto_bonds: bool = True

    def __post_init__(self):
        if isinstance(self.positions, list):
            self.positions = np.array(self.positions)
        if isinstance(self.positions, list):
            self.cell = np.array(self.cell)
        if isinstance(self.positions, list):
            self.numbers = np.array(self.numbers)
        if isinstance(self.positions, list):
            self.colors = np.array(self.colors)
        if isinstance(self.positions, list):
            self.radii = np.array(self.radii)
        # i know its ugly, but it works
        # if i convert them without the check, the speed drops dramatically

    @classmethod
    def from_atoms(cls, atoms: ase.Atoms):
        frame = cls(**atoms.arrays)

        frame.cell = atoms.cell
        frame.pbc = atoms.pbc

        if hasattr(atoms, "connectivity"):
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
        atoms = ase.Atoms(
            positions=self.positions,
            numbers=self.numbers,
            cell=self.cell,  # self.cell funktioniert wieso auch immer nicht. TODO
            pbc=self.pbc,
        )

        # atoms.arrays["colors"] = self.colors

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

    def calc_bonds(self):
        """
        Experimental Function, currently not in use
        """
        single_bond_multiplier: float = Field(1.2, le=2, ge=0)
        double_bond_multiplier: float = Field(0.9, le=1, ge=0)
        triple_bond_multiplier: float = Field(0.0, le=1, ge=0)

        cutoffs = [
            single_bond_multiplier,
            double_bond_multiplier,
            triple_bond_multiplier,
        ]
        frame_copy = copy.deepcopy(self)
        connectivity_matrix = np.zeros((len(self), len(self)), dtype=int)
        distance_matrix = self.dist_matrix()
        np.fill_diagonal(distance_matrix, np.inf)
        for cutoff in cutoffs:
            cutoffs = np.array(natural_cutoffs(frame_copy, mult=cutoff))
            cutoffs = cutoffs[:, None] + cutoffs[None, :]
            connectivity_matrix[distance_matrix <= cutoffs] += 1
        self.connectivity = nx.from_numpy_array(connectivity_matrix)

    def dist_matrix(self):
        """
        Experimental function, currently not in use
        """
        matrix = np.zeros((len(self), len(self)))
        for i in range(1, len(self)):
            for j in range(i, len(self)):
                matrix[i, j] = np.linalg.norm(self.positions[i] - self.positions[j])
        return matrix

    def get_bonds(self):
        bonds = []
        for edge in self.connectivity.edges:
            bonds.append((edge[0], edge[1], self.connectivity.edges[edge]["weight"]))
        return bonds

    def __len__(self):
        if isinstance(self.numbers, np.ndarray):
            return self.numbers.size
        elif isinstance(self.numbers, list):
            return len(self.numbers)
        elif isinstance(self.numbers, int):
            return 1

    def frame_to_json(self):
        frame_dict = self.__dict__

        for key, value in frame_dict.items():
            if isinstance(value, np.ndarray):
                frame_dict[key] = value.tolist()
            elif isinstance(value, np.generic):
                frame_dict[key] = value.item()

        if frame_dict["colors"] is None:
            frame_dict["colors"] = [
                _rgb2hex(jmol_colors[number]) for number in frame_dict["numbers"]
            ]

        frame_dict["radii"] = [_get_radius(number) for number in frame_dict["numbers"]]

        if self.bonds:
            try:
                if self.auto_bonds:
                    ase_bond_calculator = ASEComputeBonds()
                    self.connectivity = ase_bond_calculator.build_graph(self.to_atoms())
                frame_dict["connectivity"] = self.get_bonds()
            except AttributeError:
                frame_dict["connectivity"] = []
        else:
            frame_dict["connectivity"] = []

        return frame_dict

    @classmethod
    def frame_from_json(cls, data):
        frame = cls(
            positions=np.array(data["positions"]),
            cell=np.array(data["cell"]),
            numbers=np.array(data["numbers"]),
            colors=np.array(data["colors"]),
            radii=np.array(data["radii"]),
            pbc=data["pbc"],
            calc=data["calc"],
        )

        if (
            "vecField" in data
        ):  # currently there is the vecField part in js missing. So this is useless at the moment
            frame.vecField = data["vecField"]

        if "connectivity" in data:
            frame.connectivity = nx.Graph()
            for edge in data["connectivity"]:
                frame.connectivity.add_edge(edge[0], edge[1], weight=edge[2])

        return frame
