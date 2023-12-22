import dataclasses

import ase
import networkx as nx
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data.colors import jmol_colors

from zndraw.bonds import ASEComputeBonds
from zndraw.utils import get_radius, rgb2hex


@dataclasses.dataclass
class Frame:
    """
    Primary Attributes:
    These attributes directly contain data of your frame
    that will be displayed in the visualizer.
    -------------
    positions : np.ndarray
        contains the positions of each particle in 3 dimensions
    cell : np.ndarray
        contains the cell size of the frame
    numbers : t.Union[np.ndarray, list, int] #TODO update in case int gets removed
        contains the number of each individual atom
    colors : np.ndarray
        contains the hexadecimal color representation of each atom.
    radii : np.ndarray
        contains the radius of each atom that is displayed in the visualizer.
    pbc : bool
        determines periodic boundary conditions
    connectivity : nx.Graph()
        contains the bonds between singular atoms
    calc : dict
        contains properties of the frame, such as energy,
        that can be viewed using the analyze function.
    vector_field : dict
        WIP: contains a flowfield that will be displayed in the simulation box.
        will contain box-length, numbers of vectors per dimension and the directional
        vectors themself.#
    -------------

    Secondary Attributes:
    These attributes influence the usage of the primary attributes, such as if
    bonds are displayed or in what way bonds are calculated.
    -------------
    bonds : bool
        determines if bonds are drawn
    auto_bonds : bool
        if true uses module ase to calculate chemically accurate bonds of atoms
        using the positions and numbers (e.g. number = 1 = Hydrogen)
    """

    positions: np.ndarray = None
    cell: np.ndarray = dataclasses.field(default_factory=lambda: np.array([0.0, 0.0, 0.0]))
    numbers: np.ndarray = None
    colors: np.ndarray = None
    radii: np.ndarray = None
    momenta: np.ndarray = None
    forces: np.ndarray = None
    pbc: bool = False
    connectivity: nx.Graph() = dataclasses.field(default_factory=nx.empty_graph)
    calc: dict = None
    vector_field: dict = None

    bonds: bool = True
    auto_bonds: bool = True

    def __post_init__(self):
        """
        Converts all lists to np.ndarray
        """
        for item in ["positions", "numbers", "colors", "radii"]:
            if isinstance(getattr(self, item), list):
                setattr(self, item, np.array(getattr(self, item)))

        if not isinstance(self.cell, np.ndarray):
            self.cell = np.array(self.cell)

    @classmethod
    def from_atoms(cls, atoms: ase.Atoms):
        """
        Creates an instance of the frame class from an ase.Atoms object
        """
        data = atoms.arrays

        frame = cls(positions=data["positions"], numbers=data["numbers"])

        frame.cell = np.array(atoms.cell)
        frame.pbc = atoms.pbc

        if hasattr(atoms, "connectivity"):
            frame.connectivity = atoms.connectivity

        try:
            calc_data = {}
            for key, value in atoms.calc.results.items():
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                calc_data[key] = value
            frame.calc = calc_data
        except (
            RuntimeError,
            AttributeError,
        ):  # This exception happens, when there is no calc-attribute given.
            pass

        return frame

    def to_atoms(self) -> ase.Atoms:
        """
        Creates an ase.Atoms object from a Frame instance
        """
        atoms = ase.Atoms(
            positions=self.positions, numbers=self.numbers, cell=self.cell, pbc=self.pbc
        )

        # atoms.arrays["colors"] = self.colors # TODO: see https://github.com/zincware/ZnDraw/issues/279

        atoms.connectivity = self.connectivity

        if self.calc is not None:
            atoms.calc = SinglePointCalculator(atoms)
            atoms.calc.results = {
                key: np.array(val) if isinstance(val, list) else val
                for key, val in self.calc.items()
            }

        return atoms

    # TODO: use instead of ASEComputeBonds() in to_dict()
    # def calc_bonds(self):
    #     """
    #     Experimental Function, currently not in use
    #     """
    #     single_bond_multiplier: float = Field(1.2, le=2, ge=0)
    #     double_bond_multiplier: float = Field(0.9, le=1, ge=0)
    #     triple_bond_multiplier: float = Field(0.0, le=1, ge=0)

    #     cutoffs = [
    #         single_bond_multiplier,
    #         double_bond_multiplier,
    #         triple_bond_multiplier,
    #     ]
    #     frame_copy = copy.deepcopy(self)
    #     connectivity_matrix = np.zeros((len(self), len(self)), dtype=int)
    #     distance_matrix = self.dist_matrix()
    #     np.fill_diagonal(distance_matrix, np.inf)
    #     for cutoff in cutoffs:
    #         cutoffs = np.array(natural_cutoffs(frame_copy, mult=cutoff))
    #         cutoffs = cutoffs[:, None] + cutoffs[None, :]
    #         connectivity_matrix[distance_matrix <= cutoffs] += 1
    #     self.connectivity = nx.from_numpy_array(connectivity_matrix)

    # def dist_matrix(self):
    #     """
    #     Experimental function, currently not in use
    #     """
    #     matrix = np.zeros((len(self), len(self)))
    #     for i in range(1, len(self)):
    #         for j in range(i, len(self)):
    #             matrix[i, j] = np.linalg.norm(self.positions[i] - self.positions[j])
    #     return matrix

    def get_bonds(self) -> list:
        """
        Returns a list than contains all bonds
        """
        bonds = []
        for edge in self.connectivity.edges:
            bonds.append((edge[0], edge[1], self.connectivity.edges[edge]["weight"]))
        return bonds

    def __len__(self):
        if isinstance(self.numbers, np.ndarray):
            return self.numbers.size
        elif isinstance(self.numbers, int):
            return 1

    def __eq__(self, other):
        """
        Check for identity of two frame objects.
        """
        try:
            return (
                len(self) == len(other)
                and (self.positions == other.positions).all()
                and (self.numbers == other.numbers).all()
                and (self.cell == other.cell).all()
                and (self.pbc == other.pbc).all()
                and (self.connectivity == other.connectivity)
            )
        except AttributeError:
            return NotImplemented

    def to_dict(self) -> dict:
        """
        Creates a dictionary than contains all the relevant information of the Frame object
        """
        frame_dict = {}

        for field in dataclasses.fields(self):
            frame_dict[field.name] = getattr(self, field.name)

        for key, value in frame_dict.items():
            if isinstance(value, np.ndarray):
                frame_dict[key] = value.tolist()

        if frame_dict["colors"] is None:
            frame_dict["colors"] = [
                rgb2hex(jmol_colors[number]) for number in frame_dict["numbers"]
            ]

        if frame_dict["radii"] is None:
            frame_dict["radii"] = [
                get_radius(number) for number in frame_dict["numbers"]
            ]

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
    def from_dict(cls, data):
        """
        Creates an instance of the Frame class from a dictionary
        """
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
            "vector_field" in data
        ):  # currently there is the vector_field part in js missing. So this is useless at the moment
            frame.vector_field = data["vector_field"]

        if "connectivity" in data:
            frame.connectivity = nx.Graph()
            for edge in data["connectivity"]:
                frame.connectivity.add_edge(edge[0], edge[1], weight=edge[2])

        return frame
