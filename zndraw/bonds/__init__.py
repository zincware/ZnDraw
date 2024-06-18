import ase
import networkx as nx
import numpy as np
from ase.neighborlist import natural_cutoffs
from networkx.exception import NetworkXError
from pydantic import BaseModel, Field


class ASEComputeBonds(BaseModel):
    single_bond_multiplier: float = Field(1.2, le=2, ge=0)
    double_bond_multiplier: float = Field(0.9, le=1, ge=0)
    triple_bond_multiplier: float = Field(0.0, le=1, ge=0)

    def build_graph(self, atoms: ase.Atoms):
        cutoffs = [
            self.single_bond_multiplier,
            self.double_bond_multiplier,
            self.triple_bond_multiplier,
        ]
        atoms_copy = atoms.copy()
        connectivity_matrix = np.zeros((len(atoms_copy), len(atoms_copy)), dtype=int)
        atoms_copy.pbc = False
        distance_matrix = atoms_copy.get_all_distances(mic=False)
        np.fill_diagonal(distance_matrix, np.inf)
        for cutoff in cutoffs:
            cutoffs = np.array(natural_cutoffs(atoms_copy, mult=cutoff))
            cutoffs = cutoffs[:, None] + cutoffs[None, :]
            connectivity_matrix[distance_matrix <= cutoffs] += 1
        G = nx.from_numpy_array(connectivity_matrix)
        return G

    def update_graph_using_modifications(self, atoms: ase.Atoms):
        modifications = atoms.info.get("modifications", {})
        graph = atoms.connectivity
        for key in modifications:
            atom_1, atom_2 = key
            weight = modifications[key]
            if weight == 0:
                self.remove_edge(graph, atom_1, atom_2)
            else:
                graph.add_edge(atom_1, atom_2, weight=weight)

    @staticmethod
    def remove_edge(graph, atom_1, atom_2):
        try:
            graph.remove_edge(atom_1, atom_2)
        except NetworkXError:
            pass

    def get_bonds(self, atoms: ase.Atoms, graph: nx.Graph = None):
        if graph is None:
            graph = self.build_graph(atoms)
        bonds = []
        for edge in graph.edges:
            bonds.append((edge[0], edge[1], graph.edges[edge]["weight"]))
        return bonds

    def update_bond_order(self, atoms: ase.Atoms, particles: list[int], order: int):
        if len(particles) != 2:
            raise ValueError("Exactly two particles must be selected")
        modifications = atoms.info.get("modifications", {})
        sorted_particles = tuple(sorted(particles))
        modifications[sorted_particles] = order
        atoms.info["modifications"] = modifications
