import functools

import ase
import networkx as nx
import numpy as np
from ase.data.colors import jmol_colors
from ase.neighborlist import build_neighbor_list, natural_cutoffs
from pydantic import BaseModel, Field

from zndraw import shared


def _rgb2hex(data):
    r, g, b = np.array(data * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)


class ASEComputeBonds(BaseModel):
    single_bond_multiplier: float = Field(1.0, le=10, ge=0)
    double_bond_multiplier: float = Field(0.85, le=10, ge=0)
    triple_bond_multiplier: float = Field(0.7, le=10, ge=0)

    def get_frame(self, step: int):
        atoms = shared.config.get_atoms(step=int(step))
        atoms.info["graph_representation"] = self.build_graph(
            atoms
        )  ## TODO: Can probably add a try except here to avoid rebuilding the graph every time
        self.update_graph_using_modifications(atoms)
        return {
            "particles": [
                {
                    "id": idx,
                    "x": atom.position[0],
                    "y": atom.position[1],
                    "z": atom.position[2],
                    "color": _rgb2hex(jmol_colors[atom.number]),
                    "radius": 0.25 * (2 - np.exp(-0.2 * atom.number)),
                    # "species": atom.species,
                }
                for idx, atom in enumerate(atoms)
            ],
            "bonds": self.get_bonds(atoms),
        }

    def build_graph(self, atoms: ase.Atoms):
        cutoffs = [
            self.single_bond_multiplier,
            self.double_bond_multiplier,
            self.triple_bond_multiplier,
        ]
        connectivity_matrix = np.zeros((len(atoms), len(atoms)), dtype=int)
        for cutoff in cutoffs:
            nl = build_neighbor_list(
                atoms, cutoffs=natural_cutoffs(atoms, cutoff), self_interaction=False
            )
            connectivity_matrix += nl.get_connectivity_matrix(sparse=False)
        G = nx.from_numpy_array(connectivity_matrix)
        return G

    def update_graph_using_modifications(self, atoms: ase.Atoms):
        modifications = atoms.info.get("modifications", {})
        graph = atoms.info["graph_representation"]
        for key in modifications:
            atom_1, atom_2 = key
            weight = modifications[key]
            graph.add_edge(atom_1, atom_2, weight=weight)

    def get_bonds(self, atoms: ase.Atoms):
        graph = atoms.info["graph_representation"]
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
