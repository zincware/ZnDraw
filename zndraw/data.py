import dataclasses
import functools
import random

import ase
import dask.dataframe as dd
import networkx as nx
import numpy as np
import pandas as pd
import znh5md
from ase.data.colors import jmol_colors
from ase.neighborlist import build_neighbor_list, natural_cutoffs
from dask.distributed import Client, Lock, Variable
from networkx.exception import NetworkXError
from pydantic import BaseModel, Field


def _rgb2hex(value):
    r, g, b = np.array(value * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)


def _get_radius(value):
    return (0.25 * (2 - np.exp(-0.2 * value)),)


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
        connectivity_matrix = np.zeros((len(atoms), len(atoms)), dtype=int)
        atoms.pbc = False
        distance_matrix = atoms.get_all_distances(mic=False)
        np.fill_diagonal(distance_matrix, np.inf)
        for cutoff in cutoffs:
            cutoffs = np.array(natural_cutoffs(atoms, mult=cutoff))
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

    def get_bonds(self, atoms: ase.Atoms):
        graph = atoms.connectivity
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


@dataclasses.dataclass
class DataHandler:
    client: Client

    ase_bond_calculator: ASEComputeBonds = dataclasses.field(
        default_factory=ASEComputeBonds
    )

    def create_dataset(self, filename):
        # TODO this should happen on a worker
        atoms_list = znh5md.ASEH5MD(filename).get_atoms_list()

        for atoms in atoms_list:
            atoms.connectivity = self.ase_bond_calculator.build_graph(atoms)

        df = dd.DataFrame.from_dict({"atoms": atoms_list}, npartitions=10)
        self.client.persist(df)
        self.client.publish_dataset(atoms=df)

    def __len__(self):
        return len(self.get_dataset())

    def get_dataset(self):
        return self.client.get_dataset("atoms")

    def get_atoms(self, _slice) -> dict[str, ase.Atoms]:
        df = self.get_dataset()
        return df.loc[_slice]["atoms"].compute().to_dict()

    def get_atoms_json(self, _slice) -> dict[str, ase.Atoms]:
        data = {}
        for idx, atoms in self.get_atoms(_slice).items():
            atoms_dict = atoms.todict()
            for key in list(atoms_dict):
                # includes ['numbers', 'positions', 'cell', 'pbc']
                if isinstance(atoms_dict[key], np.ndarray):
                    atoms_dict[key] = atoms_dict[key].tolist()

            atoms_dict["colors"] = [
                _rgb2hex(jmol_colors[number]) for number in atoms_dict["numbers"]
            ]
            atoms_dict["radii"] = [
                _get_radius(number) for number in atoms_dict["numbers"]
            ]
            atoms_dict["connectivity"] = self.ase_bond_calculator.get_bonds(atoms)
            data[idx] = atoms_dict

        return data
