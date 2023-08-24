import dataclasses
import functools
import pathlib
import collections.abc

import ase.io
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
class DataHandler(collections.abc.MutableSequence):
    client: Client

    ase_bond_calculator: ASEComputeBonds = dataclasses.field(
        default_factory=ASEComputeBonds
    )

    def __delitem__(self, index):
        pass

    def __getitem__(self, index) -> list[ase.Atoms]:
        df = self.get_dataset()
        value = list(df.loc[index]["atoms"].compute().to_dict().values())
        return value if len(value) > 1 else value[0]

    def __len__(self):
        return len(self.get_dataset())

    def __setitem__(self, index, value):
        df = self.get_dataset().compute()

        if isinstance(index, slice):
            indices = list(range(index.stop))[index]
            for idx, atoms in zip(indices, value):
                df.loc[idx, "atoms"] = atoms
        elif isinstance(index, int):
            df.loc[index, "atoms"] = value
        else:
            raise TypeError(f"Index must be int or slice not {type(index)}")
        
        self.client.unpublish_dataset("atoms")
        self.client.restart()
        df = dd.from_pandas(df, npartitions=10)
        self.client.persist(df)
        self.client.publish_dataset(atoms=df)

    def insert(self, index, value):
        pass

    def create_dataset(self, filename):
        # TODO this should happen on a worker

        if pathlib.Path(filename).suffix == ".h5":
            # Read file using znh5md and convert to list[ase.Atoms]
            atoms_list = znh5md.ASEH5MD(filename).get_atoms_list()
        else:
            atoms_list = list(ase.io.iread(filename))

        for atoms in atoms_list:
            atoms.connectivity = self.ase_bond_calculator.build_graph(atoms)

        df = dd.DataFrame.from_dict({"atoms": atoms_list}, npartitions=10)
        # move dataframe to client and keep it there
        self.client.persist(df)
        # make dataset available with the id "atoms"
        self.client.publish_dataset(atoms=df)  # = **{"atoms": df}

    def get_dataset(self):
        return self.client.get_dataset("atoms")
    
    def index(self):
        return self.get_dataset().index.compute()

    def get_atoms_json(self, item) -> dict[str, ase.Atoms]:
        data = {}
        
        for idx, atoms in zip(self.index()[item], self[item]):
            atoms_dict = atoms.todict()
            for key in list(atoms_dict):
                # includes ['numbers', 'positions', 'cell', 'pbc']
                if isinstance(atoms_dict[key], np.ndarray):
                    atoms_dict[key] = atoms_dict[key].tolist()
                elif isinstance(atoms_dict[key], np.generic):
                    atoms_dict[key] = atoms_dict[key].item()

            # remove info if available # currently not used
            atoms_dict.pop("info", None)

            atoms_dict["colors"] = [
                _rgb2hex(jmol_colors[number]) for number in atoms_dict["numbers"]
            ]
            atoms_dict["radii"] = [
                _get_radius(number) for number in atoms_dict["numbers"]
            ]
            try:
                atoms_dict["connectivity"] = self.ase_bond_calculator.get_bonds(atoms)
            except AttributeError:
                atoms_dict["connectivity"] = []
            data[idx] = atoms_dict

        return data
