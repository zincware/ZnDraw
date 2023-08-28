import random

import ase
import networkx as nx
from pydantic import BaseModel, Field


class SelectionBase(BaseModel):
    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        raise NotImplementedError()


class NoneSelection(SelectionBase):
    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return []


class All(SelectionBase):
    """Select all atoms."""

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return list(range(len(atoms)))


class Invert(SelectionBase):
    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return list(set(range(len(atoms))) - set(selected_ids))


class Range(SelectionBase):
    start: int = Field(..., description="Start index")
    end: int = Field(..., description="End index")
    step: int = Field(1, description="Step size")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return list(range(self.start, self.end, self.step))


class Random(SelectionBase):
    count: int = Field(..., description="Number of atoms to select")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return random.sample(range(len(atoms)), self.count)


class IdenticalSpecies(SelectionBase):
    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        selected_ids = set(selected_ids)
        for idx in tuple(selected_ids):
            selected_symbol = atoms[idx].symbol
            selected_ids.update(
                idx for idx, atom in enumerate(atoms) if atom.symbol == selected_symbol
            )
        return list(selected_ids)


class ConnectedParticles(SelectionBase):
    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        total_ids = []
        try:
            graph = atoms.connectivity
        except AttributeError:
            return selected_ids

        for node_id in selected_ids:
            total_ids += list(nx.node_connected_component(graph, node_id))

        return list(set(total_ids))
