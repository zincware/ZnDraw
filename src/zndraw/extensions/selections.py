import random
import typing as t

import networkx as nx
import numpy as np
from pydantic import Field

from zndraw.extensions.abc import Extension, Category


class Selection(Extension):
    """The base class for all selection extensions."""

    category: t.ClassVar[Category] = Category.SELECTION


class NoneSelection(Selection):
    def run(self, vis) -> None:
        vis.selection = []


class All(Selection):
    """Select all atoms."""

    def run(self, vis) -> None:
        atoms = vis[vis.step]
        vis.selection = list(range(len(atoms)))


class Invert(Selection):
    def run(self, vis) -> None:
        atoms = vis[vis.step]
        selected_ids = vis.selection
        vis.selection = list(set(range(len(atoms))) - set(selected_ids))


class Range(Selection):
    start: int = Field(0, description="Start index")
    end: int = Field(5, description="End index")
    step: int = Field(1, description="Step size")

    def run(self, vis) -> None:
        vis.selection = list(range(self.start, self.end, self.step))


class Random(Selection):
    count: int = Field(..., description="Number of atoms to select")
    seed: int = Field(42, description="Random seed for reproducibility")

    def run(self, vis) -> None:
        atoms = vis[vis.step]
        random.seed(self.seed)
        vis.selection = random.sample(range(len(atoms)), self.count)


class IdenticalSpecies(Selection):
    def run(self, vis) -> None:
        atoms = vis[vis.step]
        selected_ids = vis.selection
        selected_ids = set(selected_ids)
        for idx in tuple(selected_ids):
            selected_symbol = atoms[idx].symbol
            selected_ids.update(
                idx for idx, atom in enumerate(atoms) if atom.symbol == selected_symbol
            )
        vis.selection = list(selected_ids)


class ConnectedParticles(Selection):
    def run(self, vis) -> None:
        atoms = vis.atoms
        selected_ids = vis.selection
        total_ids = []
        try:
            edges = atoms.info.get("connectivity", [])
            graph = nx.Graph()
            for edge in edges:
                node_a, node_b, weight = edge
                graph.add_edge(node_a, node_b, weight=weight)
        except (AttributeError, KeyError):
            return selected_ids

        for node_id in selected_ids:
            total_ids += list(nx.node_connected_component(graph, node_id))

        vis.selection = [int(x) for x in set(total_ids)]


class Neighbour(Selection):
    """Select the nth order neighbours of the selected atoms."""

    order: int = Field(1, description="Order of neighbour")

    def run(self, vis) -> None:
        total_ids = []
        atoms = vis[vis.step]
        selected_ids = vis.selection
        try:
            connectivity = atoms.info.get("connectivity", [])
            graph = nx.Graph()
            for u, v, weight in connectivity:
                graph.add_edge(u, v, weight=weight)
        except (AttributeError, KeyError):
            return selected_ids

        for node_id in selected_ids:
            total_ids += list(
                nx.single_source_shortest_path_length(graph, node_id, self.order).keys()
            )

        vis.selection = [int(x) for x in set(total_ids)]


class UpdateSelection(Selection):
    """Reload Selection."""

    def run(self, vis) -> None:
        vis.selection = vis.selection


selections: dict[str, t.Type[Selection]] = {
    ConnectedParticles.__name__: ConnectedParticles,
    NoneSelection.__name__: NoneSelection,
    All.__name__: All,
    Invert.__name__: Invert,
    Range.__name__: Range,
    Random.__name__: Random,
    IdenticalSpecies.__name__: IdenticalSpecies,
    Neighbour.__name__: Neighbour,
    UpdateSelection.__name__: UpdateSelection,
}
