import random
import typing as t

import networkx as nx
import numpy as np
from pydantic import Field

from zndraw.base import Extension
from zndraw.zndraw import ZnDraw

try:
    from zndraw.select import mda  # noqa: F401
except ImportError:
    # mdanalysis is not installed
    pass


class NoneSelection(Extension):
    def run(self, vis) -> None:
        vis.selection = []


class All(Extension):
    """Select all atoms."""

    def run(self, vis) -> None:
        atoms = vis[vis.step]
        vis.selection = list(range(len(atoms)))


class Invert(Extension):
    def run(self, vis) -> None:
        atoms = vis[vis.step]
        selected_ids = vis.selection
        vis.selection = list(set(range(len(atoms))) - set(selected_ids))


class Range(Extension):
    start: int = Field(0, description="Start index")
    end: int = Field(5, description="End index")
    step: int = Field(1, description="Step size")

    def run(self, vis) -> None:
        vis.selection = list(range(self.start, self.end, self.step))


class Random(Extension):
    count: int = Field(..., description="Number of atoms to select")

    def run(self, vis) -> None:
        atoms = vis[vis.step]
        vis.selection = random.sample(range(len(atoms)), self.count)


class IdenticalSpecies(Extension):
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


class ConnectedParticles(Extension):
    def run(self, vis) -> None:
        atoms = vis.atoms
        selected_ids = vis.selection
        total_ids = []
        try:
            edges = atoms.connectivity
            graph = nx.Graph()
            for edge in edges:
                node_a, node_b, weight = edge
                graph.add_edge(node_a, node_b, weight=weight)
        except AttributeError:
            return selected_ids

        for node_id in selected_ids:
            total_ids += list(nx.node_connected_component(graph, node_id))
            total_ids = np.array(total_ids)

        vis.selection = [x.item() for x in set(total_ids)]


class Neighbour(Extension):
    """Select the nth order neighbours of the selected atoms."""

    order: int = Field(1, description="Order of neighbour")

    def run(self, vis) -> None:
        total_ids = []
        atoms = vis[vis.step]
        selected_ids = vis.selection
        try:
            graph = atoms.connectivity
        except AttributeError:
            return selected_ids

        for node_id in selected_ids:
            total_ids += list(
                nx.single_source_shortest_path_length(graph, node_id, self.order).keys()
            )

        vis.selection = list(set(total_ids))


class UpdateSelection(Extension):
    """Reload Selection."""

    def run(self, vis: ZnDraw) -> None:
        vis.selection = vis.selection


selections: dict[str, t.Type[Extension]] = {
    ConnectedParticles.__name__: ConnectedParticles,
    NoneSelection.__name__: NoneSelection,
    All.__name__: All,
    Invert.__name__: Invert,
    Range.__name__: Range,
    Random.__name__: Random,
    IdenticalSpecies.__name__: IdenticalSpecies,
    Neighbour.__name__: Neighbour,
}
