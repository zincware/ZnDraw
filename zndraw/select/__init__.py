import random
import typing as t
from typing import Any

import networkx as nx
import numpy as np
from pydantic import BaseModel, Field

try:
    from zndraw.select import mda  # noqa: F401
except ImportError:
    # mdanalysis is not installed
    pass

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class SelectionBase(BaseModel):
    def run(self, vis: "ZnDraw") -> None:
        raise NotImplementedError()


class NoneSelection(SelectionBase):
    discriminator: t.Literal["NoneSelection"] = Field("NoneSelection")

    def run(self, vis) -> None:
        vis.selection = []


class All(SelectionBase):
    """Select all atoms."""

    discriminator: t.Literal["All"] = Field("All")

    def run(self, vis) -> None:
        atoms = vis[vis.step]
        vis.selection = list(range(len(atoms)))


class Invert(SelectionBase):
    discriminator: t.Literal["Invert"] = Field("Invert")

    def run(self, vis) -> None:
        atoms = vis[vis.step]
        selected_ids = vis.selection
        vis.selection = list(set(range(len(atoms))) - set(selected_ids))


class Range(SelectionBase):
    discriminator: t.Literal["Range"] = Field("Range")

    start: int = Field(..., description="Start index")
    end: int = Field(..., description="End index")
    step: int = Field(1, description="Step size")

    def run(self, vis) -> None:
        vis.selection = list(range(self.start, self.end, self.step))


class Random(SelectionBase):
    discriminator: t.Literal["Random"] = Field("Random")

    count: int = Field(..., description="Number of atoms to select")

    def run(self, vis) -> None:
        atoms = vis[vis.step]
        vis.selection = random.sample(range(len(atoms)), self.count)


class IdenticalSpecies(SelectionBase):
    discriminator: t.Literal["IdenticalSpecies"] = Field("IdenticalSpecies")

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


class ConnectedParticles(SelectionBase):
    discriminator: t.Literal["ConnectedParticles"] = Field("ConnectedParticles")

    def run(self, vis) -> None:
        atoms = vis[vis.step]
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
            total_ids = np.array(total_ids, dtype=int)

        vis.selection = [x.item() for x in set(total_ids)]


class Neighbour(SelectionBase):
    """Select the nth order neighbours of the selected atoms."""

    discriminator: t.Literal["Neighbour"] = Field("Neighbour")

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


def get_selection_class(methods):
    class Selection(SelectionBase):
        method: methods = Field(
            ..., description="Selection method", discriminator="discriminator"
        )

        def run(self, vis: "ZnDraw") -> None:
            self.method.run(vis)

        @classmethod
        def model_json_schema(cls, *args, **kwargs) -> dict[str, Any]:
            schema = super().model_json_schema(*args, **kwargs)
            for prop in [x.__name__ for x in t.get_args(methods)]:
                schema["$defs"][prop]["properties"]["discriminator"]["options"] = {
                    "hidden": True
                }
                schema["$defs"][prop]["properties"]["discriminator"]["type"] = "string"

            return schema

    return Selection
