import random
import typing as t
from typing import Any

import ase
import networkx as nx
from pydantic import BaseModel, Field


class SelectionBase(BaseModel):
    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        raise NotImplementedError()


class NoneSelection(SelectionBase):
    method: t.Literal["NoneSelection"] = Field("NoneSelection")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return []


class All(SelectionBase):
    """Select all atoms."""

    method: t.Literal["All"] = Field("All")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return list(range(len(atoms)))


class Invert(SelectionBase):
    method: t.Literal["Invert"] = Field("Invert")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return list(set(range(len(atoms))) - set(selected_ids))


class Range(SelectionBase):
    method: t.Literal["Range"] = Field("Range")

    start: int = Field(..., description="Start index")
    end: int = Field(..., description="End index")
    step: int = Field(1, description="Step size")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return list(range(self.start, self.end, self.step))


class Random(SelectionBase):
    method: t.Literal["Random"] = Field("Random")

    count: int = Field(..., description="Number of atoms to select")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        return random.sample(range(len(atoms)), self.count)


class IdenticalSpecies(SelectionBase):
    method: t.Literal["IdenticalSpecies"] = Field("IdenticalSpecies")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        selected_ids = set(selected_ids)
        for idx in tuple(selected_ids):
            selected_symbol = atoms[idx].symbol
            selected_ids.update(
                idx for idx, atom in enumerate(atoms) if atom.symbol == selected_symbol
            )
        return list(selected_ids)


class ConnectedParticles(SelectionBase):
    method: t.Literal["ConnectedParticles"] = Field("ConnectedParticles")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        total_ids = []
        try:
            graph = atoms.connectivity
        except AttributeError:
            return selected_ids

        for node_id in selected_ids:
            total_ids += list(nx.node_connected_component(graph, node_id))

        return list(set(total_ids))


class Neighbour(SelectionBase):
    """Select the nth order neighbours of the selected atoms."""

    method: t.Literal["Neighbour"] = Field("Neighbour")

    order: int = Field(1, description="Order of neighbour")

    def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
        total_ids = []
        try:
            graph = atoms.connectivity
        except AttributeError:
            return selected_ids

        for node_id in selected_ids:
            total_ids += list(
                nx.single_source_shortest_path_length(graph, node_id, self.order).keys()
            )

        return list(set(total_ids))


def get_selection_class(methods):
    class Selection(SelectionBase):
        method: methods = Field(
            ..., description="Selection method", discriminator="method"
        )

        def get_ids(self, atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
            return self.method.get_ids(atoms, selected_ids)

        @classmethod
        def model_json_schema(cls, *args, **kwargs) -> dict[str, Any]:
            schema = super().model_json_schema(*args, **kwargs)
            for prop in [x.__name__ for x in t.get_args(methods)]:
                schema["$defs"][prop]["properties"]["method"]["options"] = {
                    "hidden": True
                }
                schema["$defs"][prop]["properties"]["method"]["type"] = "string"

            return schema

    return Selection
