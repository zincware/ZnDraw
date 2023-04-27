"""Tools used for selecting objects."""
import ase
import networkx as nx

from zndraw import io


def select_identical_species(atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
    """Select all atoms with the same species as the selected atoms."""
    selected_ids = set(selected_ids)
    for idx in tuple(selected_ids):
        selected_symbol = atoms[idx].symbol
        selected_ids.update(
            idx for idx, atom in enumerate(atoms) if atom.symbol == selected_symbol
        )
    return list(selected_ids)


def select_connected(atoms: ase.Atoms, selected_ids: list[int]) -> list[int]:
    """Select all atoms connected to the selected atoms."""
    graph = io.get_graph(atoms)

    total_ids = []

    for node_id in selected_ids:
        total_ids += list(nx.node_connected_component(graph, node_id))

    return list(set(total_ids))
