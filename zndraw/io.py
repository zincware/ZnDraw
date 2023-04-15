import ase.io

import networkx as nx
from ase.neighborlist import build_neighbor_list
import numpy as np


def read_file(filename: str) -> ase.Atoms:
    """Reads a file and returns a networkx graph."""
    return ase.io.read(filename)


colors = {"H": "white", "C": "grey", "N": "blue", "O": "red", "F": "green"}


def get_graph(atoms: ase.Atoms) -> nx.Graph:
    """Returns a networkx graph from an ase.Atoms object."""
    atoms.pbc = False
    nl = build_neighbor_list(atoms, self_interaction=False)
    cm = nl.get_connectivity_matrix(sparse=False)
    G = nx.from_numpy_array(cm)

    node_data = {}

    for node in G.nodes:
        node_data[node] = {
            "symbol": atoms[node].symbol,
            "number": atoms[node].number.item(),
            "color": colors.get(atoms[node].symbol, "skyblue"),
            "radius": 0.25 * (2 - np.exp(-0.2 * atoms[node].number)),
            "position": atoms[node].position.tolist(),
        }

    nx.set_node_attributes(G, node_data)

    return G
