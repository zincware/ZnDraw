import pathlib

import ase.io
import networkx as nx
import numpy as np
import znh5md
from ase.data.colors import jmol_colors
from ase.neighborlist import build_neighbor_list


def read_file(filename: str) -> ase.Atoms:
    """Reads a file and returns a networkx graph."""
    if pathlib.Path(filename).suffix == ".h5":
        return znh5md.ASEH5MD(filename).get_atoms_list()
    return ase.io.read(filename)


def _rgb2hex(data):
    r, g, b = np.array(data * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)


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
            "color": _rgb2hex(jmol_colors[atoms[node].number]),
            "radius": 0.25 * (2 - np.exp(-0.2 * atoms[node].number)),
            "position": atoms[node].position.tolist(),
        }

    nx.set_node_attributes(G, node_data)

    return G
