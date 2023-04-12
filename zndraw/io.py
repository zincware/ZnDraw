import ase.io
import pandas as pd
import plotly.express as px

import networkx as nx
import numpy as np
from ase.neighborlist import build_neighbor_list


def read_file(filename: str) -> nx.Graph:
    atoms = ase.io.read(filename)
    return get_graph(atoms)


def get_graph(atoms: ase.Atoms) -> nx.Graph:
    nl = build_neighbor_list(atoms, self_interaction=False)
    cm = nl.get_connectivity_matrix(sparse=False)
    G = nx.from_numpy_array(cm)

    node_data = {}

    for node in G.nodes:
        node_data[node] = {
            "symbol": atoms[node].symbol,
            "number": atoms[node].number,
            "x": atoms[node].position[0],
            "y": atoms[node].position[1],
            "z": atoms[node].position[2],
        }

    nx.set_node_attributes(G, node_data)

    return G
