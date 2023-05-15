import functools

import ase
import networkx as nx
import numpy as np
from ase.data.colors import jmol_colors
from ase.neighborlist import build_neighbor_list

from zndraw import shared


def _rgb2hex(data):
    r, g, b = np.array(data * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)


def serialize_atoms(start: int, stop: int):
    result = {"position": [], "force": [], "box": []}
    try:
        for step in range(start, stop):
            atoms = shared.config.get_atoms(step=int(step))
            result["position"].append(atoms.get_positions().tolist())
            result["box"].append(atoms.get_cell().diagonal().tolist())
            # TODO MAKE THIS OPTIONAL!!, also energy, etc.
            # try:
            #     result["force"].append(atoms.get_forces().tolist())
            # except:
            #     result["force"].append(np.zeros_like(atoms.get_positions()).tolist())
        return result
    except KeyError:
        return result


@functools.lru_cache(maxsize=16)
def get_bonds(step: int):
    atoms = shared.config.get_atoms(step=step)
    atoms.pbc = False
    nl = build_neighbor_list(atoms, self_interaction=False)
    cm = nl.get_connectivity_matrix(sparse=False)
    G = nx.from_numpy_array(cm)
    return list(G.edges)


def serialize_frame(step):
    atoms = shared.config.get_atoms(step=int(step))
    return {
        "particles": [
            {
                "id": idx,
                "x": atom.position[0],
                "y": atom.position[1],
                "z": atom.position[2],
                "color": _rgb2hex(jmol_colors[atom.number]),
                "radius": 0.25 * (2 - np.exp(-0.2 * atom.number)),
                # "species": atom.species,
            }
            for idx, atom in enumerate(atoms)
        ],
        "bonds": get_bonds(int(step)),
    }
