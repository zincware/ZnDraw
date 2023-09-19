import ase.io
import networkx as nx
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data.colors import jmol_colors

from zndraw.bonds import ASEComputeBonds


def _rgb2hex(value):
    r, g, b = np.array(value * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)


def _get_radius(value):
    return (0.25 * (2 - np.exp(-0.2 * value)),)


def atoms_to_json(atoms: ase.Atoms) -> dict:
    ase_bond_calculator = ASEComputeBonds()
    if not hasattr(atoms, "connectivity"):
        # TODO this should be removed over time.
        atoms.connectivity = ase_bond_calculator.build_graph(atoms)

    atoms_dict = atoms.todict()
    for key in list(atoms_dict):
        # includes ['numbers', 'positions', 'cell', 'pbc']
        if isinstance(atoms_dict[key], np.ndarray):
            atoms_dict[key] = atoms_dict[key].tolist()
        elif isinstance(atoms_dict[key], np.generic):
            atoms_dict[key] = atoms_dict[key].item()

    # remove info if available # currently not used
    atoms_dict.pop("info", None)

    atoms_dict["colors"] = [
        _rgb2hex(jmol_colors[number]) for number in atoms_dict["numbers"]
    ]
    atoms_dict["radii"] = [_get_radius(number) for number in atoms_dict["numbers"]]

    try:
        calc_data = {}
        for key in atoms.calc.results:
            value = atoms.calc.results[key]
            if isinstance(value, np.ndarray):
                value = value.tolist()
            calc_data[key] = value

        atoms_dict["calc"] = calc_data
    except (RuntimeError, AttributeError):
        pass

    try:
        atoms_dict["connectivity"] = ase_bond_calculator.get_bonds(atoms)
    except AttributeError:
        atoms_dict["connectivity"] = []

    return atoms_dict


def atoms_from_json(data: dict) -> ase.Atoms:
    atoms = ase.Atoms(
        numbers=data["numbers"],
        cell=data["cell"],
        pbc=data["pbc"],
        positions=data["positions"],
    )

    if "calc" in data:
        atoms.calc = SinglePointCalculator(atoms)
        atoms.calc.results = {
            key: np.array(val) if isinstance(val, list) else val
            for key, val in data["calc"].items()
        }

    if "connectivity" in data:
        atoms.connectivity = nx.Graph()
        for edge in data["connectivity"]:
            atoms.connectivity.add_edge(edge[0], edge[1], weight=edge[2])

    return atoms
