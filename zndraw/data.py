import ase.io
import networkx as nx
import numpy as np
import typing as t
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data.colors import jmol_colors

from zndraw.bonds import ASEComputeBonds
from zndraw.frame import Frame

def _rgb2hex(value):
    r, g, b = np.array(value * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)

def _get_radius(value):
    return (0.25 * (2 - np.exp(-0.2 * value)),)

def atoms_to_json(frame: t.Union[ase.Atoms, Frame]):
    if isinstance(frame, ase.Atoms):
        frame = Frame.from_atoms(frame)

    ase_bond_calculator = ASEComputeBonds()

    if not hasattr(frame, "connectivity"):
        # TODO ase_bond_calculator for Frame Class
        #frame.connectivity = ase_bond_calculator.build_graph(frame)
        pass

    frame_dict = frame.__dict__

    for key, value in frame_dict.items():
        if isinstance(value, np.ndarray):
            frame_dict[key] = value.tolist()
        elif isinstance(value, np.generic):
            frame_dict[key] = value.item()

    if frame_dict["colors"] is None:
        frame_dict["colors"] = [
            _rgb2hex(jmol_colors[number]) for number in frame_dict["numbers"]
        ]

    frame_dict["radii"] = [_get_radius(number) for number in frame_dict["numbers"]]


    try:
        calc_data = {}
        for key in frame.calc.results:
            value = frame.calc.results[key]
            if isinstance(value, np.ndarray):
                value = value.tolist()
            calc_data[key] = value

        frame_dict["calc"] = calc_data
    except (RuntimeError, AttributeError):
        frame_dict["calc"] = None

    try:
        frame_dict["connectivity"] = []
        #frame_dict["connectivity"] = ase_bond_calculator.get_bonds(frame)
        # ase_bond_calculator für Frame anpassen TODO
    except AttributeError:
        frame_dict["connectivity"] = []
        
    return frame_dict

def atoms_from_json(data: dict) -> Frame:
    
    frame = Frame(positions=np.array(data["positions"]),
                  cell=np.array(data["cell"]),
                  numbers=np.array(data["numbers"]),
                  colors=np.array(data["colors"]),
                  radii=np.array(data["radii"]),
                  pbc=data["pbc"])

    #if data["calc"] is not None:
    #    frame.calc = SinglePointCalculator(atoms)
    #    frame.calc.results = {
    #        key: np.array(val) if isinstance(val, list) else val
    #        for key, val in data["calc"].items()
    #    }

    if "connectivity" in data:
        frame.connectivity = nx.Graph()
        for edge in data["connectivity"]:
            frame.connectivity.add_edge(edge[0], edge[1], weight=edge[2])

    return frame
