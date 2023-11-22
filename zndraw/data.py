import ase.io
import networkx as nx
import numpy as np
import typing as t
from ase.data.colors import jmol_colors

def _rgb2hex(value):
    r, g, b = np.array(value * 255, dtype=int)
    return "#%02x%02x%02x" % (r, g, b)

def _get_radius(value):
    return (0.25 * (2 - np.exp(-0.2 * value)),)
