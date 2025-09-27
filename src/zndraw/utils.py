import functools
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
import ase
import numpy as np


@functools.lru_cache(maxsize=128)
def get_scaled_radii() -> np.ndarray:
    """Scale down the covalent radii to visualize bonds better."""
    radii = covalent_radii
    # shift the values such that they are in [0.3, 1.3]
    radii = radii - np.min(radii)
    radii = radii / np.max(radii)
    radii = radii + 0.3
    return radii

def update_colors_and_radii(atoms: ase.Atoms) -> None:
    """Update the colors and radii of the atoms in-place."""
    if "colors" not in atoms.arrays:
        colors = np.array([jmol_colors[atom.number] for atom in atoms], dtype=np.float32)
        atoms.set_array("colors", colors)
    if "radii" not in atoms.arrays:
         radii = np.array([covalent_radii[atom.number] for atom in atoms], dtype=np.float32)
         atoms.set_array("radii", radii)
