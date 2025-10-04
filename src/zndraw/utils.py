import functools

import ase
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.data import covalent_radii
from ase.data.colors import jmol_colors


def _convert_numpy_scalars(obj):
    """Convert numpy scalar types to native Python types.
    
    This function recursively converts numpy scalars (int64, float64, bool_, etc.)
    to their Python equivalents. Numpy arrays are left unchanged.
    """
    if isinstance(obj, np.generic):
        # Use .item() to convert any numpy scalar to native Python type
        return obj.item()
    elif isinstance(obj, dict):
        return {key: _convert_numpy_scalars(value) for key, value in obj.items()}
    elif isinstance(obj, (list, tuple)):
        return type(obj)(_convert_numpy_scalars(item) for item in obj)
    else:
        return obj


def atoms_to_dict(atoms: ase.Atoms) -> dict:
    if not atoms.calc:
        result = atoms.todict()
    else:
        result = atoms.todict()
        result["<SinglePointCalculator>"] = atoms.calc.results
    # return result
    # Convert any numpy scalars to Python native types to avoid JSON serialization errors
    return _convert_numpy_scalars(result)


def atoms_from_dict(d: dict) -> ase.Atoms:
    if "<SinglePointCalculator>" not in d:
        return ase.Atoms.fromdict(d)
    calc_results = d.pop("<SinglePointCalculator>")
    atoms = ase.Atoms.fromdict(d)
    atoms.calc = SinglePointCalculator(atoms)
    atoms.calc.results = calc_results
    return atoms


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
        colors = np.array(
            [
                jmol_colors[atom.number]
                if atom.number < len(jmol_colors)
                else [0.0, 0.0, 0.0]
                for atom in atoms
            ],
            dtype=np.float32,
        )
        atoms.set_array("colors", colors)
    if "radii" not in atoms.arrays:
        radii = np.array(
            [
                covalent_radii[atom.number]
                if atom.number < len(covalent_radii)
                else 0.77
                for atom in atoms
            ],
            dtype=np.float32,
        )
        atoms.set_array("radii", radii)
