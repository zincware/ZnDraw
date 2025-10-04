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
    result = atoms.todict()

    # Create a new dict with properly prefixed keys
    prefixed_result = {}

    # Get the list of keys that are actually in atoms.arrays
    array_keys = set(atoms.arrays.keys())

    # Process all keys from todict
    for key, value in result.items():
        if key in array_keys:
            # Keys in atoms.arrays get 'arrays.' prefix
            prefixed_result[f"arrays.{key}"] = value
        elif key == "info":
            # Expand info dict with 'info.' prefix
            # Convert scalar values to numpy arrays to preserve dtype information
            for info_key, info_value in value.items():
                if isinstance(info_value, (int, float, bool, np.number)):
                    prefixed_result[f"info.{info_key}"] = np.array(info_value)
                else:
                    prefixed_result[f"info.{info_key}"] = info_value
        else:
            # cell, pbc, and other non-array fields keep their original keys
            prefixed_result[key] = value

    # Handle calculator results with 'calc.' prefix
    # Convert scalar values to numpy arrays to preserve dtype information
    if atoms.calc:
        for calc_key, calc_value in atoms.calc.results.items():
            if isinstance(calc_value, (int, float, bool, np.number)):
                prefixed_result[f"calc.{calc_key}"] = np.array(calc_value)
            else:
                prefixed_result[f"calc.{calc_key}"] = calc_value

    # Convert any numpy scalars to Python native types to avoid JSON serialization errors
    return _convert_numpy_scalars(prefixed_result)


def atoms_from_dict(d: dict) -> ase.Atoms:
    # Reconstruct the original dict structure from prefixed keys
    reconstructed = {}
    info_dict = {}
    calc_dict = {}

    for key, value in d.items():
        if key.startswith("arrays."):
            # Remove 'arrays.' prefix for reconstruction
            array_key = key[7:]  # len("arrays.") = 7
            reconstructed[array_key] = value
        elif key.startswith("info."):
            # Collect info keys
            info_key = key[5:]  # len("info.") = 5
            info_dict[info_key] = value
        elif key.startswith("calc."):
            # Collect calc keys
            calc_key = key[5:]  # len("calc.") = 5
            calc_dict[calc_key] = value
        else:
            reconstructed[key] = value

    # Add info dict if it has any keys
    if info_dict:
        reconstructed["info"] = info_dict

    # Create atoms from reconstructed dict
    atoms = ase.Atoms.fromdict(reconstructed)

    # Add calculator if calc_dict has any keys
    if calc_dict:
        atoms.calc = SinglePointCalculator(atoms)
        atoms.calc.results = calc_dict

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
