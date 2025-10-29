import functools

import ase
import ase.constraints
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


def _constraint_to_dict(constraint) -> dict | None:
    """Serialize an ASE constraint to a dictionary using ASE's todict method.

    Supports all constraints that inherit from ase.constraints.FixConstraint
    and implement the todict() method.

    Parameters
    ----------
    constraint
        An ASE constraint object.

    Returns
    -------
    dict or None
        Dictionary representation of the constraint using ASE's format,
        or None if the constraint type is not supported.
    """
    if isinstance(constraint, ase.constraints.FixConstraint):
        return constraint.todict()

    # Silently ignore unsupported constraint types
    return None


def _constraint_from_dict(constraint_dict: dict):
    """Deserialize an ASE constraint from a dictionary using ASE's dict2constraint.

    Supports all constraints that can be reconstructed via dict2constraint.

    Parameters
    ----------
    constraint_dict : dict
        Dictionary representation of a constraint in ASE's format.

    Returns
    -------
    constraint
        An ASE constraint object, or None if unsupported.
    """
    try:
        return ase.constraints.dict2constraint(constraint_dict)
    except (KeyError, ValueError, AttributeError):
        # Silently ignore unsupported constraint types
        return None


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

    # Handle constraints
    if atoms.constraints:
        constraints_list = []
        for constraint in atoms.constraints:
            constraint_dict = _constraint_to_dict(constraint)
            if constraint_dict is not None:
                constraints_list.append(constraint_dict)

        if constraints_list:
            prefixed_result["constraints"] = constraints_list

    # Convert any numpy scalars to Python native types to avoid JSON serialization errors
    return _convert_numpy_scalars(prefixed_result)


def atoms_from_dict(d: dict) -> ase.Atoms:
    # Reconstruct the original dict structure from prefixed keys
    reconstructed = {}
    info_dict = {}
    calc_dict = {}
    constraints_list = None

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
        elif key == "constraints":
            # Store constraints for later processing
            constraints_list = value
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

    # Restore constraints
    if constraints_list:
        constraint_objects = []
        for constraint_dict in constraints_list:
            constraint = _constraint_from_dict(constraint_dict)
            if constraint is not None:
                constraint_objects.append(constraint)

        if constraint_objects:
            atoms.set_constraint(constraint_objects)

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
    """Update the colors and radii of the atoms in-place.

    Colors are stored as hex strings instead of RGB arrays.
    """
    if "colors" not in atoms.arrays:
        # Get RGB colors from jmol_colors
        rgb_colors = np.array(
            [
                jmol_colors[atom.number]
                if atom.number < len(jmol_colors)
                else [0.0, 0.0, 0.0]
                for atom in atoms
            ],
            dtype=np.float32,
        )

        # Convert RGB (0-1 range) to hex strings
        hex_colors = []
        for rgb in rgb_colors:
            r, g, b = (np.clip(rgb, 0, 1) * 255).astype(int)
            hex_colors.append(f"#{r:02x}{g:02x}{b:02x}")

        # Store as numpy array of object dtype (required for ASE)
        atoms.set_array("colors", np.array(hex_colors, dtype=object))

    if "radii" not in atoms.arrays:
        radii = covalent_radii[atoms.numbers].astype(np.float32)
        atoms.set_array("radii", radii)


def generate_room_name(base_name: str, redis_client, max_length: int = 20) -> str:
    """Generate a unique room name from a base name.

    The room name is truncated to max_length characters. If a room with that name
    already exists, a random UUID suffix is appended to ensure uniqueness.

    Parameters
    ----------
    base_name : str
        The original filename or path to create a room name from.
    redis_client
        Redis client to check for existing room names.
    max_length : int
        Maximum length for the room name (default: 20).

    Returns
    -------
    str
        A unique room name.

    Examples
    --------
    >>> generate_room_name("very_long_filename.xyz", redis_client)
    'very_long_filename.'  # if unique
    >>> generate_room_name("structure.xyz", redis_client)  # if collision
    'structure.xyz_a3f2'  # random UUID suffix added
    """
    import uuid

    # Truncate to max length
    truncated = base_name[:max_length]

    # If no redis client provided, skip collision check
    if redis_client is None:
        return truncated

    # Check if this room name already exists
    room_exists = False
    for key in redis_client.scan_iter(match=f"room:{truncated}:*", count=1):
        room_exists = True
        break

    if not room_exists:
        return truncated

    # Room exists, add random UUID suffix for guaranteed uniqueness
    # Use first 4 characters of UUID hex for short, unique suffix
    uuid_suffix = uuid.uuid4().hex[:4]

    # Truncate further to make room for UUID suffix
    suffix_length = 5  # underscore + 4 hex chars
    if len(truncated) + suffix_length > max_length:
        truncated = truncated[: max_length - suffix_length]

    return f"{truncated}_{uuid_suffix}"
