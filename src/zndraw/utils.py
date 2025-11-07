import functools

import ase
import numpy as np
from ase.data import covalent_radii
from ase.data.colors import jmol_colors


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

        # Store as numpy array with Unicode string dtype (not object dtype)
        # object dtype uses pickle encoding which is not supported in JavaScript
        atoms.set_array("colors", np.array(hex_colors, dtype="U7"))
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
