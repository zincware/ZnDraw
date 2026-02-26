"""Server-side frame enrichment: colors, radii, and connectivity.

Ensures frames have the derived arrays the frontend needs for rendering
(``arrays.colors``, ``arrays.radii``, ``info.connectivity``).
"""

import ase
import numpy as np
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from asebytes import decode, encode

from zndraw.connectivity import add_connectivity
from zndraw.storage.base import RawFrame

# Pre-computed hex color lookup table indexed by atomic number.
# Avoids per-atom Python loops — a single numpy fancy-index suffices.
_HEX_COLORS_LUT = np.array(
    [
        f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"
        for r, g, b in jmol_colors
    ],
    dtype="U7",
)

_KEY_COLORS = b"arrays.colors"
_KEY_RADII = b"arrays.radii"
_KEY_CONNECTIVITY = b"info.connectivity"


def add_colors(atoms: ase.Atoms) -> None:
    """Add colors array to *atoms* in-place if not present."""
    if "colors" not in atoms.arrays:
        safe_numbers = np.clip(atoms.numbers, 0, len(_HEX_COLORS_LUT) - 1)
        atoms.set_array("colors", _HEX_COLORS_LUT[safe_numbers])


def add_radii(atoms: ase.Atoms) -> None:
    """Add radii array to *atoms* in-place if not present."""
    if "radii" not in atoms.arrays:
        atoms.set_array("radii", covalent_radii[atoms.numbers].astype(np.float32))


def enrich_raw_frame(frame: RawFrame, connectivity_threshold: int = 100) -> RawFrame:
    """Ensure a raw frame has colors, radii, and connectivity.

    Returns the frame unchanged (no decode/encode round-trip) when all
    three keys are already present.

    Parameters
    ----------
    frame
        Raw frame data (``dict[bytes, bytes]``).
    connectivity_threshold
        Skip connectivity computation for frames with more atoms than this.
    """
    needs_colors = _KEY_COLORS not in frame
    needs_radii = _KEY_RADII not in frame
    needs_connectivity = _KEY_CONNECTIVITY not in frame

    if not (needs_colors or needs_radii or needs_connectivity):
        return frame

    atoms: ase.Atoms = decode(frame)

    if needs_colors:
        add_colors(atoms)

    if needs_radii:
        add_radii(atoms)

    if needs_connectivity and len(atoms) < connectivity_threshold:
        add_connectivity(atoms)

    return encode(atoms)
