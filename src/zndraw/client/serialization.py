"""Frame serialization helpers for the ZnDraw client."""

from __future__ import annotations

import base64
from typing import Any

import ase
from asebytes import decode, encode

from zndraw.enrichment import add_colors, add_radii

# =============================================================================
# Frame Serialization Helpers
# =============================================================================


def atoms_to_json_dict(
    atoms: ase.Atoms, connectivity_threshold: int = 1000
) -> dict[str, Any]:
    """Convert ase.Atoms to a JSON-compatible dictionary.

    Ensures colors, radii, and connectivity are present, then uses asebytes.encode()
    to get dict[bytes, bytes], then converts:
    - Byte keys to base64-encoded strings with "b64:" prefix
    - Byte values to base64-encoded strings

    Parameters
    ----------
    atoms
        The ASE Atoms object to convert.
    connectivity_threshold
        Maximum number of atoms for automatic connectivity calculation.
        Connectivity is only computed if the number of atoms is below this
        threshold and connectivity is not already present. Default: 1000.
    """
    add_colors(atoms)
    add_radii(atoms)

    if len(atoms) < connectivity_threshold and "connectivity" not in atoms.info:
        from zndraw.connectivity import add_connectivity

        add_connectivity(atoms)

    encoded = encode(atoms)
    result: dict[str, Any] = {}
    for key, value in encoded.items():
        # Encode key as base64 string
        key_str = "b64:" + base64.b64encode(key).decode("ascii")
        # Encode value as base64 string
        value_str = base64.b64encode(value).decode("ascii")
        result[key_str] = value_str
    return result


def json_dict_to_atoms(data: dict[str, Any]) -> ase.Atoms:
    """Convert a JSON dictionary back to ase.Atoms.

    Reverses atoms_to_json_dict(): decodes base64 keys/values back to bytes,
    then uses asebytes.decode().
    """
    encoded: dict[bytes, bytes] = {}
    for key_str, value_str in data.items():
        # Decode key from base64 (strip "b64:" prefix if present)
        if key_str.startswith("b64:"):
            key = base64.b64decode(key_str[4:])
        else:
            # Legacy format or string key - convert to bytes
            key = key_str.encode("utf-8") if isinstance(key_str, str) else key_str
        # Decode value from base64
        value = base64.b64decode(value_str) if isinstance(value_str, str) else value_str
        encoded[key] = value
    return decode(encoded)


def raw_frame_to_atoms(frame: dict[bytes, bytes]) -> ase.Atoms:
    """Convert a raw msgpack frame (dict[bytes, bytes]) to ase.Atoms.

    This is used when frames come directly from msgpack (server response).
    """
    return decode(frame)


# Target byte size per upload chunk. Frames are accumulated until adding
# the next frame would exceed this limit, then the chunk is flushed.
_TARGET_CHUNK_BYTES = 2_000_000  # 2 MB

# Hard cap on frames per chunk -- must not exceed the server's
# FrameCreateRequest.max_length (1000).
_MAX_CHUNK_FRAMES = 1000


def _estimate_frame_size(frame: dict[str, Any]) -> int:
    """Estimate serialized size of a JSON-encoded frame dict.

    Uses the sum of base64 value lengths as a cheap proxy -- no
    extra serialization needed since we already have the strings.
    Intentionally ignores key lengths (~20-30% undercount) since this
    is only used for chunking heuristics, not exact measurement.
    """
    return sum(len(v) for v in frame.values() if isinstance(v, str))
