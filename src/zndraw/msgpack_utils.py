"""Utilities for msgpack serialization of ZnDraw frame data."""

import msgpack
import msgpack_numpy as m


def dict_to_msgpack(data: dict) -> dict[bytes, bytes]:
    """Convert Python dict to msgpack format for storage/transmission.

    Parameters
    ----------
    data : dict
        Frame dictionary with string keys (e.g., "arrays.positions", "info.energy")

    Returns
    -------
    dict[bytes, bytes]
        Msgpack-serialized dictionary with bytes keys and msgpack-packed values
    """
    result = {}
    for key, value in data.items():
        key_bytes = key.encode()
        value_bytes = msgpack.packb(value, default=m.encode)
        result[key_bytes] = value_bytes
    return result


def msgpack_to_dict(data: dict[bytes, bytes]) -> dict:
    """Convert msgpack format to Python dict.

    Parameters
    ----------
    data : dict[bytes, bytes]
        Msgpack-serialized dictionary with bytes keys and msgpack-packed values

    Returns
    -------
    dict
        Frame dictionary with string keys and unpacked Python values
    """
    result = {}
    for key_bytes, value_bytes in data.items():
        key = key_bytes.decode()
        value = msgpack.unpackb(value_bytes, object_hook=m.decode)
        result[key] = value
    return result
