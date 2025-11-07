"""Test msgpack encoding/decoding for frame data.

This test verifies that the frame data endpoint returns properly encoded
msgpack data that can be decoded by the frontend (msgpack-numpy-js).

Expected format:
- Backend returns: dict[bytes, bytes] where values are msgpack-encoded with numpy support
- Frontend should decode to: dict[str, Any] where Any is numpy arrays or JSON types
"""

import msgpack
import msgpack_numpy as m
import numpy as np
import requests
import asebytes

from zndraw import ZnDraw


def test_msgpack_encoding(server, s22):
    """Test that a single frame is properly encoded as dict[bytes, bytes]."""
    vis = ZnDraw(url=server, room="msgpack-test", user="user1")
    vis.extend(s22)

    response = requests.get(
        f"{server}/api/rooms/msgpack-test/frames",
        params={"indices": "0,1", "keys": "arrays.positions,arrays.numbers"},
    )

    assert response.status_code == 200
    assert response.headers["Content-Type"] == "application/msgpack"
    raw_data = response.content
    unpacked_data = msgpack.unpackb(raw_data, raw=True)
    assert isinstance(unpacked_data, list)
    assert len(unpacked_data) == 2
    
    ## Key 0 ##
    assert unpacked_data[0].keys() == {b'arrays.positions', b'arrays.numbers'}
    positions = msgpack.unpackb(unpacked_data[0][b'arrays.positions'], object_hook=m.decode)
    assert isinstance(positions, np.ndarray)
    assert positions.shape == (8, 3)
    # compare ASE Atoms objects
    atoms = asebytes.decode(unpacked_data[0])
    assert atoms == s22[0]

    ## Key 1 ##
    assert unpacked_data[1].keys() == {b'arrays.positions', b'arrays.numbers'}
    positions = msgpack.unpackb(unpacked_data[1][b'arrays.positions'], object_hook=m.decode)
    assert isinstance(positions, np.ndarray)
    assert positions.shape == (6, 3)
    # compare ASE Atoms objects
    atoms = asebytes.decode(unpacked_data[1])
    assert atoms == s22[1]
