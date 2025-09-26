import json

import numpy as np
import pytest
import zarr

from zndraw_communication.storage import (
    ZarrStorageSequence,
    create_zarr,
    decode_data,
    encode_data,
    extend_zarr,
    read_zarr,
)


def assert_equal(original, reconstructed):
    if isinstance(original, np.ndarray):
        assert isinstance(reconstructed, np.ndarray)
        assert original.shape == reconstructed.shape
        assert original.dtype == reconstructed.dtype
        np.testing.assert_array_equal(original, reconstructed)
    elif isinstance(original, dict):
        assert isinstance(reconstructed, dict)
        assert set(original.keys()) == set(reconstructed.keys())
        for k in original:
            assert_equal(original[k], reconstructed[k])
    else:
        assert original == reconstructed


@pytest.fixture
def sample_data():
    return {
        "numbers": np.array([6, 1, 1, 1, 1]),
        "positions": np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
            ]
        ),
        "tags": np.array([1, 2, 3, 4, 5]),
        "initial_charges": np.array(
            [0.68062226, 0.77952437, 0.98666274, 0.12193437, 0.94664167]
        ),
        "momenta": np.array(
            [
                [0.11615927, 0.03696042, 0.46672605],
                [0.15866194, 0.11610305, 0.83340719],
                [0.59826544, 0.67848578, 0.74957987],
                [0.14949151, 0.68352953, 0.92811033],
                [0.40530313, 0.61190158, 0.43241552],
            ]
        ),
        "name": np.array(["Carbon", "Hydrogen", "Hydrogen", "Hydrogen", "Hydrogen"]),
        "info": {
            "string": "Lorem Ipsum",
            "float": 3.14,
            "list": [1, 2, 3],
            "array": np.array([1.0, 2.0, 3.0]),
            "dict": {"a": 1, "b": 2},
            "d2": {"a": [1, 2], "b": [3, 4]},
            "REPEATED_KEY": "info",
        },
        "REPEATED_KEY": np.array([0, 1, 2, 3, 4]),
        "cell": np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]),
        "pbc": np.array([True, False, True]),
        "celldisp": np.array([0.1, 0.2, 0.3]),
        "<SinglePointCalculator>": {
            "energy": 1.0,
            "forces": np.array(
                [
                    [0.22055347, 0.10920619, 0.27019844],
                    [0.72417325, 0.71259895, 0.77200743],
                    [0.07070797, 0.14246075, 0.31810531],
                    [0.34846727, 0.33212284, 0.09877173],
                    [0.46943827, 0.29961628, 0.43061745],
                ]
            ),
            "string-calc": "this is a string from calc",
        },
    }


def test_encode_data(sample_data):
    encoded = encode_data(sample_data)

    assert encoded == {
        "<SinglePointCalculator>": {
            "energy": 1.0,
            "forces": {
                "data": b"/\xacU\x9a\x18;\xcc?y\x1d\x92\xd6\xef\xf4\xbb?\xe5\xbb\xcee\xeeJ\xd1?\xc2j,am,\xe7?\xfd>-P\x9c\xcd\xe6?\xed\xce\x03\xedH\xb4\xe8?\xc3o\xb7\xe2\xea\x19\xb2?KX\x1bc'<\xc2?\x80\x92\xc8_\xd6[\xd4?\xf8\x16\x18\xaaIM\xd6?;\x7f\x03(\x80A\xd5?4\x8d\x1e\xa6\x1aI\xb9?)\x03I\xd0F\x0b\xde?Y\xbf\xfc\xc2\xe9,\xd3?\x05\x905~<\x8f\xdb?",
                "dtype": "float64",
                "shape": (
                    5,
                    3,
                ),
            },
            "string-calc": "this is a string from calc",
        },
        "REPEATED_KEY": {
            "data": b"\x00\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x03\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00",
            "dtype": "int64",
            "shape": (5,),
        },
        "cell": {
            "data": b"\x00\x00\x00\x00\x00\x00$@\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00$@\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00$@",
            "dtype": "float64",
            "shape": (
                3,
                3,
            ),
        },
        "celldisp": {
            "data": b"\x9a\x99\x99\x99\x99\x99\xb9?\x9a\x99\x99\x99\x99\x99\xc9?333333\xd3?",
            "dtype": "float64",
            "shape": (3,),
        },
        "info": {
            "REPEATED_KEY": "info",
            "array": {
                "data": b"\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00@\x00\x00\x00\x00\x00\x00\x08@",
                "dtype": "float64",
                "shape": (3,),
            },
            "d2": {
                "a": [
                    1,
                    2,
                ],
                "b": [
                    3,
                    4,
                ],
            },
            "dict": {
                "a": 1,
                "b": 2,
            },
            "float": 3.14,
            "list": [
                1,
                2,
                3,
            ],
            "string": "Lorem Ipsum",
        },
        "initial_charges": {
            "data": b"\xc2%tU\xa8\xc7\xe5?Y\xb8r\x17\xdd\xf1\xe8?\x81j\x0f\xbd\xbd\x92\xef?\x83\x8ahC\x177\xbf?P\xc9\xb5x\xe3J\xee?",
            "dtype": "float64",
            "shape": (5,),
        },
        "momenta": {
            "data": b"\xcd\xf8\xc6)\x9d\xbc\xbd?!\xf8J\xe6y\xec\xa2?\x9e=<\xf0\xd6\xde\xdd?\xc0\xbf\xb5\xd1\x08O\xc4?!B\xb7\xf2\xed\xb8\xbd?\xe8\xa3)\x8eE\xab\xea?\xcc\x10d\x90\xfd$\xe3?i\xd5|\xcf'\xb6\xe5?Rzl\xec"
            b'\x8e\xfc\xe7?fe=\xad\x89"\xc3?uh&Ry\xdf\xe5?\xa8\xc0Mo\x14\xb3\xed?\x7f@\x14\x8a|\xf0\xd9?3\x10O\x9f\xb2\x94\xe3?\x8d\xb3+%\xb2\xac\xdb?',
            "dtype": "float64",
            "shape": (
                5,
                3,
            ),
        },
        "name": {
            "data": b"C\x00\x00\x00a\x00\x00\x00r\x00\x00\x00b\x00\x00\x00o\x00\x00\x00n\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00H\x00\x00\x00y\x00\x00\x00d\x00\x00\x00r\x00\x00\x00o\x00\x00\x00g\x00\x00\x00e\x00\x00\x00n\x00\x00\x00H\x00\x00\x00y\x00\x00\x00d\x00\x00\x00r\x00\x00\x00o\x00\x00\x00g\x00\x00\x00e\x00\x00\x00n\x00\x00\x00H\x00\x00\x00y\x00\x00\x00d\x00\x00\x00r\x00\x00\x00o\x00\x00\x00g\x00\x00\x00e\x00\x00\x00n\x00\x00\x00H\x00\x00\x00y\x00\x00\x00d\x00\x00\x00r\x00\x00\x00o\x00\x00\x00g\x00\x00\x00e\x00\x00\x00n\x00\x00\x00",
            "dtype": "<U8",
            "shape": (5,),
        },
        "numbers": {
            "data": b"\x06\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00",
            "dtype": "int64",
            "shape": (5,),
        },
        "pbc": {
            "data": b"\x01\x00\x01",
            "dtype": "bool",
            "shape": (3,),
        },
        "positions": {
            "data": b"\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0?\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xf0\xbf\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00",
            "dtype": "float64",
            "shape": (
                5,
                3,
            ),
        },
        "tags": {
            "data": b"\x01\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x03\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x05\x00\x00\x00\x00\x00\x00\x00",
            "dtype": "int64",
            "shape": (5,),
        },
    }


def test_decode_data(sample_data):
    encoded = encode_data(sample_data)
    decoded = decode_data(encoded)

    assert set(decoded.keys()) == set(sample_data.keys())

    assert_equal(sample_data, decoded)


def test_create_zarr(tmp_path, sample_data):
    root = zarr.group(store=tmp_path / "test.zarr")
    create_zarr(root, sample_data)
    result = read_zarr(root, index=0)
    assert_equal(sample_data, result)


def test_read_zarr(tmp_path, sample_data):
    root = zarr.group(store=tmp_path / "test.zarr")
    create_zarr(root, sample_data)
    result = read_zarr(root, index=0, keys=["numbers", "positions", "info"])
    assert_equal(
        {
            "numbers": sample_data["numbers"],
            "positions": sample_data["positions"],
            "info": sample_data["info"],
        },
        result,
    )
    with pytest.raises(KeyError):
        read_zarr(root, index=0, keys=["nonexistent"])
    with pytest.raises(IndexError):
        read_zarr(root, index=1)


def test_extend_zarr(tmp_path, sample_data):
    root = zarr.group(store=tmp_path / "test.zarr")
    create_zarr(root, sample_data)

    sample_data["<SinglePointCalculator>"]["energy"] = 2.0
    extend_zarr(root, [sample_data])

    assert read_zarr(root, index=0)["<SinglePointCalculator>"]["energy"] == 1.0
    assert read_zarr(root, index=1)["<SinglePointCalculator>"]["energy"] == 2.0


def test_zarr_storage_sequence(tmp_path, sample_data):
    root = zarr.group(store=tmp_path / "test.zarr")
    create_zarr(root, sample_data)
    sequence = ZarrStorageSequence(root)
    assert_equal(sample_data, sequence[0])
    assert len(sequence) == 1

    sequence.append(sample_data)
    assert len(sequence) == 2
    assert_equal(sample_data, sequence[1])


def test_zarr_storage_sequence_create(tmp_path, sample_data):
    root = zarr.group(store=tmp_path / "test.zarr")
    sequence = ZarrStorageSequence(root)
    sequence.append(sample_data)
    assert len(sequence) == 1
    assert_equal(sample_data, sequence[0])
    sequence.get(0, keys=["numbers", "positions"])
    assert len(sequence) == 1
    assert_equal(
        {
            "numbers": sample_data["numbers"],
            "positions": sample_data["positions"],
        },
        sequence.get(0, keys=["numbers", "positions"]),
    )

    assert_equal(sequence.get(0, keys=["info"]), {"info": sample_data["info"]})
