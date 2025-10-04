import numpy as np
import numpy.testing as npt
import pytest
import zarr
from zarr.storage import MemoryStore

from zndraw.storage import (
    ZarrStorageSequence,
    create_zarr,
    decode_data,
    encode_data,
    extend_zarr,
    read_zarr,
)
from zndraw.utils import atoms_to_dict


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


def test_create_zarr(sample_data):
    root = zarr.group(store=MemoryStore())
    create_zarr(root, sample_data)
    result = read_zarr(root, index=0)
    assert_equal(sample_data, result)


def test_read_zarr(sample_data):
    root = zarr.group(store=MemoryStore())
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


def test_extend_zarr(sample_data):
    root = zarr.group(store=MemoryStore())
    create_zarr(root, sample_data)

    sample_data["<SinglePointCalculator>"]["energy"] = 2.0
    extend_zarr(root, [sample_data])

    assert read_zarr(root, index=0)["<SinglePointCalculator>"]["energy"] == 1.0
    assert read_zarr(root, index=1)["<SinglePointCalculator>"]["energy"] == 2.0


def test_zarr_storage_sequence(sample_data):
    root = zarr.group(store=MemoryStore())
    create_zarr(root, sample_data)
    sequence = ZarrStorageSequence(root)
    assert_equal(sample_data, sequence[0])
    assert len(sequence) == 1

    sequence.append(sample_data)
    assert len(sequence) == 2
    assert_equal(sample_data, sequence[1])


def test_zarr_storage_sequence_create(sample_data):
    root = zarr.group(store=MemoryStore())
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


def test_zarr_storage_sequence_add_new():
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)
    store.append({"a": np.array([1, 2])})
    store.append({"a": np.array([3, 4]), "b": np.array([5, 6])})
    store.append({"b": np.array([9, 10])})

    assert len(store) == 3
    assert_equal(store[0], {"a": np.array([1, 2])})
    assert_equal(store.get(0, keys=["a"]), {"a": np.array([1, 2])})
    assert_equal(store[1], {"a": np.array([3, 4]), "b": np.array([5, 6])})
    assert_equal(store.get(1, keys=["a"]), {"a": np.array([3, 4])})
    assert_equal(store.get(1, keys=["b"]), {"b": np.array([5, 6])})
    assert_equal(
        store.get(1, keys=["a", "b"]), {"a": np.array([3, 4]), "b": np.array([5, 6])}
    )
    assert_equal(store[2], {"b": np.array([9, 10])})
    assert_equal(store.get(2, keys=["b"]), {"b": np.array([9, 10])})

    with pytest.raises(KeyError):
        store.get(0, keys=["b"])
    with pytest.raises(KeyError):
        store.get(2, keys=["a"])
    with pytest.raises(KeyError):
        store.get(1, keys=["c"])


def test_zarr_storage_sequence_add_shape():
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)
    store.append({"a": np.array([1, 2])})
    store.append({"a": np.array([3, 4, 5, 6])})
    store.append({"b": np.array([10, 11, 12])})
    store.append({"a": np.array([3, 4, 5])})
    store.append({"b": np.array([13])})

    assert len(store) == 5

    assert_equal(store[0], {"a": np.array([1, 2])})
    assert_equal(store[1], {"a": np.array([3, 4, 5, 6])})
    assert_equal(store[2], {"b": np.array([10, 11, 12])})
    assert_equal(store[3], {"a": np.array([3, 4, 5])})
    assert_equal(store[4], {"b": np.array([13])})


def test_zarr_storage_sequence_add_shape_2d():
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)
    store.append({"a": np.arange(4).reshape(2, 2)})
    store.append({"a": np.arange(6).reshape(3, 2)})

    assert len(store) == 2

    assert_equal(store[0], {"a": np.arange(4).reshape(2, 2)})
    assert_equal(store[1], {"a": np.arange(6).reshape(3, 2)})

    with pytest.raises(ValueError):
        store.append({"a": np.arange(6).reshape(2, 3)})
    with pytest.raises(ValueError):
        store.append(
            {
                "a": np.arange(5).reshape(
                    5,
                )
            }
        )


def test_zarr_storage_large_shape_mismatch():
    """Test for shape mismatch bug when appending arrays with different first dimension sizes."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    # Create initial entry with large array
    store.append({"positions": np.random.rand(3000, 3)})

    # This should work - appending a smaller array
    store.append({"positions": np.random.rand(9, 3)})

    assert len(store) == 2
    assert store[0]["positions"].shape == (3000, 3)
    assert store[1]["positions"].shape == (9, 3)


def test_zarr_storage_variable_sizes():
    """Test that atoms objects with different particle counts are stored and retrieved correctly."""
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    dicts = [
        {"x": np.array([1, 2, 3])},
        {"x": np.array([1, 2, 3, 4])},
        {"x": np.array([1, 2])},
        {"x": np.array([1, 2, 3, 4, 5, 6])},
        {"x": np.array([1])},
    ]
    store.extend(dicts)

    assert len(store) == len(dicts)
    npt.assert_array_equal(store[0]["x"], dicts[0]["x"])
    npt.assert_array_equal(store[1]["x"], dicts[1]["x"])
    npt.assert_array_equal(store[2]["x"], dicts[2]["x"])
    npt.assert_array_equal(store[3]["x"], dicts[3]["x"])
    npt.assert_array_equal(store[4]["x"], dicts[4]["x"])


def test_zarr_storage_variable_sized_atoms(s22):
    dict = [atoms_to_dict(atoms) for atoms in s22]
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)
    store.extend(dict)

    assert len(store) == len(s22)
    for i in range(len(s22)):
        assert_equal(dict[i], store[i])


def test_zarr_storage_empty():
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    assert len(store) == 0
    with pytest.raises(IndexError):
        store[0]


def test_zarr_storage_append_empty_dict():
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    store.append({"1D": np.random.rand(3), "2D": np.random.rand(4, 3)})
    store.append({"1D": np.random.rand(0), "2D": np.random.rand(0, 3)})
    store.append({"1D": np.random.rand(2), "2D": np.random.rand(2, 3)})

    assert len(store) == 3
    assert store[0]["1D"].shape == (3,)
    assert store[1]["1D"].shape == (0,)
    assert store[2]["1D"].shape == (2,)
    assert store[0]["2D"].shape == (4, 3)
    assert store[1]["2D"].shape == (0, 3)
    assert store[2]["2D"].shape == (2, 3)


def test_zarr_slicing():
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)
    for idx in range(10):
        store.append({"1D": np.array(idx)})
    assert len(store) == 10
    assert store[0]["1D"] == 0
    assert store[-1]["1D"] == 9
    assert store[np.array(-1)]["1D"] == 9
    npt.assert_array_equal(store[2:5]["1D"], np.array([2, 3, 4]))
    npt.assert_array_equal(store[:3]["1D"], np.array([0, 1, 2]))
    npt.assert_array_equal(store[7:]["1D"], np.array([7, 8, 9]))
    npt.assert_array_equal(store[::2]["1D"], np.array([0, 2, 4, 6, 8]))
    npt.assert_array_equal(store[::-1]["1D"], np.array([9, 8, 7, 6, 5, 4, 3, 2, 1, 0]))
    npt.assert_array_equal(store[-3:]["1D"], np.array([7, 8, 9]))
    npt.assert_array_equal(store[-5:-2]["1D"], np.array([5, 6, 7]))
    npt.assert_array_equal(store[-1:-6:-1]["1D"], np.array([9, 8, 7, 6, 5]))

    npt.assert_array_equal(store[[1, 3, 5]]["1D"], np.array([1, 3, 5]))
    npt.assert_array_equal(store[[1, 3, -2]]["1D"], np.array([1, 3, 8]))
    npt.assert_array_equal(store[np.array([1, 3, -2])]["1D"], np.array([1, 3, 8]))

    with pytest.raises(IndexError):
        store[10]
    with pytest.raises(IndexError):
        store[-11]

    with pytest.raises(IndexError):
        store[np.array([1, 2, 100])]


def test_zarr_storage_mask_resize_on_matching_shape():
    """Test that mask arrays are properly resized when appending data with matching shapes.

    This test captures a bug where the mask array wasn't resized when appending
    entries with the same shape as the current max shape, causing index out of bounds errors.
    """
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    # First, create entries with different sizes to create a mask array
    store.append({"x": np.array([1, 2, 3])})
    store.append({"x": np.array([4, 5])})  # Different size, creates mask

    # Now append more entries with size 3 (matches current array max shape)
    # Before the fix, this would fail with: zarr.errors.BoundsCheckError:
    # index out of bounds for dimension with length 2
    store.append({"x": np.array([6, 7, 8])})  # Size 3 - matches current array shape
    store.append({"x": np.array([9, 10, 11])})  # Size 3 - should still work

    # Verify all entries are stored correctly
    assert len(store) == 4
    npt.assert_array_equal(store[0]["x"], np.array([1, 2, 3]))
    npt.assert_array_equal(store[1]["x"], np.array([4, 5]))
    npt.assert_array_equal(store[2]["x"], np.array([6, 7, 8]))
    npt.assert_array_equal(store[3]["x"], np.array([9, 10, 11]))





def test_numpy_scalar_types_in_storage():
    """Test that numpy scalar types (int64, float64, bool_) are properly handled in zarr storage.

    This test verifies the fix for the 'can not serialize numpy.int64 object' error
    that occurred when reading files with ASE, which produces numpy scalar types.
    """
    root = zarr.group(store=MemoryStore())
    store = ZarrStorageSequence(root)

    # Create data with numpy scalar types (as produced by ASE's atoms.todict())
    data_with_numpy_scalars = {
        "positions": np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
        "numbers": np.array([1, 6]),
        "info": {
            "energy": np.float64(123.456),  # numpy scalar
            "step": np.int64(100),  # numpy scalar
            "converged": np.bool_(True),  # numpy scalar
            "iterations": np.int32(42),  # numpy scalar
            "nested": {
                "value": np.int64(999),  # nested numpy scalar
                "factor": np.float32(2.5),  # numpy scalar
            },
        },
        "metadata": {
            "temperature": np.float64(300.0),  # numpy scalar in different dict
            "pressure": np.float32(1.0),  # numpy scalar
        },
    }

    # This should not raise TypeError about serializing numpy.int64
    store.append(data_with_numpy_scalars)

    # Verify we can read it back
    retrieved = store[0]

    # Verify arrays are preserved
    npt.assert_array_equal(retrieved["positions"], data_with_numpy_scalars["positions"])
    npt.assert_array_equal(retrieved["numbers"], data_with_numpy_scalars["numbers"])

    # Verify scalar values are preserved (converted to Python types but values match)
    assert retrieved["info"]["energy"] == 123.456
    assert retrieved["info"]["step"] == 100
    assert retrieved["info"]["converged"] is True
    assert retrieved["info"]["iterations"] == 42
    assert retrieved["info"]["nested"]["value"] == 999
    assert retrieved["info"]["nested"]["factor"] == pytest.approx(2.5)
    assert retrieved["metadata"]["temperature"] == 300.0
    assert retrieved["metadata"]["pressure"] == pytest.approx(1.0)

    # Verify multiple frames work
    data2 = {
        "positions": np.array([[2.0, 2.0, 2.0]]),
        "numbers": np.array([8]),
        "info": {
            "energy": np.float64(456.789),
            "step": np.int64(200),
        },
    }
    store.append(data2)

    assert len(store) == 2
    retrieved2 = store[1]
    assert retrieved2["info"]["energy"] == 456.789
    assert retrieved2["info"]["step"] == 200
