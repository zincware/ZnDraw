import numpy as np
import requests
from ase.calculators.singlepoint import SinglePointCalculator

from zndraw import ZnDraw


def test_frame_keys_s22(server, s22):
    """Test the lightweight /keys endpoint."""
    vis = ZnDraw(url=server, room="s22-0", user="user1")
    vis.extend(s22)

    response = requests.get(f"{server}/api/rooms/s22-0/frames/0/keys")
    assert response.status_code == 200
    data = response.json()

    # Check structure
    assert "frameId" in data
    assert "sourceRoom" in data
    assert "keys" in data
    assert "metadata" not in data  # Keys endpoint doesn't include metadata

    # Check values
    assert data["frameId"] == 0
    assert data["sourceRoom"] == "s22-0"
    assert set(data["keys"]) == {
        "arrays.numbers",
        "arrays.positions",
        "arrays.colors",
        "arrays.radii",
        "cell",
        "pbc",
        "info.connectivity",
    }

    # Test last frame
    response = requests.get(f"{server}/api/rooms/s22-0/frames/21/keys")
    assert response.status_code == 200
    data = response.json()
    assert data["frameId"] == 21
    assert set(data["keys"]) == {
        "arrays.numbers",
        "arrays.positions",
        "arrays.colors",
        "arrays.radii",
        "cell",
        "pbc",
        "info.connectivity",
    }

    # Test non-existing frame
    response = requests.get(f"{server}/api/rooms/s22-0/frames/22/keys")
    assert response.status_code == 404
    assert response.json() == {
        "error": "Invalid frame index 22, valid range: 0-21",
        "type": "IndexError",
    }


def test_metadata_s22(server, s22):
    vis = ZnDraw(url=server, room="s22-0", user="user1")
    vis.extend(s22)

    response = requests.get(f"{server}/api/rooms/s22-0/frames/0/metadata")
    assert response.status_code == 200
    metadata = response.json()

    # Check structure and keys (order doesn't matter)
    assert metadata["frameId"] == 0
    assert metadata["sourceRoom"] == "s22-0"
    assert set(metadata["keys"]) == {
        "arrays.numbers",
        "arrays.positions",
        "arrays.colors",
        "arrays.radii",
        "cell",
        "pbc",
        "info.connectivity",
    }

    # Check metadata details
    expected_metadata = {
            "cell": {
                "dtype": "float64",
                "shape": [
                    3,
                    3,
                ],
                "type": "array",
            },
            "arrays.colors": {
                "dtype": "json",
                "type": "json",
            },
            "arrays.numbers": {
                "dtype": "int64",
                "shape": [
                    8,
                ],
                "type": "array",
            },
            "pbc": {
                "dtype": "bool",
                "shape": [
                    3,
                ],
                "type": "array",
            },
            "arrays.positions": {
                "dtype": "float64",
                "shape": [
                    8,
                    3,
                ],
                "type": "array",
            },
            "arrays.radii": {
                "dtype": "float32",
                "shape": [
                    8,
                ],
                "type": "array",
            },
        }

    # info.connectivity is added automatically, check it exists and has correct type
    assert "info.connectivity" in metadata["metadata"]
    assert metadata["metadata"]["info.connectivity"]["dtype"] == "int32"
    assert metadata["metadata"]["info.connectivity"]["type"] == "array"
    assert len(metadata["metadata"]["info.connectivity"]["shape"]) == 2  # [n_bonds, 3]

    # Check other metadata matches expected
    for key, value in expected_metadata.items():
        assert metadata["metadata"][key] == value, f"Metadata mismatch for {key}"

    # last frame
    response = requests.get(f"{server}/api/rooms/s22-0/frames/21/metadata")

    assert response.status_code == 200
    metadata = response.json()

    assert metadata["frameId"] == 21
    assert metadata["sourceRoom"] == "s22-0"
    assert set(metadata["keys"]) == {
        "arrays.numbers",
        "arrays.positions",
        "arrays.colors",
        "arrays.radii",
        "cell",
        "pbc",
        "info.connectivity",
    }

    expected_metadata_frame21 = {
            "cell": {
                "dtype": "float64",
                "shape": [
                    3,
                    3,
                ],
                "type": "array",
            },
            "arrays.colors": {
                "dtype": "json",
                "type": "json",
            },
            "arrays.numbers": {
                "dtype": "int64",
                "shape": [
                    26,
                ],
                "type": "array",
            },
            "pbc": {
                "dtype": "bool",
                "shape": [
                    3,
                ],
                "type": "array",
            },
            "arrays.positions": {
                "dtype": "float64",
                "shape": [
                    26,
                    3,
                ],
                "type": "array",
            },
            "arrays.radii": {
                "dtype": "float32",
                "shape": [
                    26,
                ],
                "type": "array",
            },
        }

    # info.connectivity is added automatically, check it exists and has correct type
    assert "info.connectivity" in metadata["metadata"]
    assert metadata["metadata"]["info.connectivity"]["dtype"] == "int32"
    assert metadata["metadata"]["info.connectivity"]["type"] == "array"
    assert len(metadata["metadata"]["info.connectivity"]["shape"]) == 2  # [n_bonds, 3]

    # Check other metadata matches expected
    for key, value in expected_metadata_frame21.items():
        assert metadata["metadata"][key] == value, f"Metadata mismatch for {key}"

    # non-existing frame
    response = requests.get(f"{server}/api/rooms/s22-0/frames/22/metadata")
    assert response.status_code == 404
    assert response.json() == {
        "error": "Invalid frame index 22, valid range: 0-21",
        "type": "IndexError",
    }


def test_metadata_s22_arrays(server, s22):
    vis = ZnDraw(url=server, room="s22-0", user="user1")
    for atoms in s22:
        atoms.arrays["forces"] = np.random.rand(len(atoms), 3)
    vis.extend(s22)

    response = requests.get(f"{server}/api/rooms/s22-0/frames/0/metadata")
    assert response.status_code == 200
    metadata = response.json()

    assert set(metadata["keys"]) == {
        "arrays.numbers",
        "arrays.positions",
        "arrays.colors",
        "arrays.radii",
        "cell",
        "pbc",
        "info.connectivity",
        "arrays.forces",
    }

    assert metadata["metadata"]["arrays.forces"] == {
        "dtype": "float64",
        "shape": [8, 3],
        "type": "array",
    }


def test_metadata_s22_info(server, s22):
    vis = ZnDraw(url=server, room="s22-0", user="user1")
    for atoms in s22:
        atoms.info["energy"] = np.random.rand()
    vis.extend(s22)

    response = requests.get(f"{server}/api/rooms/s22-0/frames/0/metadata")
    assert response.status_code == 200
    metadata = response.json()

    assert set(metadata["keys"]) == {
        "arrays.numbers",
        "arrays.positions",
        "arrays.colors",
        "arrays.radii",
        "cell",
        "pbc",
        "info.connectivity",
        "info.energy",
    }

    assert metadata["metadata"]["info.energy"] == {
        "dtype": "float64",
        "shape": [],
        "type": "array",
    }


def test_metadata_s22_calc(server, s22):
    vis = ZnDraw(url=server, room="s22-0", user="user1")
    for atoms in s22:
        atoms.calc = SinglePointCalculator(
            atoms, energy=np.random.rand(), forces=np.random.rand(len(atoms), 3)
        )
    vis.extend(s22)

    response = requests.get(f"{server}/api/rooms/s22-0/frames/0/metadata")
    assert response.status_code == 200
    metadata = response.json()

    assert set(metadata["keys"]) == {
        "arrays.numbers",
        "arrays.positions",
        "arrays.colors",
        "arrays.radii",
        "cell",
        "pbc",
        "info.connectivity",
        "calc.energy",
        "calc.forces",
    }

    assert metadata["metadata"]["calc.energy"] == {
        "dtype": "float64",
        "shape": [],
        "type": "array",
    }

    assert metadata["metadata"]["calc.forces"] == {
        "dtype": "float64",
        "shape": [8, 3],
        "type": "array",
    }


def test_metadata_s22_info_arrays_calc(server, s22):
    vis = ZnDraw(url=server, room="s22-0", user="user1")
    for atoms in s22:
        atoms.calc = SinglePointCalculator(
            atoms, energy=np.random.rand(), forces=np.random.rand(len(atoms), 3)
        )
        atoms.arrays["forces"] = np.random.rand(len(atoms), 3)
        atoms.info["energy"] = np.random.rand()
    vis.extend(s22)
    response = requests.get(f"{server}/api/rooms/s22-0/frames/0/metadata")
    assert response.status_code == 200
    metadata = response.json()

    assert set(metadata["keys"]) == {
        "arrays.numbers",
        "arrays.positions",
        "arrays.colors",
        "arrays.radii",
        "cell",
        "pbc",
        "info.connectivity",
        "arrays.forces",
        "info.energy",
        "calc.energy",
        "calc.forces",
    }

    assert metadata["metadata"]["arrays.forces"] == {
        "dtype": "float64",
        "shape": [8, 3],
        "type": "array",
    }
    assert metadata["metadata"]["info.energy"] == {
        "dtype": "float64",
        "shape": [],
        "type": "array",
    }
    assert metadata["metadata"]["calc.energy"] == {
        "dtype": "float64",
        "shape": [],
        "type": "array",
    }
    assert metadata["metadata"]["calc.forces"] == {
        "dtype": "float64",
        "shape": [8, 3],
        "type": "array",
    }
