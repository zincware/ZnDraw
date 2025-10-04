import requests

from zndraw import ZnDraw


def test_metadata_s22(server, s22):
    vis = ZnDraw(url=server, room="s22-0", user="user1")
    vis.extend(s22)

    response = requests.get(f"{server}/api/rooms/s22-0/frames/0/metadata")
    assert response.status_code == 200
    metadata = response.json()
    assert metadata == {
        "frameId": 0,
        "keys": [
            "numbers",
            "positions",
            "colors",
            "radii",
            "cell",
            "pbc",
        ],
        "metadata": {
            "cell": {
                "dtype": "float64",
                "shape": [
                    3,
                    3,
                ],
                "type": "array",
            },
            "colors": {
                "dtype": "float32",
                "shape": [
                    8,
                    3,
                ],
                "type": "array",
            },
            "numbers": {
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
            "positions": {
                "dtype": "float64",
                "shape": [
                    8,
                    3,
                ],
                "type": "array",
            },
            "radii": {
                "dtype": "float32",
                "shape": [
                    8,
                ],
                "type": "array",
            },
        },
        "sourceRoom": "s22-0",
    }

    # last frame
    response = requests.get(f"{server}/api/rooms/s22-0/frames/21/metadata")

    assert response.status_code == 200
    metadata = response.json()

    assert metadata == {
        "frameId": 21,
        "keys": [
            "numbers",
            "positions",
            "colors",
            "radii",
            "cell",
            "pbc",
        ],
        "metadata": {
            "cell": {
                "dtype": "float64",
                "shape": [
                    3,
                    3,
                ],
                "type": "array",
            },
            "colors": {
                "dtype": "float32",
                "shape": [
                    26,
                    3,
                ],
                "type": "array",
            },
            "numbers": {
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
            "positions": {
                "dtype": "float64",
                "shape": [
                    26,
                    3,
                ],
                "type": "array",
            },
            "radii": {
                "dtype": "float32",
                "shape": [
                    26,
                ],
                "type": "array",
            },
        },
        "sourceRoom": "s22-0",
    }

    # non-existing frame
    response = requests.get(f"{server}/api/rooms/s22-0/frames/22/metadata")
    assert response.status_code == 404
    assert response.json() == {
        "error": "Invalid frame index 22, valid range: 0-21",
        "type": "IndexError",
    }
