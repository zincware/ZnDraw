import pytest
import numpy as np
import numpy.testing as npt

from zndraw import Box, Sphere, ZnDraw
from zndraw.utils import direction_to_euler, euler_to_direction

def test_geometries(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    assert vis.geometries == []
    vis.geometries = [Box(position=[1, 2, 3]), Sphere()]
    vis.socket.sleep(0.1)

    assert len(vis.geometries) == 2
    assert isinstance(vis.geometries[0], Box)
    assert isinstance(vis.geometries[1], Sphere)

    assert vis.geometries[0].position == [1, 2, 3]
    assert vis.geometries[1].position == [0, 0, 0]

    with pytest.raises(ValueError):
        vis.geometries = ["Geometries are not string!"]

    with pytest.raises(ValueError):
        vis.geometries = "Geometries are not string!"

def test_conversion_utils():
    """ Test conversion functions """
    direction = np.array([1,2,3])
    direction = direction / np.linalg.norm(direction)

    euler = direction_to_euler(direction)
    new_direction = euler_to_direction(euler)

    npt.assert_allclose(direction, new_direction, atol=1e-6)