import numpy as np
import numpy.testing as npt

from zndraw import ZnDraw


def test_points_segments(server):
    """Test the server fixture."""
    zndraw = ZnDraw(url=server, token="test_token")

    assert isinstance(zndraw.points, np.ndarray)
    assert isinstance(zndraw.segments, np.ndarray)

    npt.assert_array_equal(zndraw.points, [])
    npt.assert_array_equal(zndraw.segments, [])

    zndraw.points = np.array([[0, 0, 0], [1, 1, 1]])
    zndraw.socket.sleep(0.5)
    assert zndraw.points.shape == (2, 3)
    assert zndraw.segments.shape == (100, 3)

    npt.assert_array_equal(zndraw.points, [[0, 0, 0], [1, 1, 1]])
