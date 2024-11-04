import pytest

from zndraw import ZnDraw, Figure
import matplotlib.pyplot as plt


def test_figures(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    figure = plt.plot([1, 2, 3], [1, 2, 3])
    vis.figures["mtpl"] = figure

    assert len(vis.figures) == 1

