import matplotlib.pyplot as plt

from zndraw import Figure, ZnDraw


def test_figures(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [1, 2, 3])

    vis.figures["mtpl"] = Figure(figure=fig)

    assert len(vis.figures) == 1
