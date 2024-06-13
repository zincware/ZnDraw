from zndraw import ZnDraw
import pytest


def test_selection(server, s22):
    """Test the server fixture."""
    vis = ZnDraw(url=server, token="test_token")
    vis.extend(s22)

    assert len(vis.atoms) == 8

    with pytest.raises(ValueError):
        vis.selection = ["Hello"]

    with pytest.raises(IndexError):
        vis.selection = [0, 1, 2, 25]
    
    with pytest.raises(IndexError):
        vis.selection = [0, 1, 2, -10]

    vis.selection = [0, 1, 2, 3]
    assert vis.selection == [0, 1, 2, 3]

# TODO: worker is probably not running
# def test_run_selection(server, s22):
#     """Test the server fixture."""
#     vis = ZnDraw(url=server, token="test_token")
#     vis.extend(s22)
#     vis.selection = [0]
 
#     vis.socket.emit("selection:run", {"method": {"discriminator": "ConnectedParticles"}})
#     vis.socket.sleep(2)
#     assert vis.selection == [0, 1, 2, 3]
