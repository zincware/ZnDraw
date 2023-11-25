
import pytest
from ase.build import molecule

from zndraw import ZnDraw


# @pytest.mark.chrome
@pytest.mark.usefixtures("setup")
class TestZnDraw:
    def test_title(self, server):
        self.driver.get(server)
        assert self.driver.title == "ZnDraw"

    def test_vis(self, server):
        self.driver.get(server)
        vis = ZnDraw(url=server)
        assert vis.socket.connected
        for idx in range(10):
            vis[idx] = molecule("H2O")

        assert len(vis) == 10
