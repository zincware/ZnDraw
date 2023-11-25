
import pytest
from ase.build import molecule

from zndraw import ZnDraw


# @pytest.mark.chrome
@pytest.mark.usefixtures("setup")
class TestZnDraw:
    def test_title(self, server):
        self.driver.get(server)
        assert self.driver.title == "ZnDraw"

    def test_vis_connection(self, server):
        self.driver.get(server)
        vis = ZnDraw(url=server)
        assert vis.socket.connected
    
    def test_vis_len(self, server):
        self.driver.get(server)
        time.sleep(1) 
        # we need to wait for all the data to be loaded.
        # This includes jsonschemas and atoms.
        vis = ZnDraw(url=server)
        assert len(vis) == 1

