import multiprocessing as mp
import time

import pytest
from ase.build import molecule

from zndraw import ZnDraw
from zndraw.zndraw import ZnDrawDefault


# @pytest.mark.chrome
@pytest.mark.usefixtures("setup")
class TestZnDraw:
    def test_title(self, server):
        self.driver.get(server)
        assert self.driver.title == "ZnDraw"

    def test_vis(self, server):
        self.driver.get(server)
        proc = mp.Process(
            target=ZnDrawDefault,
            kwargs={"url": server, "token": "default"},
        )
        proc.start()
        time.sleep(5)
        try:
            vis = ZnDraw(url=server)
            # vis.socket.sleep(5)
            assert vis.socket.connected
            for idx in range(10):
                vis[idx] = molecule("H2O")

            assert len(vis) == 10

        finally:
            proc.terminate()
            proc.join()

        # vis[0] = molecule("H2O")
        # assert vis.socket.connected
        # assert vis[0] == molecule("H2O")

        # time.sleep(5)

        # logs = self.driver.get_log("browser")
        # for log_entry in logs:
        #     print(log_entry)
        # raise ValueError(log_entry)

        # wait = WebDriverWait(self.driver, 60)

        # class WaitUntil:
        #     def __init__(self) -> None:
        #         pass

        #     def __call__(self, driver: Any) -> bool:
        #         driver.get(server)
        #         # assert driver.title == "ZnDraw"
        #         # time.sleep(5)
        #         # return True
        #         vis = ZnDraw(url=server)
        #         assert vis.socket.connected
        #         vis[0] = molecule("H2O")
        #         # assert vis.socket.connected
        #         try:
        #             assert vis[0] == molecule("H2O")
        #             return True
        #         except:
        #             return False
        #         # return True

        # wait.until(WaitUntil())

        # assert vis.socket.connected
        # vis[0] = molecule("H2O")
        # assert vis.socket.connected
        # assert vis[0] == molecule("H2O")
