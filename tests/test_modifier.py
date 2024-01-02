import time
import typing as t

import pytest
from ase.build import molecule

from zndraw import ZnDraw
from zndraw.modify import UpdateScene


def send_raw(vis, event, data):
    msg = {
        "event": event,
        "data": data,
    }
    vis.socket.emit("debug", msg)
    vis.socket.sleep(0.5)


class CustomModifier(UpdateScene):
    discriminator: t.Literal["CustomModifier"] = "CustomModifier"

    def run(self, vis: ZnDraw) -> None:
        # raise ValueError("This is a test")
        vis.append(molecule("H2O"))


@pytest.mark.usefixtures("setup")
class TestZnDrawModifier:
    def test_register_custom_modifier(self, server):
        self.driver.get(server)
        time.sleep(1)
        # we need to wait for all the data to be loaded.
        # This includes jsonschemas and atoms.
        vis = ZnDraw(url=server)
        vis[0] = molecule("H2O")
        assert vis[0] == molecule("H2O")
        assert len(vis) == 1

        vis.register_modifier(CustomModifier, default=True)

        send_raw(
            vis,
            "modifier:run",
            {"params": {"method": {"discriminator": "CustomModifier"}}, "url": server},
        )

        assert len(vis) == 2
