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
    vis.socket.sleep(1.5)


class CustomModifier(UpdateScene):
    discriminator: t.Literal["CustomModifier"] = "CustomModifier"

    def run(self, vis: ZnDraw) -> None:
        # raise ValueError("This is a test")
        vis.append(molecule("H2O"))

class CustomModifierRunKwargs(UpdateScene):
    discriminator: t.Literal["CustomModifierRunKwargs"] = "CustomModifierRunKwargs"

    def run(self, vis: ZnDraw, structure) -> None:
        # raise ValueError("This is a test")
        vis.append(molecule(structure))

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

    def test_register_custom_modifier_run_kwargs(self, server):
        self.driver.get(server)
        time.sleep(1)
        # we need to wait for all the data to be loaded.
        # This includes jsonschemas and atoms.
        vis = ZnDraw(url=server)
        vis[0] = molecule("H2O")
        assert vis[0] == molecule("H2O")
        assert len(vis) == 1

        vis.register_modifier(CustomModifierRunKwargs, default=True, run_kwargs={"structure": "CH4"})
        assert vis._modifiers["CustomModifierRunKwargs"]["run_kwargs"] == {"structure": "CH4"}
        assert vis._modifiers["CustomModifierRunKwargs"]["cls"] == CustomModifierRunKwargs

        send_raw(
            vis,
            "modifier:run",
            {"params": {"method": {"discriminator": "CustomModifierRunKwargs"}}, "url": server},
        )

        assert len(vis) == 2
        assert vis[0] == molecule("H2O")
        assert vis[1] == molecule("CH4")
