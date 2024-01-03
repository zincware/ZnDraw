import time
import typing as t

import pytest
from ase.build import molecule

from zndraw import ZnDraw
from zndraw.modify import UpdateScene
from pydantic import BaseModel, Field


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


class Option1(BaseModel):
    discriminator: t.Literal["Option1"] = "Option1"

class Option2(BaseModel):
    discriminator: t.Literal["Option2"] = "Option2"

class RunType1(UpdateScene):
    discriminator: t.Literal["RunType1"] = Field("RunType1")
    options: t.Union[Option1, Option2] = Option1()

class RunType2(UpdateScene):
    discriminator: t.Literal["RunType2"] = Field("RunType2")

class RunType3(UpdateScene):
    discriminator: t.Literal["RunType3"] = Field("RunType3")

class NestedModifier(UpdateScene):
    discriminator: t.Literal["NestedModifier"] = "NestedModifier"
    run_type: t.Union[RunType1, RunType2, RunType3] = Field(discriminator="discriminator")

    def run(self, vis: ZnDraw) -> None:
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

    def test_register_nested_custom_modifier(self, server):
        self.driver.get(server)
        time.sleep(1)
        # we need to wait for all the data to be loaded.
        # This includes jsonschemas and atoms.
        vis = ZnDraw(url=server)
        vis[0] = molecule("H2O")
        assert vis[0] == molecule("H2O")
        assert len(vis) == 1

        vis.register_modifier(NestedModifier, default=True)

        send_raw(
            vis,
            "modifier:run",
            {"params": {"method": {"discriminator": "NestedModifier"}}, "url": server},
        )

        assert len(vis) == 2
