import time
import typing as t

import pytest
from ase.build import molecule
from pydantic import BaseModel, Field

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


class PerAngstrom(BaseModel):
    discriminator: t.Literal["PerAngstrom"] = "PerAngstrom"
    atoms_per_angstrom: float = Field(
        1.2,
        ge=0,
        le=3.0,
        description="Num atoms added = atoms_per_angstrom * curve_length",
    )


class FixedNumber(BaseModel):
    discriminator: t.Literal["FixedNumber"] = "FixedNumber"
    number_of_atoms: int = Field(
        5, ge=1, le=30, description="Number of atoms to generate"
    )


class Generate(UpdateScene):
    discriminator: t.Literal["Generate"] = Field("Generate")
    num_steps: int = Field(
        50, le=100, ge=20, description="Number of steps in the generation."
    )
    atom_number: t.Union[FixedNumber, PerAngstrom] = FixedNumber(number_of_atoms=5)
    guiding_force_multiplier: float = Field(
        1.0,
        ge=1.0,
        le=10.0,
        description="Multiplier for guiding force. Default value should be enough for simple geometries.",
    )


class Relax(UpdateScene):
    discriminator: t.Literal["Relax"] = Field("Relax")
    max_steps: int = Field(50, ge=1)


class Hydrogenate(UpdateScene):
    discriminator: t.Literal["Hydrogenate"] = Field("Hydrogenate")
    max_steps: int = Field(30, ge=1)


run_types = t.Union[Generate, Relax, Hydrogenate]


class DiffusionModelling(UpdateScene):
    discriminator: t.Literal["DiffusionModelling"] = "DiffusionModelling"
    run_type: run_types = Field(discriminator="discriminator")
    client_address: str = Field("http://127.0.0.1:5000/run")
    path: str = Field(
        "/home/rokas/Programming/MACE-Models",
        description="Path to the repo holding the required models",
    )


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

        vis.register_modifier(DiffusionModelling, default=True)
