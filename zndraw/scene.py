import enum

import ase
import znframe
import znsocket
from flask import current_app, session
from pydantic import BaseModel, Field
from redis import Redis


class Material(str, enum.Enum):
    MeshBasicMaterial = "MeshBasicMaterial"
    MeshLambertMaterial = "MeshLambertMaterial"
    MeshMatcapMaterial = "MeshMatcapMaterial"
    MeshPhongMaterial = "MeshPhongMaterial"
    MeshPhysicalMaterial = "MeshPhysicalMaterial"
    MeshStandardMaterial = "MeshStandardMaterial"
    MeshToonMaterial = "MeshToonMaterial"


# create a class for the material, resolution, etc.
class Scene(BaseModel):
    fps: int = Field(30, ge=1, le=120, description="Maximum frames per second")
    # material: Material = Field(Material.MeshPhongMaterial, description="Material")
    # resolution: int = Field(10, ge=1, le=50, description="Resolution")
    # particle_size: float = Field(1.0, ge=0.1, le=5, description="Particle Size")
    # bonds_size: float = Field(1.0, ge=0.1, le=5, description="Bonds Size")
    # wireframe: bool = Field(False, description="Wireframe")
    loop: bool = Field(
        False,
        alias="Animation Loop",
        description="Automatically restart animation when finished.",
    )
    simulation_box: bool = Field(
        False,
        description="Show the simulation box.",
    )

    vectors: str = Field("", description="Visualize vectorial property")
    # bonds: bool = Field(
    #     True,
    #     description="Show bonds.",
    # )
    # line_label: bool = Field(
    #     True,
    #     description="Show the length of the line.",
    # )
    # label_offset: int = Field(
    #     0,
    #     ge=-7,
    #     le=7,
    #     description="Move the label to the left or right (keypress i).",
    # )

    @staticmethod
    def _get_atoms() -> ase.Atoms:
        # TODO: move into utils
        room = session["token"]
        r: Redis = current_app.config["redis"]
        step = r.get(f"room:{room}:step")
        key = (
            f"room:{room}:frames"
            if r.exists(f"room:{room}:frames")
            else "room:default:frames"
        )
        lst = znsocket.List(r, key)
        try:
            frame_json = lst[int(step)]
            return znframe.Frame.from_json(frame_json).to_atoms()
        except TypeError:
            # step is None
            return ase.Atoms()
        except IndexError:
            return ase.Atoms()

    @classmethod
    def updated_schema(cls) -> dict:
        schema = cls.model_json_schema()

        atoms = cls._get_atoms()
        if atoms.calc is not None and "forces" in atoms.calc.results:
            schema["properties"]["vectors"]["enum"] = ["", "forces"]
            schema["properties"]["vectors"]["default"] = ""

        # schema["properties"]["wireframe"]["format"] = "checkbox"
        schema["properties"]["Animation Loop"]["format"] = "checkbox"
        schema["properties"]["simulation_box"]["format"] = "checkbox"
        # schema["properties"]["resolution"]["format"] = "range"
        # schema["properties"]["label_offset"]["format"] = "range"
        # schema["properties"]["particle_size"]["format"] = "range"
        schema["properties"]["fps"]["format"] = "range"
        # schema["properties"]["particle_size"]["step"] = 0.1
        # schema["properties"]["bonds_size"]["format"] = "range"
        # schema["properties"]["bonds_size"]["step"] = 0.1
        # schema["properties"]["bonds"]["format"] = "checkbox"
        # schema["properties"]["line_label"]["format"] = "checkbox"

        return schema
