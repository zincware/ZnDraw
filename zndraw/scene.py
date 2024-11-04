import enum
import typing as t

import ase
import numpy as np
from flask import current_app, session
from pydantic import BaseModel, Field

from zndraw.utils import emit_with_retry

if t.TYPE_CHECKING:
    from zndraw import ZnDraw


class Material(str, enum.Enum):
    MeshBasicMaterial = "MeshBasicMaterial"
    # MeshLambertMaterial = "MeshLambertMaterial"
    # MeshMatcapMaterial = "MeshMatcapMaterial"
    # MeshPhongMaterial = "MeshPhongMaterial"
    MeshPhysicalMaterial = "MeshPhysicalMaterial"
    MeshStandardMaterial = "MeshStandardMaterial"
    MeshToonMaterial = "MeshToonMaterial"


class Controls(str, enum.Enum):
    OrbitControls = "OrbitControls"
    TrackballControls = "TrackballControls"


class Camera(str, enum.Enum):
    PerspectiveCamera = "PerspectiveCamera"
    OrthographicCamera = "OrthographicCamera"


# create a class for the material, resolution, etc.
class Scene(BaseModel):
    fps: int = Field(30, ge=1, le=120, description="Maximum frames per second")
    material: Material = Field(Material.MeshStandardMaterial, description="Material")
    # resolution: int = Field(10, ge=1, le=50, description="Resolution")
    particle_size: float = Field(1.0, ge=0.1, le=5, description="Particle Size")
    bond_size: float = Field(1.0, ge=0.1, le=5, description="Bonds Size")
    # wireframe: bool = Field(False, description="Wireframe")
    loop: bool = Field(
        False,
        alias="Animation Loop",
        description="Automatically restart animation when finished.",
    )
    simulation_box: bool = Field(
        True,
        description="Show the simulation box.",
    )
    vectorfield: bool = Field(True, description="Show vectorfield.")

    controls: Controls = Field(Controls.OrbitControls, description="Controls")

    vectors: str = Field("", description="Visualize vectorial property")
    vector_scale: float = Field(1.0, ge=0.1, le=5, description="Rescale Vectors")
    selection_color: str = Field("#ffa500", description="Selection color")
    camera: Camera = Field(Camera.PerspectiveCamera, description="Camera")
    camera_near: float = Field(
        0.1, ge=0.1, le=100, description="Camera near rendering plane"
    )
    camera_far: float = Field(
        500, ge=1, le=1000, description="Camera far rendering plane"
    )
    frame_update: bool = Field(
        True,
        description="Jump to updated frames.",
    )
    crosshair: bool = Field(
        False,
        description="Show camera controls target.",
    )
    floor: bool = Field(
        False,
        description="Show the floor.",
    )
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
        from zndraw import ZnDraw

        room = session["token"]

        vis = ZnDraw(
            r=current_app.extensions["redis"],
            url=current_app.config["SERVER_URL"],
            token=room,
        )
        atoms = vis.atoms
        vis.socket.disconnect()
        # vis.socket.sleep(1)
        return atoms

    @classmethod
    def get_updated_schema(cls) -> dict:
        schema = cls.model_json_schema()

        atoms = cls._get_atoms()
        array_props = [""]
        for key in atoms.calc.results.keys():
            if (
                np.array(atoms.calc.results[key]).ndim == 2
                and np.array(atoms.calc.results[key]).shape[1] == 3
            ):
                array_props.append(key)
        for key in atoms.arrays.keys():
            if (
                np.array(atoms.arrays[key]).ndim == 2
                and np.array(atoms.arrays[key]).shape[1] == 3
            ):
                array_props.append(key)
        schema["properties"]["vectors"]["enum"] = array_props
        schema["properties"]["vectors"]["default"] = ""

        # schema["properties"]["wireframe"]["format"] = "checkbox"
        schema["properties"]["Animation Loop"]["format"] = "checkbox"
        schema["properties"]["simulation_box"]["format"] = "checkbox"
        schema["properties"]["vectorfield"]["format"] = "checkbox"
        schema["properties"]["frame_update"]["format"] = "checkbox"
        schema["properties"]["crosshair"]["format"] = "checkbox"
        schema["properties"]["floor"]["format"] = "checkbox"
        # schema["properties"]["resolution"]["format"] = "range"
        # schema["properties"]["label_offset"]["format"] = "range"
        schema["properties"]["particle_size"]["format"] = "range"
        schema["properties"]["fps"]["format"] = "range"
        schema["properties"]["selection_color"]["format"] = "color"
        schema["properties"]["particle_size"]["step"] = 0.1
        schema["properties"]["bond_size"]["format"] = "range"
        schema["properties"]["bond_size"]["step"] = 0.1
        schema["properties"]["vector_scale"]["format"] = "range"
        schema["properties"]["vector_scale"]["step"] = 0.05

        schema["properties"]["camera_near"]["format"] = "range"
        schema["properties"]["camera_near"]["step"] = 0.1
        schema["properties"]["camera_far"]["format"] = "range"
        schema["properties"]["camera_far"]["step"] = 1
        # schema["properties"]["bonds"]["format"] = "checkbox"
        # schema["properties"]["line_label"]["format"] = "checkbox"

        return schema

    def __setattr__(self, name: str, value) -> None:
        super().__setattr__(name, value)
        if hasattr(self, "_vis") and self._vis is not None:
            # do we still need to check for `vis` ?
            if name != "vis" and name in self.model_fields:
                data = {
                    "scene": self.model_dump(),
                }
                emit_with_retry(
                    self._vis.socket,
                    "room:config:set",
                    data,
                    retries=self._vis.timeout["emit_retries"],
                )
                self._vis.socket.sleep(0.1)  # maybe use call?

    def set_vis(self, vis: "ZnDraw") -> None:
        self._vis = vis
