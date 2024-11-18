import enum
import typing as t

import ase
import numpy as np
from pydantic import BaseModel, Field

if t.TYPE_CHECKING:
    pass


HSLColor = t.Tuple[float, float, float]


class SettingsBase(BaseModel):
    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        return cls.model_json_schema()


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


class CameraEnum(str, enum.Enum):
    PerspectiveCamera = "PerspectiveCamera"
    OrthographicCamera = "OrthographicCamera"


class Particle(SettingsBase):
    particle_size: float = Field(1.0, ge=0.1, le=5, description="Particle Size")
    bond_size: float = Field(1.0, ge=0.1, le=5, description="Bonds Size")
    material: Material = Field(Material.MeshStandardMaterial, description="Material")
    selection_color: str = Field("#ffa500", description="Selection color")
    hover_opacity: float = Field(0.8, ge=0.0, le=1.0, description="Hover opacity")
    selection_opacity: float = Field(0.5, ge=0.0, le=1.0, description="Selection opacity")

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()
        schema["properties"]["particle_size"]["format"] = "range"
        schema["properties"]["particle_size"]["step"] = 0.1
        schema["properties"]["bond_size"]["format"] = "range"
        schema["properties"]["bond_size"]["step"] = 0.1
        schema["properties"]["selection_color"]["format"] = "color"
        schema["properties"]["hover_opacity"]["format"] = "range"
        schema["properties"]["hover_opacity"]["step"] = 0.05
        schema["properties"]["selection_opacity"]["format"] = "range"
        schema["properties"]["selection_opacity"]["step"] = 0.05
        return schema


class VectorDisplay(SettingsBase):
    vectorfield: bool = Field(True, description="Show vectorfield.")
    vectors: list = Field("", description="Visualize vectorial property")
    vector_scale: float = Field(1.0, ge=0.1, le=5, description="Rescale Vectors")

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()
        array_props = []
        if atoms.calc is not None:
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
        # remove "positions" from the list
        array_props = [x for x in array_props if x != "positions"]
        schema["properties"]["vectors"]["items"] = {"type": "string", "enum": array_props}
        schema["properties"]["vectors"]["uniqueItems"] = True
        schema["properties"]["vectorfield"]["format"] = "checkbox"
        schema["properties"]["vector_scale"]["format"] = "range"
        schema["properties"]["vector_scale"]["step"] = 0.05
        return schema


class Visualization(SettingsBase):
    simulation_box: bool = Field(
        True,
        description="Show the simulation box.",
    )
    floor: bool = Field(False, description="Show the floor.")
    frame_update: bool = Field(
        True,
        description="Jump to updated frames.",
    )
    animation_loop: bool = Field(
        False,
        description="Automatically restart animation when finished.",
    )

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()
        schema["properties"]["simulation_box"]["format"] = "checkbox"
        schema["properties"]["frame_update"]["format"] = "checkbox"
        schema["properties"]["animation_loop"]["format"] = "checkbox"
        schema["properties"]["floor"]["format"] = "checkbox"
        return schema


class Camera(SettingsBase):
    camera: CameraEnum = Field(CameraEnum.PerspectiveCamera)
    camera_near: float = Field(
        0.1, ge=0, le=100, description="Camera near rendering plane"
    )
    camera_far: float = Field(
        300, ge=1, le=1000, description="Camera far rendering plane"
    )
    crosshair: bool = Field(
        False,
        description="Show camera controls target.",
    )
    synchronize_camera: bool = Field(
        True,
        description="Synchronize camera with other room members.",
    )
    fps: int = Field(30, ge=1, le=120, description="Maximum frames per second")
    controls: Controls = Field(Controls.OrbitControls, description="Controls")

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()
        schema["properties"]["fps"]["format"] = "range"
        schema["properties"]["fps"]["step"] = 1
        schema["properties"]["camera_near"]["format"] = "range"
        schema["properties"]["camera_near"]["step"] = 0.1
        schema["properties"]["camera_far"]["format"] = "range"
        schema["properties"]["camera_far"]["step"] = 1

        schema["properties"]["crosshair"]["format"] = "checkbox"
        schema["properties"]["synchronize_camera"]["format"] = "checkbox"
        return schema


class Arrows(SettingsBase):
    """Experimental vector color settings."""

    colormap: list[HSLColor] = ((0.5, 0.9, 0.5), (1.0, 0.9, 0.5))
    normalize: bool = True
    colorrange: tuple[float, float] = (0, 1.0)
    scale_vector_thickness: bool = False
    opacity: float = Field(1.0, ge=0.0, le=1.0, description="Opacity")

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()
        schema["properties"]["colormap"] = {
            "type": "array",
            "description": "Defines the colormap as a list of HSL tuples (Hue, Saturation, Lightness).",
            "items": {
                "type": "array",
                "minItems": 3,
                "maxItems": 3,
                "headertemplate": "{{ i1 }}",
                "items": [
                    {
                        "type": "number",
                        "description": "Hue (0.0 - 1.0)",
                        "format": "range",
                        "step": 0.01,
                        "minimum": 0,
                        "maximum": 1,
                    },
                    {
                        "type": "number",
                        "description": "Saturation (0.0 - 1.0)",
                        "format": "range",
                        "step": 0.01,
                        "minimum": 0,
                        "maximum": 1,
                    },
                    {
                        "type": "number",
                        "description": "Lightness (0.0 - 1.0)",
                        "format": "range",
                        "step": 0.01,
                        "minimum": 0,
                        "maximum": 1,
                    },
                ],
            },
        }

        # Enhance "colorrange"
        schema["properties"]["colorrange"] = {
            "type": "array",
            "description": "Specifies the range of values for colors, defined as [min, max].",
            "minItems": 2,
            "maxItems": 2,
            "items": [
                {
                    "type": "number",
                    "description": "Minimum value of the range.",
                    "format": "range",
                    "step": 0.01,
                },
                {
                    "type": "number",
                    "description": "Maximum value of the range.",
                    "format": "range",
                    "step": 0.01,
                },
            ],
        }

        schema["properties"]["normalize"]["format"] = "checkbox"
        schema["properties"]["scale_vector_thickness"]["format"] = "checkbox"
        schema["properties"]["opacity"]["format"] = "range"
        schema["properties"]["opacity"]["step"] = 0.05
        return schema


class EnvironmentPreset(str, enum.Enum):
    none = "none"
    apartment = "apartment"
    city = "city"
    dawn = "dawn"
    forest = "forest"
    lobby = "lobby"
    night = "night"
    park = "park"
    studio = "studio"
    sunset = "sunset"
    warehouse = "warehouse"


class PathTracer(SettingsBase):
    """Experimental path tracer settings."""

    enabled: bool = False
    environment: EnvironmentPreset = EnvironmentPreset.studio
    metalness: float = Field(0.7, ge=0.0, le=1.0, description="Metalness")
    roughness: float = Field(0.2, ge=0.0, le=1.0, description="Roughness")
    clearcoat: float = Field(0.0, ge=0.0, le=1.0, description="Clearcoat")
    clearcoatRoughness: float = Field(
        0.0, ge=0.0, le=1.0, description="Clearcoat Roughness"
    )

    @classmethod
    def model_json_schema_from_atoms(cls, atoms: ase.Atoms) -> dict:
        schema = cls.model_json_schema()
        schema["properties"]["enabled"]["format"] = "checkbox"
        # make all of them sliders
        schema["properties"]["metalness"]["format"] = "range"
        schema["properties"]["roughness"]["format"] = "range"
        schema["properties"]["clearcoat"]["format"] = "range"
        schema["properties"]["clearcoatRoughness"]["format"] = "range"
        # also set the step
        schema["properties"]["metalness"]["step"] = 0.05
        schema["properties"]["roughness"]["step"] = 0.05
        schema["properties"]["clearcoat"]["step"] = 0.05
        schema["properties"]["clearcoatRoughness"]["step"] = 0.05
        return schema


SETTINGS = {
    Particle.__name__: Particle,
    Visualization.__name__: Visualization,
    Camera.__name__: Camera,
    PathTracer.__name__: PathTracer,
    VectorDisplay.__name__: VectorDisplay,
    Arrows.__name__: Arrows,
}
