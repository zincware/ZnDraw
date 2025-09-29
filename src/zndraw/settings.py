import enum
import typing as t

import ase
import numpy as np
from pydantic import BaseModel, Field

from pydantic import field_validator, ValidationInfo, ConfigDict
from pydantic.json_schema import SkipJsonSchema

HSLColor = t.Tuple[float, float, float]


class SettingsBase(BaseModel):
    model_config = ConfigDict(validate_assignment=True)

    callback: SkipJsonSchema[t.Callable[[], None] | None] = Field(
        default=None,
        exclude=True,  # 🚀 excludes from model_dump and schema
    )

    @field_validator("*")
    @classmethod
    def trigger_callback_on_change(cls, v: t.Any, info: ValidationInfo) -> t.Any:
        if "callback" in info.data and callable(info.data["callback"]):
            new_data = info.data
            new_data[info.field_name] = v
            info.data["callback"]({k: v for k, v in new_data.items() if k != "callback"})
        return v


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
    show_bonds: bool = Field(True, description="Show bonds")
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
        schema["properties"]["show_bonds"]["format"] = "checkbox"
        schema["properties"]["selection_color"]["format"] = "color"
        schema["properties"]["hover_opacity"]["format"] = "range"
        schema["properties"]["hover_opacity"]["step"] = 0.05
        schema["properties"]["selection_opacity"]["format"] = "range"
        schema["properties"]["selection_opacity"]["step"] = 0.05
        return schema


class VectorDisplay(SettingsBase):
    vectorfield: bool = Field(True, description="Show vectorfield.")
    vectors: list[str] = Field(
        default_factory=list, description="Visualize vectorial property"
    )
    vector_scale: float = Field(1.0, ge=0.1, le=5, description="Rescale Vectors")
    vector_colors: dict[str, str] = Field(
        default_factory=dict, description="Color for each vector"
    )

    # Arrow configuration properties
    normalize: bool = Field(True, description="Normalize vector lengths")
    scale_vector_thickness: bool = Field(
        False, description="Scale arrow thickness based on vector magnitude"
    )
    opacity: float = Field(1.0, ge=0.0, le=1.0, description="Arrow transparency")
    colorrange: tuple[float, float] = Field(
        (0.0, 1.0), description="Min and max values for color mapping"
    )
    default_colormap: list[tuple[float, float, float]] = Field(
        default=[(0.66, 1.0, 0.5), (0.0, 1.0, 0.5)],
        description="Default HSL colormap for vector fields (blue to red)",
    )

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

        # Configure vector_colors as color pickers for each available vector
        schema["properties"]["vector_colors"] = {
            "type": "object",
            "description": "Color for each vector",
            "additionalProperties": {
                "type": "string",
                "format": "color",
                "default": "#ff0000",
            },
        }
        if array_props:
            schema["properties"]["vector_colors"]["properties"] = {
                prop: {"type": "string", "format": "color", "default": "#ff0000"}
                for prop in array_props
            }

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


settings = {
    "particle": Particle,
    "visualization": Visualization,
    "camera": Camera,
    "path_tracer": PathTracer,
    "vector_display": VectorDisplay,
}


class RoomConfig(SettingsBase):
    """ZnDraw room configuration combining all settings sections."""

    particle: Particle = Particle()
    visualization: Visualization = Visualization()
    camera: Camera = Camera()
    path_tracer: PathTracer = PathTracer()
    vector_display: VectorDisplay = VectorDisplay()



if __name__ == "__main__":
    import json
    import subprocess
    from pathlib import Path

    # Create output directories if they don't exist
    backend_output_dir = Path(__file__).parent / "generated"
    frontend_output_dir = Path(__file__).parent.parent / "app" / "src" / "types"
    backend_output_dir.mkdir(exist_ok=True)
    frontend_output_dir.mkdir(exist_ok=True)

    # Generate schema using pydantic model (includes default values)
    schema = RoomConfig.model_json_schema()

    # Export JSON schema
    schema_file = backend_output_dir / "config.json"
    with open(schema_file, "w") as f:
        json.dump(schema, f, indent=2)

    # Convert to TypeScript using bunx quicktype with JSON schema
    ts_file = frontend_output_dir / "room-config.ts"
    subprocess.run(
        [
            "bunx",
            "quicktype",
            "--lang",
            "typescript",
            "--src-lang",
            "schema",
            "--src",
            str(schema_file),
            "--out",
            str(ts_file),
            "--top-level",
            "RoomConfig",
            "--prefer-unions",
            "--just-types",
        ],
        check=True,
    )

    # Create default values object from actual pydantic instances
    room_config_instance = RoomConfig(
        Particle=Particle(),
        Visualization=Visualization(),
        Camera=Camera(),
        PathTracer=PathTracer(),
        VectorDisplay=VectorDisplay(),
    )
    defaults = room_config_instance.model_dump()

    # Add header comment and default values to the TypeScript file
    with open(ts_file, "r") as f:
        content = f.read()

    # Add header comment at the top
    header_comment = """/**
 * ZnDraw room configuration combining all settings sections.
 *
 * This file is automatically generated by running: python -m zndraw.config
 * Do not edit manually - changes will be overwritten.
 */
"""

    # Replace the content with header + existing content + defaults
    with open(ts_file, "w") as f:
        f.write(header_comment)
        f.write(content)
        f.write("\n\n// Default configuration values from pydantic models\n")
        f.write("export const DEFAULT_ROOM_CONFIG: RoomConfig = ")
        f.write(json.dumps(defaults, indent=2))
        f.write(";\n")

    print(f"Generated JSON schema: {schema_file}")
    print(f"Generated TypeScript: {ts_file}")
    print("Frontend types available at: app/src/types/room-config.ts")