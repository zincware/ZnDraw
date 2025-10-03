import enum
import typing as t

import ase
import numpy as np
from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator
from pydantic.json_schema import SkipJsonSchema

HSLColor = t.Tuple[float, float, float]


class SettingsBase(BaseModel):
    model_config = ConfigDict(validate_assignment=True)

    callback: SkipJsonSchema[t.Callable[[], None] | None] = Field(
        default=None,
        exclude=True,  # 🚀 excludes from model_dump and schema
        repr=False,
    )

    @field_validator("*")
    @classmethod
    def trigger_callback_on_change(cls, v: t.Any, info: ValidationInfo) -> t.Any:
        if "callback" in info.data and callable(info.data["callback"]):
            new_data = info.data
            new_data[info.field_name] = v
            info.data["callback"](
                {k: v for k, v in new_data.items() if k != "callback"}
            )
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


class Representation(SettingsBase):
    """Controls the visual style of atoms and bonds."""
    particle_size: float = Field(1.0, ge=0.1, le=5, description="Atom radius scaling factor")
    show_bonds: bool = Field(True, description="Render bonds between atoms")
    bond_size: float = Field(1.0, ge=0.1, le=5, description="Bond radius scaling factor")
    material: Material = Field(Material.MeshStandardMaterial, description="Atom and bond material")

class Scene(SettingsBase):
    """Controls the background, lighting, and environment."""
    background_color: str = Field("#ffffff", description="Scene background color")
    show_simulation_box: bool = Field(True, description="Show the periodic cell boundary")
    show_floor: bool = Field(False, description="Display a reflective grid floor")
    environment: EnvironmentPreset = Field(
        EnvironmentPreset.studio, 
        description="HDR environment for reflections and lighting"
    )
    # You might add light controls here later, e.g., ambient light intensity

class Playback(SettingsBase):
    """Controls for trajectory animation."""
    animation_loop: bool = Field(False, description="Loop animation when it ends")
    playback_speed: float = Field(1.0, ge=0.1, le=10.0, description="Animation speed multiplier")
    # This was 'frame_update', which is a bit ambiguous. Renaming for clarity.
    sync_with_updates: bool = Field(True, description="Automatically jump to newly added frames")

class Camera(SettingsBase):
    """Defines the camera projection and user interaction controls."""
    camera_type: CameraEnum = Field(
        CameraEnum.PerspectiveCamera, 
        alias="camera", 
        description="Camera projection type"
    )
    controls: Controls = Field(Controls.OrbitControls, description="Mouse interaction mode")
    near_plane: float = Field(0.1, ge=0, le=100, description="Camera near rendering plane")
    far_plane: float = Field(300, ge=1, le=1000, description="Camera far rendering plane")
    show_crosshair: bool = Field(False, description="Show a crosshair at the camera's focal point")
    # This is a collaboration feature, which could even be in its own category if you add more.
    synchronize_view: bool = Field(
        True, 
        description="Synchronize camera with other users in the room"
    )


class PathTracerSettings(SettingsBase):
    """Settings for the experimental path tracer. Overrides standard materials."""
    enabled: bool = False
    metalness: float = Field(0.7, ge=0.0, le=1.0, description="Global metalness")
    roughness: float = Field(0.2, ge=0.0, le=1.0, description="Global roughness")
    clearcoat: float = Field(0.0, ge=0.0, le=1.0, description="Global clearcoat")
    clearcoatRoughness: float = Field(0.0, ge=0.0, le=1.0, description="Global clearcoat roughness")


class Rendering(SettingsBase):
    """Controls for rendering quality and visual effects."""
    max_fps: int = Field(30, ge=1, le=120, description="Maximum frames per second")
    antialiasing: bool = Field(True, description="Enable multisample anti-aliasing (MSAA)")
    ambient_occlusion: bool = Field(False, description="Enable Screen-Space Ambient Occlusion (SSAO)")
    shadows: bool = Field(False, description="Enable dynamic shadows")
    path_tracer: PathTracerSettings = PathTracerSettings()


class Interaction(SettingsBase):
    """Controls for selection, hover effects, and measurements."""
    selection_color: str = Field("#ffa500", description="Highlight color for selected atoms")
    selection_opacity: float = Field(0.5, ge=0.0, le=1.0, description="Opacity of non-selected atoms")
    hover_opacity: float = Field(0.8, ge=0.0, le=1.0, description="Opacity of non-hovered atoms")


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

settings = {
    "camera": Camera,
    "representation": Representation,
    "scene": Scene,
    "playback": Playback,
    "rendering": Rendering,
    "interaction": Interaction,
    "vector_display": VectorDisplay,
}


class RoomConfig(SettingsBase):
    """ZnDraw room configuration combining all settings sections."""

    representation: Representation = Representation()
    scene: Scene = Scene()
    playback: Playback = Playback()
    camera: Camera = Camera()
    rendering: Rendering = Rendering()
    interaction: Interaction = Interaction()
    vector_display: VectorDisplay = VectorDisplay() # Or nest under a general DataOverlays model


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
