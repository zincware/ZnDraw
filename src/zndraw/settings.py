import enum
import typing as t

from pydantic import BaseModel, ConfigDict, Field, ValidationInfo, field_validator
from pydantic.json_schema import SkipJsonSchema

HSLColor = t.Tuple[float, float, float]


class SettingsBase(BaseModel):
    model_config = ConfigDict(validate_assignment=True)

    callback: SkipJsonSchema[t.Callable[[dict], None] | None] = Field(
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


class StudioLighting(SettingsBase):
    """Controls for the neutral studio lighting setup."""

    background_color: str = Field(
        default="default",
        description="Neutral background color of the scene",
        json_schema_extra={"format": "color"},
    )
    key_light: float = Field(
        default=0.7,
        ge=0.0,
        le=3.0,
        description="Intensity of the main light attached to the camera",
        json_schema_extra={"format": "range", "step": 0.01},
    )
    fill_light: float = Field(
        default=0.4,
        ge=0.0,
        le=3.0,
        description="Intensity of the soft global light that lifts shadows",
        json_schema_extra={"format": "range", "step": 0.01},
    )
    rim_light: float = Field(
        default=0.5,
        ge=0.0,
        le=5.0,
        description="Intensity of the back light that creates highlights",
        json_schema_extra={"format": "range", "step": 0.01},
    )
    hemisphere_light: float = Field(
        default=0.3,
        ge=0.0,
        le=3.0,
        description="Intensity of the ambient light from above",
        json_schema_extra={"format": "range", "step": 0.01},
    )
    ambient_light: float = Field(
        default=0.35,
        ge=0.0,
        le=3.0,
        description="Intensity of the ambient light that fills the scene",
        json_schema_extra={"format": "range", "step": 0.01},
    )
    contact_shadow: bool = Field(
        default=False, description="Show contact shadow below the model"
    )


class PropertyInspector(SettingsBase):
    """Property Inspector settings for per-particle and global property display."""

    enabled_properties: list[str] = Field(
        default_factory=list,
        description="Selected property keys to display in the inspector table",
        json_schema_extra={"x-custom-type": "property-inspector"},
    )


class PathTracing(SettingsBase):
    """GPU Path Tracing settings for high-quality physically-based rendering."""

    enabled: bool = Field(default=False, description="Enable GPU path tracing renderer")

    min_samples: int = Field(
        default=1,
        ge=1,
        description="Minimum samples before displaying result",
        json_schema_extra={"format": "range"},
    )

    samples: int = Field(
        default=256,
        ge=1,
        le=10000,
        description="Maximum samples to render",
        json_schema_extra={"format": "range"},
    )

    bounces: int = Field(
        default=3,
        ge=1,
        le=32,
        description="Number of light bounces for global illumination",
        json_schema_extra={"format": "range"},
    )

    tiles: int = Field(
        default=1,
        ge=1,
        le=8,
        description="Rendering tile count (higher = less memory, slower)",
        json_schema_extra={"format": "range"},
    )

    environment_preset: EnvironmentPreset = Field(
        default=EnvironmentPreset.studio,
        description="HDRI environment preset for scene lighting",
    )

    environment_intensity: float = Field(
        default=1.0,
        ge=0.0,
        le=10.0,
        description="Environment map brightness multiplier",
        json_schema_extra={"format": "range", "step": 0.1},
    )

    environment_blur: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
        description="Environment background blur amount",
        json_schema_extra={"format": "range", "step": 0.01},
    )

    environment_background: bool = Field(
        default=False, description="Show environment as visible background"
    )


class RoomConfig(SettingsBase):
    """ZnDraw room configuration combining all settings sections."""

    model_config = ConfigDict(validate_assignment=True, populate_by_name=True)

    studio_lighting: StudioLighting = Field(default_factory=StudioLighting)
    property_inspector: PropertyInspector = Field(default_factory=PropertyInspector)
    pathtracing: PathTracing = Field(default_factory=PathTracing)


if __name__ == "__main__":
    import json
    import subprocess
    from pathlib import Path

    # Create output directories if they don't exist
    backend_output_dir = Path(__file__).parent / "generated"
    frontend_output_dir = Path(__file__).parent.parent.parent / "app" / "src" / "types"
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
    room_config_instance = RoomConfig()
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
        # f.write("\n\n// Default configuration values from pydantic models\n")
        # f.write("export const DEFAULT_ROOM_CONFIG: RoomConfig = ")
        # f.write(json.dumps(defaults, indent=2))
        # f.write(";\n")

    print(f"Generated JSON schema: {schema_file}")
    print(f"Generated TypeScript: {ts_file}")
    print("Frontend types available at: app/src/types/room-config.ts")
