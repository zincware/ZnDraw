"""Settings models for ZnDraw configuration.

These Pydantic models define the structure for room-level settings,
generating JSON schemas for the frontend settings panel.
"""

import enum

from pydantic import BaseModel, ConfigDict, Field


class EnvironmentPreset(str, enum.Enum):
    """HDRI environment preset options for path tracing."""

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


class StudioLighting(BaseModel):
    """Controls for the neutral studio lighting setup."""

    model_config = ConfigDict(validate_assignment=True)

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


class PropertyInspector(BaseModel):
    """Property Inspector settings for per-particle and global property display."""

    model_config = ConfigDict(validate_assignment=True)

    enabled_properties: list[str] = Field(
        default_factory=list,
        description="Selected property keys to display in the inspector table",
        json_schema_extra={"x-custom-type": "property-inspector"},
    )


class PathTracing(BaseModel):
    """GPU Path Tracing settings for high-quality physically-based rendering."""

    model_config = ConfigDict(validate_assignment=True)

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


class RoomConfig(BaseModel):
    """ZnDraw room configuration combining all settings sections."""

    model_config = ConfigDict(validate_assignment=True, populate_by_name=True)

    studio_lighting: StudioLighting = Field(default_factory=StudioLighting)
    property_inspector: PropertyInspector = Field(default_factory=PropertyInspector)
    pathtracing: PathTracing = Field(default_factory=PathTracing)
