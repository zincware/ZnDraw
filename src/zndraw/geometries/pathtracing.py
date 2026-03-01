"""Path tracing renderer configuration as a scene object."""

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


class PathTracing(BaseModel):
    """GPU path tracing renderer configuration.

    When active, switches from standard WebGL renderer to GPU path tracer
    for physically-based rendering with global illumination.
    """

    model_config = ConfigDict(frozen=True)

    owner: str | None = Field(
        default=None,
        description="User ID of the pathtracing config owner.",
        json_schema_extra={"x-custom-type": "ownership-toggle"},
    )
    active: bool = Field(
        default=False,
        description="Enable GPU path tracing renderer.",
    )
    min_samples: int = Field(
        default=1,
        ge=1,
        description="Minimum samples before displaying result.",
        json_schema_extra={"format": "range"},
    )
    samples: int = Field(
        default=256,
        ge=1,
        le=10000,
        description="Maximum samples to render.",
        json_schema_extra={"format": "range"},
    )
    bounces: int = Field(
        default=3,
        ge=1,
        le=32,
        description="Number of light bounces for global illumination.",
        json_schema_extra={"format": "range"},
    )
    tiles: int = Field(
        default=1,
        ge=1,
        le=8,
        description="Rendering tile count (higher = less memory, slower).",
        json_schema_extra={"format": "range"},
    )
    environment_preset: EnvironmentPreset = Field(
        default=EnvironmentPreset.studio,
        description="HDRI environment preset for scene lighting.",
    )
    environment_intensity: float = Field(
        default=1.0,
        ge=0.0,
        le=10.0,
        description="Environment map brightness multiplier.",
        json_schema_extra={"format": "range", "step": 0.1},
    )
    environment_blur: float = Field(
        default=0.0,
        ge=0.0,
        le=1.0,
        description="Environment background blur amount.",
        json_schema_extra={"format": "range", "step": 0.01},
    )
    environment_background: bool = Field(
        default=False,
        description="Show environment as visible background.",
    )
