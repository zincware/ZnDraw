"""Light geometry types for scene illumination.

Provides directional, ambient, and hemisphere lights as scene objects.
Lights can be attached to the camera to follow view direction.
"""

from pydantic import BaseModel, ConfigDict, Field

Vec3 = tuple[float, float, float]


class LightPosition(BaseModel):
    """Position of a light - either fixed in world or attached to camera.

    When ``camera_attached`` is True the coordinates are an offset in camera
    space; when False they are absolute world coordinates.
    """

    model_config = ConfigDict(frozen=True)

    camera_attached: bool = Field(
        default=True,
        description="If true, position is relative to camera",
    )
    x: float = Field(default=5.0, description="X coordinate")
    y: float = Field(default=2.0, description="Y coordinate")
    z: float = Field(default=8.0, description="Z coordinate")


class BaseLight(BaseModel):
    """Base class for lights - simpler than BaseGeometry.

    Lights don't need instancing, materials, or selection support.
    """

    model_config = ConfigDict(frozen=True)

    owner: str | None = Field(
        default=None,
        description="User ID of the light owner.",
        json_schema_extra={"x-custom-type": "ownership-toggle"},
    )
    active: bool = Field(
        default=True,
        description="Whether this light is enabled.",
    )
    intensity: float = Field(
        default=1.0,
        ge=0.0,
        le=5.0,
        description="Light intensity multiplier.",
        json_schema_extra={"format": "range", "step": 0.01},
    )


class DirectionalLight(BaseLight):
    """Directional light with parallel rays (like sunlight).

    Position can be fixed in world or attached to camera.
    """

    position: LightPosition = Field(
        default_factory=LightPosition,
        description="Light position - fixed or camera-attached",
    )
    color: str = Field(
        default="#ffffff",
        description="Light color.",
        json_schema_extra={"format": "color"},
    )


class AmbientLight(BaseLight):
    """Uniform ambient light that fills all shadows equally."""

    intensity: float = Field(
        default=0.35,
        ge=0.0,
        le=3.0,
        description="Ambient light intensity.",
        json_schema_extra={"format": "range", "step": 0.01},
    )
    color: str = Field(
        default="#ffffff",
        description="Ambient light color.",
        json_schema_extra={"format": "color"},
    )


class HemisphereLight(BaseLight):
    """Sky/ground gradient ambient light.

    Provides soft lighting with color gradient from sky to ground.
    """

    intensity: float = Field(
        default=0.3,
        ge=0.0,
        le=3.0,
        description="Hemisphere light intensity.",
        json_schema_extra={"format": "range", "step": 0.01},
    )
    sky_color: str = Field(
        default="#e0f7ff",
        description="Color from above (sky).",
        json_schema_extra={"format": "color"},
    )
    ground_color: str = Field(
        default="#ffffff",
        description="Color from below (ground).",
        json_schema_extra={"format": "color"},
    )
