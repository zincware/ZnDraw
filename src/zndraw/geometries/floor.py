"""Floor geometry with grid, shadows, and fog support."""

import typing as t

from pydantic import Field, field_validator

from .base import BaseGeometry


class Floor(BaseGeometry):
    """A floor plane with grid, shadows, and fog support."""

    active: bool = Field(
        default=False,
        description="Whether this geometry should be rendered. Inactive geometries are hidden.",
    )

    size: float = Field(
        default=600.0,
        ge=10.0,
        le=2000.0,
        description="Floor plane size (width and depth)",
    )

    # Override base properties - Floor doesn't need dynamic data
    position: tuple[float, float, float] = Field(
        default=(0, 0, 0), description="Floor center position [x,y,z]"
    )

    @field_validator("position", mode="before")
    @classmethod
    def normalize_position(cls, v: t.Any) -> tuple[float, float, float]:
        """Convert position input to tuple[float, float, float]."""
        if not hasattr(v, "__len__") or len(v) != 3:
            raise ValueError("position must have exactly 3 elements [x, y, z]")
        return (float(v[0]), float(v[1]), float(v[2]))

    color: str = Field(
        default="default",
        description="Floor base color (uses theme colors when 'default')",
    )

    # Floor-specific properties
    height: float = Field(
        default=-5.0, le=50.0, ge=-50.0, description="Y-position of the floor plane"
    )

    grid_spacing: float = Field(
        default=10.0, ge=0.5, le=50.0, description="Spacing between grid lines"
    )

    grid_color: str = Field(
        default="default",
        description="Grid line color (uses theme colors when 'default')",
    )

    grid_opacity: float = Field(
        default=1.0, ge=0.0, le=1.0, description="Grid line opacity"
    )

    show_grid: bool = Field(default=True, description="Toggle grid visibility")

    show_shadows: bool = Field(default=True, description="Toggle contact shadows")

    shadow_opacity: float = Field(
        default=0.5, ge=0.0, le=1.0, description="Shadow opacity"
    )

    shadow_blur: float = Field(
        default=0.5, ge=0.0, le=2.0, description="Shadow blur radius"
    )

    fog_enabled: bool = Field(default=True, description="Enable distance fog effect")

    fog_color: str = Field(
        default="default",
        description="Fog color (uses theme background when 'default')",
    )

    fog_near: float = Field(
        default=180.0,
        ge=0.0,
        le=1000.0,
        description="Distance where fog starts",
    )

    fog_far: float = Field(
        default=300.0,
        ge=10.0,
        le=2000.0,
        description="Distance where fog is fully opaque",
    )

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Remove material property - Floor always uses MeshStandardMaterial
        schema["properties"].pop("material", None)

        # Color pickers
        schema["properties"]["color"]["format"] = "color"
        schema["properties"]["grid_color"]["format"] = "color"

        # Range sliders with steps
        schema["properties"]["grid_spacing"]["format"] = "range"
        schema["properties"]["grid_spacing"]["step"] = 0.5

        schema["properties"]["grid_opacity"]["format"] = "range"
        schema["properties"]["grid_opacity"]["step"] = 0.01

        schema["properties"]["shadow_opacity"]["format"] = "range"
        schema["properties"]["shadow_opacity"]["step"] = 0.01

        schema["properties"]["shadow_blur"]["format"] = "range"
        schema["properties"]["shadow_blur"]["step"] = 0.1

        schema["properties"]["height"]["format"] = "range"
        schema["properties"]["height"]["step"] = 0.5
        schema["properties"]["height"]["minimum"] = -50
        schema["properties"]["height"]["maximum"] = 50

        # Size slider
        schema["properties"]["size"]["format"] = "range"
        schema["properties"]["size"]["step"] = 10

        # Fog color picker
        schema["properties"]["fog_color"]["format"] = "color"

        # Fog distance sliders
        schema["properties"]["fog_near"]["format"] = "range"
        schema["properties"]["fog_near"]["step"] = 10

        schema["properties"]["fog_far"]["format"] = "range"
        schema["properties"]["fog_far"]["step"] = 10

        return schema
