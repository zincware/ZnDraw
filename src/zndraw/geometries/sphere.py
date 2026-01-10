"""Sphere geometry for ZnDraw."""

from pydantic import Field

from .base import (
    BaseGeometry,
    InteractionSettings,
    ScaleProp,
    SizeProp,
)


class Sphere(BaseGeometry):
    """A sphere geometry.

    By default, uses dynamic references to arrays.positions and arrays.colors.
    This makes it suitable for visualizing atoms/particles.
    """

    radius: SizeProp = Field(
        default=[1.0],
        description="Sphere radius [r, ...]. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array", "transform"],
        },
    )

    resolution: int = Field(
        default=16,
        ge=4,
        le=64,
        description="Sphere geometry resolution (number of segments).",
        json_schema_extra={"format": "range", "step": 1},
    )

    scale: ScaleProp = Field(
        default=[(1.0, 1.0, 1.0)],
        description="Scale factors [(sx, sy, sz), ...]. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Sphere opacity, between 0 (transparent) and 1 (opaque).",
        json_schema_extra={"format": "range", "step": 0.01},
    )

    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings.",
    )

    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings.",
    )
