"""Plane geometry for ZnDraw."""

from pydantic import Field

from .base import (
    BaseGeometry,
    ColorProp,
    InteractionSettings,
    PositionProp,
    RotationProp,
    ScaleProp,
    Size2DProp,
)


class Plane(BaseGeometry):
    """Flat rectangular plane geometry with customizable dimensions.

    Planes are defined by a center position and size dimensions [width, height].
    Optional rotation can be applied as Euler angles [x, y, z] in radians.
    By default, planes are oriented in the XY plane (facing along Z-axis).

    By default, creates a single 1x1 plane at the origin (0,0,0).
    """

    position: PositionProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Position coordinates [(x,y,z), ...]. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    color: ColorProp = Field(
        default=["#808080"],
        description="Color values ['#RRGGBB', ...]. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": [
                "color-picker",
                "dynamic-atom-props",
                "free-solo",
                "editable-array",
            ],
        },
    )

    size: Size2DProp = Field(
        default=[(1.0, 1.0)],
        description="Plane dimensions [(width, height), ...]. String for dynamic.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    rotation: RotationProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Rotation as Euler angles [(x, y, z), ...] in radians.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
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
        description="Plane opacity, between 0 (transparent) and 1 (opaque).",
        json_schema_extra={"format": "range", "step": 0.01},
    )

    double_sided: bool = Field(
        default=True,
        description="If true, plane is visible from both sides.",
    )

    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings.",
    )

    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings.",
    )
