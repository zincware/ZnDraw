"""Box geometry for ZnDraw."""

from pydantic import Field

from .base import (
    BaseGeometry,
    ColorProp,
    InteractionSettings,
    PositionProp,
    RotationProp,
    ScaleProp,
    Size3DProp,
)


class Box(BaseGeometry):
    """Box/cuboid geometry with customizable dimensions.

    Boxes are defined by a center position and size dimensions [width, height, depth].
    Optional rotation can be applied as Euler angles [x, y, z] in radians.

    By default, creates a single 1x1x1 box at the origin (0,0,0).
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

    size: Size3DProp = Field(
        default=[(1.0, 1.0, 1.0)],
        description="Box dimensions [(width, height, depth), ...]. String for dynamic.",
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
        description="Box opacity, between 0 (transparent) and 1 (opaque).",
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
