"""Arrow geometry for ZnDraw."""

from pydantic import Field

from .base import (
    BaseGeometry,
    InteractionSettings,
    PositionProp,
    ScaleProp,
)


class Arrow(BaseGeometry):
    """Arrow geometry with direction vector.

    Arrows are defined by a start position and a direction vector.
    The direction vector determines both the arrow's orientation and its length.
    """

    position: PositionProp = Field(
        default="arrays.positions",
        description="Arrow start positions [(x,y,z), ...]. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    direction: PositionProp = Field(
        default="calc.forces",
        description="Direction vectors [(x,y,z), ...]. Defines arrow orientation and length. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    radius: float = Field(
        default=1.0,
        description="Arrow shaft radius.",
        ge=0.1,
        le=10,
        json_schema_extra={"format": "range", "step": 0.1},
    )

    scale: ScaleProp = Field(
        default=[(1.0, 1.0, 1.0)],
        description="Scale factors [(sx, sy, sz), ...]. First component used as length multiplier.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    resolution: int = Field(
        default=16,
        ge=4,
        le=64,
        description="Arrow geometry resolution (number of segments).",
        json_schema_extra={"format": "range", "step": 1},
    )

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Arrow opacity, between 0 (transparent) and 1 (opaque).",
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
