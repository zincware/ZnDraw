"""Shape geometry for ZnDraw - 2D polygons rendered in 3D space."""

from pydantic import Field

from .base import (
    BaseGeometry,
    ColorProp,
    InteractionSettings,
    PositionProp,
    RotationProp,
    ScaleProp,
)


class Shape(BaseGeometry):
    """A 2D polygon shape rendered in 3D space.

    Vertices define the shape outline in 2D (XY plane).
    Position, rotation, scale place instances in 3D space.
    Supports instancing: one shape definition, multiple placements.

    Note: Vertex editing is not supported because vertices are 2D
    but transform controls operate in 3D.

    Example
    -------
    >>> from zndraw import ZnDraw, Shape
    >>> vis = ZnDraw()
    >>> vis.geometries["triangle"] = Shape(
    ...     vertices=[(0, 0), (1, 0), (0.5, 0.866)],
    ...     position=[(0, 0, 0)],
    ...     color=["#FF0000"],
    ... )
    """

    vertices: list[tuple[float, float]] = Field(
        default=[(0.0, 0.0), (1.0, 0.0), (0.5, 0.866)],
        description="Shape vertices in 2D [(x,y), ...]. Connected by lines, auto-closed.",
        json_schema_extra={"x-custom-type": "vertices-2d"},
    )

    position: PositionProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Position coordinates [(x,y,z), ...]. String for dynamic data key.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    rotation: RotationProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Rotation as Euler angles [(x, y, z), ...] in radians. String for dynamic data key.",
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

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Shape opacity.",
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
