"""Box geometry for ZnDraw."""

import typing as t
from pydantic import Field

from .base import BaseGeometry, RotationProp, InteractionSettings, apply_schema_feature


class Box(BaseGeometry):
    """Box/cuboid geometry with customizable dimensions.

    Boxes are defined by a center position and size dimensions [width, height, depth].
    Optional rotation can be applied as Euler angles [x, y, z] in radians.

    By default, creates a single 1x1x1 box at the origin (0,0,0).
    """

    # Override defaults for user-created geometries
    position: t.Union[str, list[tuple[float, float, float]]] = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Position coordinates. String for dynamic data key (e.g. 'arrays.positions'), list of tuples for static per-instance positions [(x,y,z), ...].",
    )

    color: t.Union[str, list[str]] = Field(
        default="#808080",
        description="Color values. String for dynamic key (e.g. 'arrays.colors') or shared hex color (e.g. '#FF0000'), list of hex colors for per-instance ['#FF0000', '#00FF00', ...].",
    )

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Apply schema features - don't force string type for tuple/array fields
        apply_schema_feature(schema, "position", ["dynamic-atom-props"])
        apply_schema_feature(schema, "size", ["dynamic-atom-props"])
        apply_schema_feature(schema, "color", ["color-picker", "dynamic-atom-props", "free-solo"])
        apply_schema_feature(schema, "rotation", ["dynamic-atom-props"])
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo"],
            definition_path="InteractionSettings"
        )

        return schema

    size: t.Union[str, tuple[float, float, float], list[tuple[float, float, float]]] = Field(
        default=(1.0, 1.0, 1.0),
        description="Box dimensions [width, height, depth]. Tuple for shared size across all instances, list of tuples for per-instance sizes, string for dynamic data key.",
    )

    rotation: RotationProp = Field(
        default=(0.0, 0.0, 0.0),
        description="Rotation as Euler angles [x, y, z] in radians. Tuple for shared rotation across all instances, list of tuples for per-instance rotations, string for dynamic data key.",
    )

    scale: float = Field(
        default=1.0,
        ge=0.0,
        description="Uniform scale factor applied to box size.",
    )

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Box opacity, between 0 (transparent) and 1 (opaque).",
    )

    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings.",
    )

    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings.",
    )
