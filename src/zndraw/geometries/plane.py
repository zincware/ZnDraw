"""Plane geometry for ZnDraw."""

import typing as t
from pydantic import Field

from .base import BaseGeometry, RotationProp, InteractionSettings, apply_schema_feature

# Custom type for 2D plane size (width, height)
# Can be shared (single tuple) or per-instance (list of tuples) or dynamic
PlaneSize = t.Union[str, tuple[float, float], list[tuple[float, float]]]


class Plane(BaseGeometry):
    """Flat rectangular plane geometry with customizable dimensions.

    Planes are defined by a center position and size dimensions [width, height].
    Optional rotation can be applied as Euler angles [x, y, z] in radians.
    By default, planes are oriented in the XY plane (facing along Z-axis).

    By default, creates a single 1x1 plane at the origin (0,0,0).
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

        # Apply schema features using helper
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

    size: PlaneSize = Field(
        default=(1.0, 1.0),
        description="Plane dimensions [width, height]. Tuple for shared size across all instances, list of tuples for per-instance sizes, string for dynamic data key.",
    )

    rotation: RotationProp = Field(
        default=(0.0, 0.0, 0.0),
        description="Rotation as Euler angles [x, y, z] in radians. Tuple for shared rotation across all instances, list of tuples for per-instance rotations, string for dynamic data key.",
    )

    scale: float = Field(
        default=1.0,
        ge=0.0,
        description="Uniform scale factor applied to plane size.",
    )

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Plane opacity, between 0 (transparent) and 1 (opaque).",
    )

    double_sided: bool = Field(
        default=True,
        description="If true, plane is visible from both sides. If false, only visible from front.",
    )

    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings.",
    )

    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings.",
    )
