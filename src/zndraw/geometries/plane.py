"""Plane geometry for ZnDraw."""

import typing as t

from pydantic import Field, field_validator

from .base import (
    BaseGeometry,
    ColorProp,
    InteractionSettings,
    PositionProp,
    RotationProp,
    Size2DProp,
    apply_schema_feature,
)


class Plane(BaseGeometry):
    """Flat rectangular plane geometry with customizable dimensions.

    Planes are defined by a center position and size dimensions [width, height].
    Optional rotation can be applied as Euler angles [x, y, z] in radians.
    By default, planes are oriented in the XY plane (facing along Z-axis).

    By default, creates a single 1x1 plane at the origin (0,0,0).
    """

    # Override defaults for user-created geometries
    position: PositionProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Position coordinates. String for dynamic data key (e.g. 'arrays.positions'), list of tuples for static per-instance positions [(x,y,z), ...].",
    )

    color: ColorProp = Field(
        default=["#808080"],
        description="Color values. String for dynamic key (e.g. 'arrays.colors') or list of hex colors ['#FF0000', '#00FF00', ...].",
    )

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Apply schema features using helper
        apply_schema_feature(
            schema, "position", ["dynamic-atom-props", "editable-array"]
        )
        apply_schema_feature(schema, "size", ["dynamic-atom-props", "editable-array"])
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo", "editable-array"],
        )
        apply_schema_feature(
            schema, "rotation", ["dynamic-atom-props", "editable-array"]
        )
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo", "editable-array"],
            definition_path="InteractionSettings",
        )

        return schema

    size: Size2DProp = Field(
        default=[(1.0, 1.0)],
        description="Plane dimensions [width, height]. List of tuples for per-instance sizes, string for dynamic data key. Single tuples are automatically normalized to lists.",
    )

    rotation: RotationProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Rotation as Euler angles [x, y, z] in radians. List of tuples for per-instance rotations, string for dynamic data key. Single tuples are automatically normalized to lists.",
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

    @field_validator("position", "rotation", "size", mode="before")
    @classmethod
    def normalize_vector_fields(cls, v):
        """Normalize vector fields to list of tuples."""
        if v is None:
            return []
        if isinstance(v, str):  # Dynamic reference
            return v
        if isinstance(v, tuple):  # Single tuple -> wrap in list
            return [v]
        return v  # Already a list

    @field_validator("color", mode="before")
    @classmethod
    def normalize_color(cls, v):
        """Normalize color to list of hex strings."""
        if v is None:
            return []
        # Dynamic reference -> pass through
        if isinstance(v, str) and not v.startswith("#"):
            return v
        # Single hex color -> wrap in list
        if isinstance(v, str) and v.startswith("#"):
            return [v]
        # Already a list -> return as is
        return v
