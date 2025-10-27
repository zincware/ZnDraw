"""Box geometry for ZnDraw."""

import typing as t

from pydantic import Field, field_validator

from .base import (
    BaseGeometry,
    ColorProp,
    InteractionSettings,
    PositionProp,
    RotationProp,
    Size3DProp,
    apply_schema_feature,
)


class Box(BaseGeometry):
    """Box/cuboid geometry with customizable dimensions.

    Boxes are defined by a center position and size dimensions [width, height, depth].
    Optional rotation can be applied as Euler angles [x, y, z] in radians.

    By default, creates a single 1x1x1 box at the origin (0,0,0).
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

        # Apply schema features - don't force string type for tuple/array fields
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

    size: Size3DProp = Field(
        default=[(1.0, 1.0, 1.0)],
        description="Box dimensions [width, height, depth]. List of tuples for per-instance sizes, string for dynamic data key. Single tuples are automatically normalized to lists.",
    )

    rotation: RotationProp = Field(
        default=[(0.0, 0.0, 0.0)],
        description="Rotation as Euler angles [x, y, z] in radians. List of tuples for per-instance rotations, string for dynamic data key. Single tuples are automatically normalized to lists.",
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

    @field_validator("rotation", "size", mode="before")
    @classmethod
    def normalize_vector_fields(cls, v):
        """Normalize vector fields to list of tuples (position inherited from BaseGeometry)."""
        if v is None:
            return []
        if isinstance(v, str):  # Dynamic reference
            return v
        if isinstance(v, tuple):  # Single tuple -> wrap in list
            return [v]
        return v  # Already a list
