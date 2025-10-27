"""Arrow geometry for ZnDraw."""

import typing as t

from pydantic import Field, field_validator

from .base import (
    BaseGeometry,
    InteractionSettings,
    PositionProp,
    apply_schema_feature,
)


class Arrow(BaseGeometry):
    """Arrow geometry with direction vector.

    Arrows are defined by a start position and a direction vector.
    The direction vector determines both the arrow's orientation and its length.
    """

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Apply schema features using helper
        apply_schema_feature(
            schema, "position", ["dynamic-atom-props", "editable-array"]
        )
        apply_schema_feature(
            schema, "direction", ["dynamic-atom-props", "editable-array"]
        )
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo", "editable-array"],
        )
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo", "editable-array"],
            definition_path="InteractionSettings",
        )

        return schema

    # Use 'start' as the primary field with 'position' alias for compatibility
    position: PositionProp = Field(
        default="arrays.positions",
        description="Arrow start position [x,y,z]. String for dynamic data key, list of tuples for static per-instance positions.",
    )

    direction: PositionProp = Field(
        default="calc.forces",
        description="Direction vector [x,y,z]. Defines arrow orientation and base length. String for dynamic data key, list of tuples for static per-instance directions.",
    )

    radius: float = Field(
        default=1,
        description="Arrow shaft radius. String for dynamic data key, float for static value.",
    )

    scale: float = Field(
        default=1.0,
        description="Length scale multiplier applied to direction vector length. String for dynamic data key, float for static value.",
    )

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Arrow opacity, between 0 (transparent) and 1 (opaque).",
    )

    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings.",
    )

    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings.",
    )

    @field_validator("direction", mode="before")
    @classmethod
    def normalize_direction(cls, v):
        """Normalize direction vector to list of tuples (inherited from BaseGeometry for position)."""
        if v is None:
            return []
        if isinstance(v, str):  # Dynamic reference
            return v
        if isinstance(v, tuple):  # Single tuple -> wrap in list
            return [v]
        return v  # Already a list
