import typing as t

from pydantic import Field, field_validator

from .base import (
    BaseGeometry,
    InteractionSettings,
    ScaleProp,
    SizeProp,
    apply_schema_feature,
    normalize_scale as _normalize_scale,
)


class Sphere(BaseGeometry):
    """A sphere geometry.

    By default, uses dynamic references to arrays.positions and arrays.colors.
    This makes it suitable for visualizing atoms/particles.
    """

    # Keep defaults from BaseGeometry (dynamic references)
    # position: default="arrays.positions"
    # color: default="arrays.colors"

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Apply schema features using helper
        apply_schema_feature(
            schema, "position", ["dynamic-atom-props", "editable-array", "transform"]
        )
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo", "editable-array"],
        )
        apply_schema_feature(
            schema, "radius", ["dynamic-atom-props", "editable-array", "transform"]
        )
        apply_schema_feature(schema, "scale", ["dynamic-atom-props", "editable-array"])
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo", "editable-array"],
            definition_path="InteractionSettings",
        )

        return schema

    radius: SizeProp = Field(
        default=1,
        description="Sphere radius. String for dynamic data key, float for shared value across all instances, list for per-instance radii.",
    )

    resolution: int = Field(
        default=16,
        ge=4,
        le=64,
        description="Sphere geometry resolution (number of segments). Higher values = smoother sphere.",
    )

    scale: ScaleProp = Field(
        default=1.0,
        description="Scale factor(s). Float for uniform, tuple [sx,sy,sz] for anisotropic, list for per-instance. String for dynamic data key.",
    )
    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Sphere opacity, between 0 (transparent) and 1 (opaque).",
    )

    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings.",
    )
    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings.",
    )

    @field_validator("scale", mode="before")
    @classmethod
    def normalize_scale(cls, v):
        """Normalize and validate scale using shared function."""
        return _normalize_scale(v, validate=True)
