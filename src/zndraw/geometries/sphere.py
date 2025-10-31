import typing as t

from pydantic import Field

from .base import BaseGeometry, InteractionSettings, SizeProp, apply_schema_feature


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

    scale: float = Field(
        default=0.7,
        ge=0.0,
        description="Uniform scale factor applied to sphere radius.",
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
