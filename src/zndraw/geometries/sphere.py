from pydantic import Field
import typing as t

from .base import BaseGeometry, SizeProp, InteractionSettings, apply_schema_feature


class Sphere(BaseGeometry):
    """A sphere geometry.

    By default, creates a single sphere at the origin (0,0,0).
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
        apply_schema_feature(schema, "color", ["color-picker", "dynamic-atom-props", "free-solo"])
        apply_schema_feature(schema, "radius", ["dynamic-atom-props"])
        apply_schema_feature(
            schema,
            "color",
            ["color-picker", "dynamic-atom-props", "free-solo"],
            definition_path="InteractionSettings"
        )

        return schema

    radius: SizeProp = Field(
        default="arrays.radii",
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
