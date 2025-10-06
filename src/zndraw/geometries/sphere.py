from pydantic import Field

from .base import BaseGeometry, DataProp, InteractionSettings




class Sphere(BaseGeometry):
    """A sphere geometry."""

    radius: DataProp = Field(
        default="arrays.radii",
        description="Sphere radius. String for dynamic data key, float for static value.",
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

    selecting: InteractionSettings | None = Field(
        default_factory=lambda: InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings."
    )
    hovering: InteractionSettings | None = Field(
        default_factory=lambda: InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings."
    )
