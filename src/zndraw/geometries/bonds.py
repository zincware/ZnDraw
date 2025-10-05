from pydantic import Field

from .base import BaseGeometry, DataProp, InteractionSettings




class Bond(BaseGeometry):
    """A bond geometry."""

    radius: DataProp = Field(
        default="arrays.radii",
        description="Bond radius. String for dynamic data key, float for static value.",
    )

    resolution: int = Field(
        default=16,
        ge=4,
        le=64,
        description="Bond geometry resolution (number of segments). Higher values = smoother bond.",
    )

    connectivity: DataProp = Field(
        default="info.connectivity",
        description="Connectivity information. String for dynamic data key, list of tuples for static value.",
    )

    scale: float = Field(
        default=0.25,
        ge=0.0,
        description="Uniform scale factor applied to bond radius.",
    )
    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Bond opacity, between 0 (transparent) and 1 (opaque).",
    )

    selecting: InteractionSettings = Field(
        default_factory=lambda: InteractionSettings(enabled=True, color="#FF6A00", opacity=0.5),
        description="Selection interaction settings."
    )
    hovering: InteractionSettings = Field(
        default_factory=lambda: InteractionSettings(enabled=True, color="#FF0000", opacity=0.5),
        description="Hover interaction settings."
    )

