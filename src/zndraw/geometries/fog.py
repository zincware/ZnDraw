"""Fog scene object for distance-based fog effect."""

from pydantic import BaseModel, ConfigDict, Field


class Fog(BaseModel):
    """Scene fog effect (distance-based).

    Creates atmospheric depth by fading objects to a color based on distance.
    Use 'default' color to match theme background automatically.
    """

    model_config = ConfigDict(frozen=True)

    owner: str | None = Field(
        default=None,
        description="User ID of the fog owner.",
        json_schema_extra={"x-custom-type": "ownership-toggle"},
    )
    active: bool = Field(
        default=True,
        description="Enable fog effect.",
    )
    color: str = Field(
        default="default",
        description="Fog color ('default' uses theme background).",
        json_schema_extra={"format": "color"},
    )
    near: float = Field(
        default=180.0,
        ge=0.0,
        le=1000.0,
        description="Distance where fog starts.",
        json_schema_extra={"format": "range", "step": 10},
    )
    far: float = Field(
        default=300.0,
        ge=10.0,
        le=2000.0,
        description="Distance where fog is fully opaque.",
        json_schema_extra={"format": "range", "step": 10},
    )
