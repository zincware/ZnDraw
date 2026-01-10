"""Bond geometry for ZnDraw."""

import typing as t

from pydantic import Field

from .base import (
    BaseGeometry,
    ConnectivityProp,
    InteractionSettings,
)


class Bond(BaseGeometry):
    """A bond geometry."""

    position: str = Field(
        default="arrays.positions",
        description="Position data key for atom positions.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    color: str = Field(
        default="arrays.colors",
        description="Color data key for bond colors.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": [
                "color-picker",
                "dynamic-atom-props",
                "free-solo",
                "editable-array",
            ],
        },
    )

    connectivity: ConnectivityProp = Field(
        default="info.connectivity",
        description="Connectivity information. String for dynamic data key, list of tuples for static value.",
        json_schema_extra={
            "x-custom-type": "dynamic-enum",
            "x-features": ["dynamic-atom-props", "editable-array"],
        },
    )

    radius: float = Field(
        default=1,
        description="Bond radius.",
        ge=0.1,
        le=10,
        json_schema_extra={"format": "range", "step": 0.1},
    )

    resolution: int = Field(
        default=16,
        ge=4,
        le=64,
        description="Bond geometry resolution (number of segments). Higher values = smoother bond.",
        json_schema_extra={"format": "range", "step": 1},
    )

    scale: float = Field(
        default=1.0,
        ge=0.0,
        le=10,
        description="Uniform scale factor applied to bond radius.",
        json_schema_extra={"format": "range", "step": 0.1},
    )

    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Bond opacity, between 0 (transparent) and 1 (opaque).",
        json_schema_extra={"format": "range", "step": 0.01},
    )

    bond_order_mode: t.Literal["parallel", "ignore"] = Field(
        default="parallel",
        description="Bond order visualization mode. 'parallel': render multiple cylinders for double/triple bonds. 'ignore': render all bonds as single cylinders.",
    )

    bond_order_offset: float = Field(
        default=3,
        ge=0.0,
        le=10,
        description="Spacing between parallel cylinders in parallel mode, as fraction of bond radius.",
        json_schema_extra={"format": "range", "step": 0.1},
    )

    bond_order_radius_scale: dict[float, float] = Field(
        default={1: 1.0, 1.5: 0.75, 2: 0.75, 3: 0.7},
        description="Radius multiplier per bond order in parallel mode. Maps bond order to scale factor.",
    )

    selecting: InteractionSettings = Field(
        default_factory=lambda: InteractionSettings(
            enabled=True, color="#FF6A00", opacity=0.5
        ),
        description="Selection interaction settings.",
    )

    hovering: InteractionSettings = Field(
        default_factory=lambda: InteractionSettings(
            enabled=True, color="#FF0000", opacity=0.5
        ),
        description="Hover interaction settings.",
    )
