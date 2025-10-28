import typing as t

from pydantic import Field

from .base import (
    BaseGeometry,
    ConnectivityProp,
    InteractionSettings,
    apply_schema_feature,
)


class Bond(BaseGeometry):
    """A bond geometry."""

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)

        # Apply schema features using helper
        apply_schema_feature(
            schema, "position", ["dynamic-atom-props", "editable-array"]
        )
        apply_schema_feature(
            schema, "connectivity", ["dynamic-atom-props", "editable-array"]
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

    radius: float = Field(
        default=1,
        description="Bond radius.",
    )

    resolution: int = Field(
        default=16,
        ge=4,
        le=64,
        description="Bond geometry resolution (number of segments). Higher values = smoother bond.",
    )

    connectivity: ConnectivityProp = Field(
        default="info.connectivity",
        description="Connectivity information. String for dynamic data key, list of tuples for static value.",
    )

    scale: float = Field(
        default=0.15,
        ge=0.0,
        description="Uniform scale factor applied to bond radius.",
    )
    opacity: float = Field(
        default=1.0,
        ge=0.0,
        le=1.0,
        description="Bond opacity, between 0 (transparent) and 1 (opaque).",
    )

    bond_order_mode: t.Literal["parallel", "ignore"] = Field(
        default="parallel",
        description="Bond order visualization mode. 'parallel': render multiple cylinders for double/triple bonds. 'ignore': render all bonds as single cylinders.",
    )

    bond_order_offset: float = Field(
        default=3,
        ge=0.0,
        description="Spacing between parallel cylinders in parallel mode, as fraction of bond radius.",
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
