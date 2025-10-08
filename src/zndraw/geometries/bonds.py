from pydantic import Field
import typing as t

from .base import BaseGeometry, DataProp, InteractionSettings




class Bond(BaseGeometry):
    """A bond geometry."""

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        schema["properties"]["position"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["position"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["position"]["type"] = "string"
        schema["properties"]["position"].pop("anyOf", None)

        schema["properties"]["connectivity"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["connectivity"]["x-features"] = ["dynamic-atom-props"]
        schema["properties"]["connectivity"]["type"] = "string"
        schema["properties"]["connectivity"].pop("anyOf", None)

        schema["properties"]["color"]["x-custom-type"] = "dynamic-enum"
        schema["properties"]["color"]["type"] = "string"
        schema["properties"]["color"]["x-features"] = ["color-picker", "dynamic-atom-props", "free-solo"]
        schema["properties"]["color"].pop("anyOf", None)

        schema["$defs"]["InteractionSettings"]["properties"]["color"][
            "x-custom-type"
        ] = "dynamic-enum"
        schema["$defs"]["InteractionSettings"]["properties"]["color"][
            "x-features"
        ] = ["color-picker", "dynamic-atom-props", "free-solo"]
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

