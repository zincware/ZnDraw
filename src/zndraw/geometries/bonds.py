from pydantic import Field
import typing as t

from .base import BaseGeometry, DataProp, InteractionSettings




class Bond(BaseGeometry):
    """A bond geometry."""

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        schema["properties"]["position"]["x-dynamic-enum"] = "AVAILABLE_ATOMS_KEYS"
        schema["properties"]["position"]["type"] = "string"
        schema["properties"]["position"].pop("anyOf", None)

        schema["properties"]["connectivity"]["x-dynamic-enum"] = "AVAILABLE_ATOMS_KEYS"
        schema["properties"]["connectivity"]["type"] = "string"
        schema["properties"]["connectivity"].pop("anyOf", None)

        schema["properties"]["color"]["x-dynamic-enum"] = "AVAILABLE_ATOMS_KEYS"
        schema["properties"]["color"]["type"] = "string"
        schema["properties"]["color"]["x-color-picker"] = True
        schema["properties"]["color"].pop("anyOf", None)

        schema["properties"]["radius"]["x-dynamic-enum"] = "AVAILABLE_ATOMS_KEYS"
        schema["properties"]["radius"]["type"] = "string"
        schema["properties"]["radius"].pop("anyOf", None)
        
        schema["$defs"]["InteractionSettings"]["properties"]["color"][
            "x-color-picker"
        ] = True
        schema["$defs"]["InteractionSettings"]["properties"]["color"][
            "x-dynamic-enum"
        ] = "AVAILABLE_ATOMS_KEYS"
        return schema

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

