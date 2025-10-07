from pydantic import Field
import typing as t

from .base import BaseGeometry, DataProp, InteractionSettings




class Sphere(BaseGeometry):
    """A sphere geometry."""

    @classmethod
    def model_json_schema(cls, **kwargs: t.Any) -> dict[str, t.Any]:
        schema = super().model_json_schema(**kwargs)
        schema["properties"]["position"]["x-dynamic-enum"] = "AVAILABLE_ATOMS_KEYS"
        schema["properties"]["position"]["type"] = "string"
        schema["properties"]["position"].pop("anyOf", None)
        schema["properties"]["color"]["x-dynamic-enum"] = "AVAILABLE_ATOMS_KEYS"
        schema["properties"]["color"]["type"] = "string"
        schema["properties"]["color"]["x-color-picker"] = True
        schema["properties"]["color"].pop("anyOf", None)
        schema["properties"]["radius"]["x-dynamic-enum"] = "AVAILABLE_ATOMS_KEYS"
        schema["properties"]["radius"]["type"] = "string"
        schema["properties"]["radius"].pop("anyOf", None)
        return schema

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

    selecting: InteractionSettings = Field(
        default=InteractionSettings(color="#FF6A00", opacity=0.5),
        description="Selection interaction settings."
    )
    hovering: InteractionSettings = Field(
        default=InteractionSettings(color="#FF0000", opacity=0.5),
        description="Hover interaction settings."
    )
